#include "image_lib.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>

using namespace std;

static void skipWhitespaceAndComments(istream& in)
{
    while (in.peek() != EOF) {
        char c = static_cast<char>(in.peek());
        if (c == '#') {
            string line;
            getline(in, line);
        } else if (isspace(static_cast<unsigned char>(c))) {
            in.get();
        } else {
            break;
        }
    }
}

Image readPPM(const string& filename)
{
    ifstream file(filename);
    if (!file.is_open()){
        throw runtime_error("Cannot open input file: " + filename);
    }
 
    skipWhitespaceAndComments(file);
    string magic;
    file >> magic;
    if (magic != "P3"){
        throw runtime_error("Unsupported format: expected P3, got " + magic);
    }
 
    Image img;
    skipWhitespaceAndComments(file);
    file >> img.width;
    skipWhitespaceAndComments(file);
    file >> img.height;
    skipWhitespaceAndComments(file);
    file >> img.maxVal;
 
    if (img.width <= 0 || img.height <= 0 || img.maxVal <= 0){
        throw runtime_error("Invalid PPM header values.");
    }
 
    img.pixels.assign(img.height,
        vector<vector<int>>(img.width,
            vector<int>(3, 0)));
 
    for (int row = 0; row < img.height; ++row) {
        for (int col = 0; col < img.width; ++col) {
            for (int ch = 0; ch < 3; ++ch) {
                if (!(file >> img.pixels[row][col][ch])){
                    throw runtime_error("Unexpected end of pixel data at row " + to_string(row) + ", col " + to_string(col));
                }
            }
        }
    }
 
    cout << "[read]  " << filename
              << "  (" << img.width << "x" << img.height
              << ", max=" << img.maxVal << ")\n";
    return img;
}

void writePPM(const string& filename, const Image& img)
{
    ofstream file(filename);
    if (!file.is_open()){
        throw runtime_error("Cannot open output file: " + filename);
    }
 
    int outMaxVal = img.maxVal;
    int multiplier = 1;
    
    if (img.maxVal == 1) {
        outMaxVal = 255;
        multiplier = 255;
    }

    file << "P3\n"
         << img.width << ' ' << img.height << '\n'
         << outMaxVal << '\n';
 
    const int tripletsPerLine = 5;
    int count = 0;
    for (int row = 0; row < img.height; ++row) {
        for (int col = 0; col < img.width; ++col) {
            file << (img.pixels[row][col][0] * multiplier) << ' '
                 << (img.pixels[row][col][1] * multiplier) << ' '
                 << (img.pixels[row][col][2] * multiplier);
            ++count;
            if (count % tripletsPerLine == 0)
                file << '\n';
            else
                file << "  ";
        }
    }
    if (count % tripletsPerLine != 0){
        file << '\n';
    }
 
    cout << "[write] " << filename << '\n';
}

Image applyGaussianBlur(const Image& src)
{
    static const int KERNEL[3][3] = {
        {1, 2, 1},
        {2, 4, 2},
        {1, 2, 1}
    };
    static const int KERNEL_SUM = 16;

    Image dst;
    dst.width  = src.width;
    dst.height = src.height;
    dst.maxVal = src.maxVal;
 
    dst.pixels.assign(dst.height,
        vector<vector<int>>(dst.width,
            vector<int>(3, 0)));
 
    for (int row = 0; row < src.height; ++row) {
        for (int col = 0; col < src.width; ++col) {
            for (int ch = 0; ch < 3; ++ch) {
 
                int acc = 0;
                for (int ki = -1; ki <= 1; ++ki) {
                    for (int kj = -1; kj <= 1; ++kj) {
                        int weight = KERNEL[ki + 1][kj + 1];
                        acc += weight * clampedPixel(src, row + ki, col + kj, ch);
                    }
                }
 
                int value = acc / KERNEL_SUM;
                dst.pixels[row][col][ch] = max(0, min(value, src.maxVal));
            }
        }
    }
 
    return dst;
}

std::pair<Image, Image> applySobelFilter(const Image& src)
{
    static const int KERNEL_X[3][3] = {
        {-1, 0, 1},
        {-2, 0, 2},
        {-1, 0, 1}
    };
    static const int KERNEL_Y[3][3] = {
        {-1, -2, -1},
        {0, 0, 0},
        {1, 2, 1}
    };
 
    Image dst_h, dst_v;
    dst_h.width = dst_v.width = src.width;
    dst_h.height = dst_v.height = src.height;
    dst_h.maxVal = dst_v.maxVal = src.maxVal;
 
    dst_h.pixels.assign(dst_h.height, vector<vector<int>>(dst_h.width, vector<int>(3, 0)));
    dst_v.pixels.assign(dst_v.height, vector<vector<int>>(dst_v.width, vector<int>(3, 0)));
 
    for (int row = 0; row < src.height; ++row) {
        for (int col = 0; col < src.width; ++col) {
            for (int ch = 0; ch < 3; ++ch) {
                int acc_x = 0;
                int acc_y = 0;
                for (int ki = -1; ki <= 1; ++ki) {
                    for (int kj = -1; kj <= 1; ++kj) {
                        int p = clampedPixel(src, row + ki, col + kj, ch);
                        acc_x += KERNEL_X[ki + 1][kj + 1] * p;
                        acc_y += KERNEL_Y[ki + 1][kj + 1] * p;
                    }
                }
                dst_h.pixels[row][col][ch] = acc_x;
                dst_v.pixels[row][col][ch] = acc_y;
            }
        }
    }
 
    return std::make_pair(dst_h, dst_v);
}

Image applyLaplaceFilter(const Image& src)
{
    static const int KERNEL[3][3] = {
        {0, -1, 0},
        {-1, 4, -1},
        {0, -1, 0}
    };

    Image dst;
    dst.width  = src.width;
    dst.height = src.height;
    dst.maxVal = src.maxVal;

    dst.pixels.assign(dst.height,
        vector<vector<int>>(dst.width,
            vector<int>(3, 0)));
 
    for (int row = 0; row < src.height; ++row) {
        for (int col = 0; col < src.width; ++col) {
            for (int ch = 0; ch < 3; ++ch) {
                int acc = 0;
                for (int ki = -1; ki <= 1; ++ki) {
                    for (int kj = -1; kj <= 1; ++kj) {
                        int weight = KERNEL[ki + 1][kj + 1];
                        acc += weight * clampedPixel(src, row + ki, col + kj, ch);
                    }
                }
                dst.pixels[row][col][ch] = acc;
            }
        }
    }

    return dst;
}