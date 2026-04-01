#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <sstream>

using namespace std;

struct Image {
    int width  = 0;
    int height = 0;
    int maxVal = 255;
    vector<vector<vector<int>>> pixels;
};

static void skipWhitespaceAndComments(istream& in)
{
    while (in.peek() != EOF) {
        char c = static_cast<char>(in.peek());
        if (c == '#') {
            string line;
            getline(in, line);   // discard comment
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
 
    file << "P3\n"
         << img.width << ' ' << img.height << '\n'
         << img.maxVal << '\n';
 
    const int tripletsPerLine = 5;
    int count = 0;
    for (int row = 0; row < img.height; ++row) {
        for (int col = 0; col < img.width; ++col) {
            file << img.pixels[row][col][0] << ' '
                 << img.pixels[row][col][1] << ' '
                 << img.pixels[row][col][2];
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

inline int clampedPixel(const Image& img, int row, int col, int ch)
{
    row = max(0, min(row, img.height - 1));
    col = max(0, min(col, img.width  - 1));
    return img.pixels[row][col][ch];
}

static const int KERNEL[3][3] = {
    {1, 2, 1},
    {2, 4, 2},
    {1, 2, 1}
};
static const int KERNEL_SUM = 16;
 
Image applyGaussianBlur(const Image& src)
{
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
                // Convolve with 3×3 kernel
                for (int ki = -1; ki <= 1; ++ki) {
                    for (int kj = -1; kj <= 1; ++kj) {
                        int weight = KERNEL[ki + 1][kj + 1];
                        acc += weight * clampedPixel(src, row + ki, col + kj, ch);
                    }
                }
 
                // Normalise and clamp to valid range
                int value = acc / KERNEL_SUM;
                dst.pixels[row][col][ch] =
                    max(0, min(value, src.maxVal));
            }
        }
    }
 
    return dst;
}

int main()
{
    const string inputFile  = "input.ppm";
    const string outputFile = "output.ppm";
 
    try {
        Image src = readPPM(inputFile);
        Image dst = applyGaussianBlur(src);
        writePPM(outputFile, dst);
        cout << "[done]  Gaussian blur applied successfully.\n";
    }
    catch (const exception& ex) {
        cerr << "[error] " << ex.what() << '\n';
        return 1;
    }
 
    return 0;
}