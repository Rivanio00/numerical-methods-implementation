#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <tuple>
#include <cmath>

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
 
tuple<Image, Image> applySobelFilter(const Image& src)
{
    Image dst_h;
    Image dst_v;

    dst_h.width  = src.width;
    dst_h.height = src.height;
    dst_h.maxVal = src.maxVal;

    dst_v.width  = src.width;
    dst_v.height = src.height;
    dst_v.maxVal = src.maxVal;
 
    dst_h.pixels.assign(dst_h.height,
        vector<vector<int>>(dst_h.width,
            vector<int>(3, 0)));

    dst_v.pixels.assign(dst_v.height,
        vector<vector<int>>(dst_v.width,
            vector<int>(3, 0)));
 
    for (int row = 0; row < src.height; ++row) {
        for (int col = 0; col < src.width; ++col) {
            for (int ch = 0; ch < 3; ++ch) {
 
                int acc = 0;
                // Convolve with 3×3 kernel
                for (int ki = -1; ki <= 1; ++ki) {
                    for (int kj = -1; kj <= 1; ++kj) {
                        int weight = KERNEL_X[ki + 1][kj + 1];
                        acc += weight * clampedPixel(src, row + ki, col + kj, ch);
                    }
                }
 
                dst_h.pixels[row][col][ch] = acc;
            }
        }
    }

    for (int row = 0; row < src.height; ++row) {
        for (int col = 0; col < src.width; ++col) {
            for (int ch = 0; ch < 3; ++ch) {
 
                int acc = 0;
                // Convolve with 3×3 kernel
                for (int ki = -1; ki <= 1; ++ki) {
                    for (int kj = -1; kj <= 1; ++kj) {
                        int weight = KERNEL_Y[ki + 1][kj + 1];
                        acc += weight * clampedPixel(src, row + ki, col + kj, ch);
                    }
                }
 
                dst_v.pixels[row][col][ch] = acc;
            }
        }
    }

    return make_tuple(dst_h, dst_v);
}

void paintPixel(int x, int y, Image& dst)
{
    for (int i = 0; i < dst.pixels[x][y].size(); i++)
    {
        dst.pixels[x][y][i] = dst.maxVal;
    }
}

Image makeSobelMatrix(const Image& src)
{
    static const int THRESHOLD = 200;
    Image dst_h, dst_v;
    tie(dst_h, dst_v) = applySobelFilter(src);

    Image dst_h_square, dst_v_square, dst_square, dst_final;

    dst_square.width  = src.width;
    dst_square.height = src.height;
    dst_square.maxVal = src.maxVal;

    dst_final.width  = src.width;
    dst_final.height = src.height;
    dst_final.maxVal = src.maxVal;

    dst_final.pixels.assign(dst_final.height,
        vector<vector<int>>(dst_final.width,
            vector<int>(3, 0)));

    dst_square.pixels.assign(dst_square.height,
        vector<vector<int>>(dst_square.width,
            vector<int>(3, 0)));

    for (int i = 0; i < dst_h.pixels.size(); i++) {
        for (int j = 0; j < dst_h.pixels[i].size(); j++) {
            for (int k = 0; k < dst_h.pixels[i][j].size(); k++) {
                dst_h.pixels[i][j][k] *= dst_h.pixels[i][j][k]; // Multiplica o valor por ele mesmo
                dst_v.pixels[i][j][k] *= dst_v.pixels[i][j][k];
                dst_square.pixels[i][j][k] = sqrt(dst_h.pixels[i][j][k] + dst_v.pixels[i][j][k]);
            }
        }
    }

    for (int row = 0; row < src.height; ++row) {
        for (int col = 0; col < src.width; ++col) {
            for (int ch = 0; ch < 3; ++ch) {
 
                // Normalise and clamp to valid range

                int value = dst_square.pixels[row][col][ch];

                if (value > THRESHOLD)
                {
                    paintPixel(row, col, dst_final);
                }
            }
        }
    }
    return dst_final;
}

int main()
{
    const string inputFile  = "input.ppm";
    const string outputFile = "output.ppm";
 
    try {
        Image src = readPPM(inputFile);
        Image sobel_matrix;
        sobel_matrix = makeSobelMatrix(src);
        writePPM(outputFile, sobel_matrix);
        cout << "[done]  Sobel filter applied successfully.\n";
    }
    catch (const exception& ex) {
        cerr << "[error] " << ex.what() << '\n';
        return 1;
    }
 
    return 0;
}