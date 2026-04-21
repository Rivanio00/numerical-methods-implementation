#pragma once

#include <vector>
#include <string>
#include <algorithm>

struct Image {
    int width  = 0;
    int height = 0;
    int maxVal = 255;
    std::vector<std::vector<std::vector<int>>> pixels;
};

Image readPPM(const std::string& filename);
void  writePPM(const std::string& filename, const Image& img);

inline int clampedPixel(const Image& img, int row, int col, int ch)
{
    row = std::max(0, std::min(row, img.height - 1));
    col = std::max(0, std::min(col, img.width  - 1));
    return img.pixels[row][col][ch];
}

Image applyGaussianBlur(const Image& src);

std::pair<Image, Image> applySobelFilter(const Image& src);

Image applyLaplaceFilter(const Image& src);