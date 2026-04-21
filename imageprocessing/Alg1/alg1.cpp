#include "../src/image_lib.h"
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

int main() {
  const string inputPpm = "input.ppm";
  const string outputPpm = "output.ppm";

  try {
    Image img = readPPM(inputPpm);

    // 1) Suavize a imagem, aplicando um filtro Gaussiano
    Image smoothed = applyGaussianBlur(img);

    // 2) Na imagem do passo 1, aplique o filtro convolucional de Gradiente
    // 2.1) um filtro de Sobel para a derivada na direção x
    // 2.2) um filtro de Sobel para a derivada na direção y
    auto [sobelX, sobelY] = applySobelFilter(smoothed);

    Image gradient;
    gradient.width = img.width;
    gradient.height = img.height;
    gradient.maxVal = img.maxVal;
    gradient.pixels.assign(
        gradient.height,
        vector<vector<int>>(gradient.width, vector<int>(3, 0)));

    // 2.3) em cada uma das matrizes, A e B, eleve ao quadrado os valores dos
    // elementos 2.4) some as duas matrizes A e B modificadas no passo 2.3 e
    // tire a raiz quadrada de cada elemento dessa matriz, C
    for (int row = 0; row < gradient.height; ++row) {
      for (int col = 0; col < gradient.width; ++col) {
        for (int ch = 0; ch < 3; ++ch) {
          double valX = sobelX.pixels[row][col][ch];
          double valY = sobelY.pixels[row][col][ch];
          double magnitude = sqrt(valX * valX + valY * valY);
          gradient.pixels[row][col][ch] = static_cast<int>(magnitude);
        }
      }
    }

    // 3)
    float threshold = 200.0f;

    // 4) gere uma matriz Final, D
    Image finalImg;
    finalImg.width = img.width;
    finalImg.height = img.height;
    finalImg.maxVal = 1;
    finalImg.pixels.assign(
        finalImg.height,
        vector<vector<int>>(finalImg.width, vector<int>(3, 0)));

    for (int row = 0; row < finalImg.height; ++row) {
      for (int col = 0; col < finalImg.width; ++col) {
        bool isEdge = false;
        for (int ch = 0; ch < 3; ++ch) {
          if (gradient.pixels[row][col][ch] > threshold) {
            isEdge = true;
            break;
          }
        }
        if (isEdge) {
          finalImg.pixels[row][col][0] = 1;
          finalImg.pixels[row][col][1] = 1;
          finalImg.pixels[row][col][2] = 1;
        }
      }
    }

    writePPM(outputPpm, finalImg);

  } catch (const exception &ex) {
    cerr << ex.what() << '\n';
    return 1;
  }

  return 0;
}
