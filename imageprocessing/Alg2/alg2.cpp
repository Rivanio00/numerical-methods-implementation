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

    // 2) Na imagem do passo 1, aplique o filtro convolucional de Laplace,
    // gerando uma imagem/matriz A
    Image laplaceMat = applyLaplaceFilter(smoothed);

    // 3) Por simplicidade, gere uma imagem/matriz, B, percorrendo a imagem A e
    // escrevendo em B
    Image finalImg;
    finalImg.width = laplaceMat.width;
    finalImg.height = laplaceMat.height;
    finalImg.maxVal = 1;
    finalImg.pixels.assign(
        finalImg.height,
        vector<vector<int>>(finalImg.width, vector<int>(3, 0)));

    const double TOLERANCE = 40.0f;

    for (int row = 0; row < finalImg.height; ++row) {
      for (int col = 0; col < finalImg.width; ++col) {
        bool isEdge = false;
        for (int ch = 0; ch < 3; ++ch) {
          if (std::abs(laplaceMat.pixels[row][col][ch]) > TOLERANCE) {
            isEdge = true;
            break;
          }
        }

        if (isEdge) {
          // 3.1) caso o pixel correspondente da matriz A seja diferente de 0
          // dentro de uma tolerância
          finalImg.pixels[row][col][0] = 1;
          finalImg.pixels[row][col][1] = 1;
          finalImg.pixels[row][col][2] = 1;
        }
        // 3.0) cado o pixel correspondente da matriz A seja igual a 0 dentro de
        // uma tolerância
      }
    }

    writePPM(outputPpm, finalImg);

  } catch (const exception &ex) {
    cerr << ex.what() << '\n';
    return 1;
  }

  return 0;
}
