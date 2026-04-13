#include <fstream>
#include <iostream>

#include "matrixCalc.h"

int main(int argc, char *argv[]) {
  if (argc != 4) {
    std::cerr << "Wrong argument!" << std::endl;
    return 1;
  }

  std::ifstream in;
  std::ofstream out;
  in.open(argv[2]);
  out.open(argv[3]);

  if (!in) {
    std::cerr << "Cannot open input file" << std::endl;
    return 1;
  }
  if (!out) {
    std::cerr << "Cannot open output file" << std::endl;
    return 1;
  }

  int rows = 0, cols = 0;
  in >> rows >> cols;

  if (rows != cols) {
    out << "error" << std::endl;
    return 1;
  }

  double **matrix;
  create2DArray(matrix, rows);
  readMatrixFromFile(matrix, in, rows);

  if (comp(argv[1], "-det")) {
    double det = detMatrixGaus(matrix, rows);
    out << det << std::endl;

    for (int i = 0; i < rows; ++i) {
      delete[] matrix[i];
    }
    delete[] matrix;
    return 0;
  }

  if (comp(argv[1], "-inv")) {
    double **extended_matrix;
    create2DExtendedArray(matrix, extended_matrix, rows);

    if (!invMatrixGausJordan(extended_matrix, rows)) {
      out << "error" << std::endl;
      for (int i = 0; i < rows; ++i) {
        delete[] extended_matrix[i];
      }
      delete[] extended_matrix;
      for (int i = 0; i < rows; ++i) {
        delete[] matrix[i];
      }
      delete[] matrix;
      return 1;
    }

    double **inverse = getInverseFromExtended(extended_matrix, rows);

    for (int i = 0; i < rows; ++i) {
      delete[] extended_matrix[i];
    }
    delete[] extended_matrix;

    printMatrixToFile(inverse, rows, out);

    for (int i = 0; i < rows; ++i) {
      delete[] matrix[i];
      delete[] inverse[i];
    }
    delete[] matrix;
    delete[] inverse;

    return 0;
  }

  out << "error" << std::endl;
  return 1;
}