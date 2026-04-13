#include "matrixCalc.h"

void create2DArray(double **&arr, int n) {
  arr = new double *[n];
  for (int i = 0; i < n; ++i) {
    arr[i] = new double[n];
  }
}

void readMatrixFromFile(double **&arr, std::ifstream &in, int n) {
  double tmp = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      in >> tmp;
      arr[i][j] = tmp;
    }
  }
}

void printMatrixToConsole(double **&arr, int n) {
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      std::cout << arr[i][j] << " ";
    }
    std::cout << std::endl;
  }
}

double detMatrixGaus(double **arr, int n) {
  double det = 1;
  for (int j = 0; j < n; ++j) {
    int pivot_row = j;
    double max_abs = myAbs(arr[j][j]);

    for (int i = j + 1; i < n; ++i) {
      double current_abs = myAbs(arr[i][j]);
      if (current_abs > max_abs) {
        pivot_row = i;
        max_abs = current_abs;
      }
    }

    if (max_abs < 1e-12) {
      return 0;
    }

    if (pivot_row != j) {
      std::swap(arr[j], arr[pivot_row]);
      det = -det;
    }

    for (int i = j + 1; i < n; ++i) {
      double factor = arr[i][j] / arr[j][j];
      for (int k = j + 1; k < n; ++k) {
        arr[i][k] -= factor * arr[j][k];
      }
    }

    det *= arr[j][j];
  }
  return det;
}

double myAbs(double n) {
  if (n < 0) {
    return -n;
  }
  return n;
}

bool comp(char *arg, const char *tmp) {
  int i = 0;

  while (arg[i] != '\0' && tmp[i] != '\0') {
    if (arg[i] != tmp[i]) {
      return false;
    }
    i++;
  }

  return (arg[i] == '\0' && tmp[i] == '\0');
}

void create2DExtendedArray(double **&arr, double **&extended_arr, int n) {
  extended_arr = new double *[n];
  for (int i = 0; i < n; ++i) {
    extended_arr[i] = new double[2 * n];
  }

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      extended_arr[i][j] = arr[i][j];
    }

    for (int j = 0; j < n; ++j) {
      if (i == j) {
        extended_arr[i][n + j] = 1.0;
      } else {
        extended_arr[i][n + j] = 0.0;
      }
    }
  }
}

bool invMatrixGausJordan(double **&extended_arr, int n) {
  for (int col = 0; col < n; ++col) {
    int pivot_row = col;
    double max_abs = myAbs(extended_arr[col][col]);

    for (int row = col + 1; row < n; ++row) {
      double current_abs = myAbs(extended_arr[row][col]);
      if (current_abs > max_abs) {
        max_abs = current_abs;
        pivot_row = row;
      }
    }

    if (max_abs < 1e-10) {
      return false;
    }

    if (pivot_row != col) {
      for (int j = 0; j < 2 * n; ++j) {
        std::swap(extended_arr[col][j], extended_arr[pivot_row][j]);
      }
    }

    double divisor = extended_arr[col][col];
    for (int j = 0; j < 2 * n; ++j) {
      extended_arr[col][j] /= divisor;
    }

    for (int row = 0; row < n; ++row) {
      if (row == col) {
        continue;
      }

      double factor = extended_arr[row][col];
      if (factor == 0) {
        continue;
      }

      for (int j = 0; j < 2 * n; ++j) {
        extended_arr[row][j] -= factor * extended_arr[col][j];
      }
    }
  }
  return true;
}

double **getInverseFromExtended(double **extended_arr, int n) {
  double **inverse = new double *[n];
  for (int i = 0; i < n; ++i) {
    inverse[i] = new double[n];
    for (int j = 0; j < n; ++j) {
      inverse[i][j] = extended_arr[i][n + j];
    }
  }
  return inverse;
}

void printMatrixToFile(double **&arr, int n, std::ofstream &out) {
  out << std::scientific;
  out.precision(10);
  out << n << " " << n << std::endl;

  for (int i = 0; i < n; ++i) {
    out << arr[i][0];
    for (int j = 1; j < n; ++j) {
      out << "   " << arr[i][j];
    }
    out << std::endl;
  }
}