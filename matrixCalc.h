#ifndef MATRIX_CALC_H
#define MATRIX_CALC_H

#include <fstream>
#include <iostream>

void create2DArray(double **&arr, int rows);
void readMatrixFromFile(double **&arr, std::ifstream &in, int rows);
void printMatrixToConsole(double **&arr, int rows);
double myAbs(double n);
double detMatrixGaus(double **arr, int n);
bool comp(char *arg, const char *tmp);
void create2DExtendedArray(double **&arr, double **&extended_arr, int rows);
bool invMatrixGausJordan(double **&extended_arr, int n);
double **getInverseFromExtended(double **extended_arr, int n);
void printMatrixToFile(double **&arr, int n, std::ofstream &out);

#endif // !MATRIX_CALC_H