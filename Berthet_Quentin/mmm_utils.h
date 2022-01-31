/**
 * author : Quentin Berthet
 * date : 13/04/2021
*/

#ifndef MMM_UTILS_H
#define MMM_UTILS_H

void print2d1d(int n);

void printMatrix(int* M, unsigned int n, int print_ij);

int* generateIdentity(unsigned int n);

void Creat_matrix_local(int *mat_Ret, int *initMat, int numX, int numY, int sizeBaseMat, int nProc, int myRank);

int* generateOneToNN(unsigned int n);

void get_Mat_Coord(int* ret, int nProc, int myRank);

int matrixEq(int* C, int* E, unsigned int n);

#endif

