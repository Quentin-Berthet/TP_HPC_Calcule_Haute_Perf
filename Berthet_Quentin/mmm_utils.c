/**
 * author : Quentin Berthet
 * date : 13/04/2021
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include "mmm_utils.h"

void print2d1d(int n)
{
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      int kc = i * n + j;
      printf("(i,j) %i, %i -> %i\n", i, j, kc);
      for (int k = 0; k < n; k++)
      {
        int ka = i * n + k;
        int kb = k * n + j;
        printf("  (i,k) %i, %i -> %i\n", i, k, ka);
        printf("  (k,j) %i, %i -> %i\n", k, j, kb);
      }
    }
  }
}

void Creat_matrix_local(int *mat_Ret, int *initMat, int numX, int numY, int sizeBaseMat, int nProc, int myRank)
{
  int jump = sqrt(sizeBaseMat);
  int nbEle = jump * jump;
  int index_X = 0;
  
  for (int i = (numX * sizeBaseMat) / jump; i < ((1 + numX) * sizeBaseMat) / jump; i++)
  {
    for (int j = (numY * sizeBaseMat) / jump; j < ((1 + numY) * sizeBaseMat) / jump; j++)
    {
      int index_in_mat = i * sizeBaseMat + j;
      mat_Ret[index_X] = initMat[index_in_mat];
      index_X += 1;
    }

  }

}

int *matMatMult(int *A, int *B, const unsigned int n)
{
  int *C = (int *)malloc(n * n * sizeof(int));
  for (size_t i = 0; i < n; i++)
  {
    for (size_t j = 0; j < n; j++)
    {
      size_t l = i * n + j;
      C[l] = 0;
      for (size_t k = 0; k < n; k++)
      {
        C[l] += A[i * n + k] * B[k * n + j];
      }
    }
  }
  return C;
}

void get_Mat_Coord(int *ret, int nProc, int myRank)
{
  int q = (int)sqrt(nProc);
  if (q * q != nProc)
  {
    printf("Nb proc is not a square: %i", nProc);
  }
  else
  {
    const int nDims = 2;
    const int nrows = q;
    const int ncols = q;

    int dims[nDims];
    dims[0] = 0;
    dims[1] = 0;

    MPI_Dims_create(nProc, nDims, dims);

    if (myRank == 0)
    {
      printf("Rank %i/%i. Grid: [%d x %d]\n", myRank, nProc, dims[0], dims[1]);
    }

    int periods[nDims];
    periods[0] = 0;
    periods[1] = 0;
    int reorder = 1;
    MPI_Comm comm2D;
    int ierr = MPI_Cart_create(MPI_COMM_WORLD, nDims, dims, periods, reorder, &comm2D);
    if (ierr != 0)
    {
      printf("Error %d while creating Cart\n", ierr);
    }

    int gridRank;
    int coord[nDims];

    MPI_Cart_coords(comm2D, myRank, nDims, coord);
    MPI_Cart_rank(comm2D, coord, &gridRank);

    ret[0] = coord[0];
    ret[1] = coord[1];
    MPI_Comm_free(&comm2D);
  }
}

void printMatrix(int *M, unsigned int n, int print_ij)
{
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      int k = i * n + j;
      if (print_ij)
      {
        printf("%i (%i, %i) ", M[k], i, j);
      }
      else
      {
        printf("%i, ", M[k]);
      }
    }
    printf("\n");
  }
}

int *generateIdentity(unsigned int n)
{
  int *I = (int *)malloc(n * n * sizeof(int));
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      int k = i * n + j;
      if (i == j)
      {
        I[k] = 1;
      }
      else
      {
        I[k] = 0;
      }
    }
  }
  return I;
}

int *generateOneToNN(unsigned int n)
{
  int *M = (int *)malloc(n * n * sizeof(int));
  for (int k = 0; k < n * n; k++)
  {
    M[k] = k;
  }
  return M;
}

int matrixEq(int *C, int *E, unsigned int n)
{
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      int k = i * n + j;
      if (C[k] != E[k])
      {
        return 0;
      }
    }
    return 1;
  }
  return 1;
}
