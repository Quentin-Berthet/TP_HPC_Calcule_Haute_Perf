/**
 * author : Quentin Berthet
 * date : 13/04/2021
*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <limits.h>

#include "para_mmm.h"
#include "mmm_utils.h"

// Ceci est requis pour l'initialisation STATIQUE
// où NN doit être connu à la compilation ...
#define NN 16
// ... avec un NS (N squared) cohérent
#define NS 4

void showStatus(int* R, int* E, const int n) {
  int success = matrixEq(R, E, n);
  if (success) {
    printf("Mat mult was successful\n");
  } else {
    printf("Mat mult failed\n");
  }
}

// NOTE: les exits "brutaux" dans ce code ne sont
// effectivement pas super clean.
int charToIn(char* inChar) {
  char* p_err;
  long int retVal = strtol(inChar, &p_err, 10);
  if (*p_err != '\0') {
    printf("Input error detected. Exiting\n");
    exit(1);
  }
  if (retVal >= INT_MIN && retVal <= INT_MAX) {
    int k = retVal;
    return k;
  } else {
    printf("LONG INT value is too large to be assigned to an INT. Exiting\n");
    exit(1);
  }
}

int getMatrixSize(char* inChar) {
  int n = charToIn(inChar);
  if (n < 1) {
    printf("Matrix size should be a positive integer. Exiting\n");
    exit(1);
  }
  return n;
}

int main(int argc, char* argv[]) {
  int M = NS;
  int static_test = 1;
  if (argc >= 2) {
    M = getMatrixSize(argv[1]);
    static_test = 0;
  }

  printf("Test type is: %i (1 is static, otherwise it's dynamic)\n", static_test);
  printf("Matrix size is: [%ix%i]\n", M, M);

  MPI_Init(&argc, &argv);
  
  int myRank;
  int nProc;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProc);
  
  // constantes de travail
  const int ROOT = 0;
  //const int TEST_TYPE = 1; // 0 ou 1
  
  if (myRank == ROOT) {
    int* C;
    if (static_test) {
      // declaration statique
      int A[NN] = {
        1, 1, 1, 1,
        2, 2, 2, 2,
        3, 3, 3, 3,
        4, 4, 4, 4
      };
      
      int B[NN] = {
        1, 1, 1, 1,
        2, 2, 2, 2,
        3, 3, 3, 3,
        4, 4, 4, 4
      };
      
      // C devrait être:
      int E[NN] = {
        10, 10, 10, 10,
        20, 20, 20, 20,
        30, 30, 30, 30,
        40, 40, 40, 40
      };
    
      C = para_mmm(A, B, M, myRank, nProc, ROOT);
      showStatus(C, E, NS);
      free(C);
    }
    else { // i.e. static_test == 0
      int* A = generateIdentity(M);
      int* B = generateOneToNN(M);
      C = para_mmm(A, B, M, myRank, nProc, ROOT);
      showStatus(C, B, M);
      free(A);
      free(B);
      free(C);
    }
  } else {
    // on ne s'intéresse pas au retour des autres proc
    int* v = para_mmm(NULL, NULL, M, myRank, nProc, ROOT);
    if (v != NULL) {
      printf(
        "Warning: proc %i/%i returned a non NULL pointer\n",
        myRank,
        nProc
      );
    }
  }
  
  MPI_Finalize();
  
  return 0;
}

