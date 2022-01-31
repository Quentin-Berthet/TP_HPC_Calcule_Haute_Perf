/**
 * author : Quentin Berthet
 * date : 13/04/2021
*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "mmm_utils.h"
#include "para_mmm.h"

// implementation de la multiplication
// matrice-matrice carrée
typedef struct sous_mat_t
{
  int *sous_mat;
  int coodMat[2];
} sous_mat_t;

void display_groupe(int groupe, int myRank, int nProc, int rRank)
{
  if (groupe == 0)
  {
    printf("WORLD GROUPE RED RANK/SIZE: %d/%d --- ROW RANK/SIZE: %d/1\n", myRank, nProc, rRank);
  }
  else
  {
    printf("WORLD GROUPE BLUE RANK/SIZE: %d/%d --- ROW RANK/SIZE: %d/2\n", myRank, nProc, rRank);
  }
  printf("I'm proc %i/%i and I recv matrix C\n", myRank, nProc);
}

void print_mat(int *mat, int sizex, int sizey)
{
  printf("+++++++++++++++++++++++++++++\n");
  for (int i = 0; i < sizex; i++)
  {
    for (int j = 0; j < sizey; j++)
    {
      int k = i * sizey + j;
      printf("%d ", mat[k]);
    }
    printf("\n");
  }
  printf("+++++++++++++++++++++++++++++\n");
}

int verif_out_of_band(int val, int size_under_mat, int myRank)
{
  int size_tot = size_under_mat * size_under_mat;

  if (val < 0)
  {
    return val + (size_tot);
  }
  else if (val > size_tot - 1)
  {
    return size_under_mat - ((size_tot)-myRank);
  }
  else
  {
    return val;
  }
}
// cette fonction est le point d'entrée de votre
// algorithme parallèle
int *para_mmm(
    int *A,
    int *B,
    const unsigned int n,
    const int myRank,
    const int nProc,
    const int root)
{

  const int DUMMY_TAG = 0;
  const unsigned int NN = n * n;
  MPI_Status status;

  if (nProc < 2)
  {
    fprintf(stderr, "Less than 2 ranks available: %i\n", nProc);
    return NULL;
  }
  else
  {

    if (myRank == root)
    {
      int *C = (int *)malloc(NN * sizeof(int));

      MPI_Request request = MPI_REQUEST_NULL;
      int *under_mat_A = (int *)malloc(n * sizeof(int));
      int *under_mat_B = (int *)malloc(n * sizeof(int));
      int *under_mat_C = (int *)malloc(n * sizeof(int));
      int size_under_mat = sqrt(n);
      int number_of_element_under_mat = size_under_mat * size_under_mat;
      int *ensemble_coord_Mat = (int *)malloc(2 * nProc * sizeof(int));
      int *coordMat = malloc(size_under_mat * sizeof(int));
      int *under_mat_C_General = (int *)malloc(size_under_mat * number_of_element_under_mat * sizeof(int));
      get_Mat_Coord(coordMat, nProc, myRank);
      //récupère coordonnée cartésienne de tout les sous matrice
      MPI_Gather(coordMat, 2, MPI_INT, ensemble_coord_Mat, 2, MPI_INT, root, MPI_COMM_WORLD);
      int index_p = 1;
      int index_coor = 2;
      int index_sender = 0;
      int tag = 1234;
      //creer sous matrice A et B du process 0
      Creat_matrix_local(under_mat_A, A, ensemble_coord_Mat[0], ensemble_coord_Mat[1], n, nProc, myRank);
      Creat_matrix_local(under_mat_B, B, ensemble_coord_Mat[0], ensemble_coord_Mat[1], n, nProc, myRank);

      //creer sous matrice A et B des autres process et l'envoie a ce dernier
      for (int i = 1; i < nProc; i++)
      {
        int *under_mat_A_tmp = (int *)malloc(n * sizeof(int));
        int *under_mat_B_tmp = (int *)malloc(n * sizeof(int));
        Creat_matrix_local(under_mat_A_tmp, A, ensemble_coord_Mat[index_coor], ensemble_coord_Mat[index_coor + 1], n, nProc, myRank);
        Creat_matrix_local(under_mat_B_tmp, B, ensemble_coord_Mat[index_coor], ensemble_coord_Mat[index_coor + 1], n, nProc, myRank);
        MPI_Send(under_mat_A_tmp, n, MPI_INT, index_p, 0, MPI_COMM_WORLD);
        MPI_Send(under_mat_B_tmp, n, MPI_INT, index_p, 0, MPI_COMM_WORLD);
        free(under_mat_A_tmp);
        free(under_mat_B_tmp);
        index_p += 1;
        index_coor += 2;
      }
      //créer un communicator pour chaque ligne de la matrice (cela servira a effectuer le boradcast de la sous matrice A)
      MPI_Comm newcomm;
      MPI_Comm_split(MPI_COMM_WORLD, coordMat[0], myRank, &newcomm);
      int rRank, rSize;
      MPI_Comm_rank(newcomm, &rRank);
      MPI_Comm_size(newcomm, &rSize);
      int test_size = size_under_mat - 1;

      for (int k = 0; k < size_under_mat; k++)
      {
        //si c'est a notre process de broadcaster la sous matrice A
        if (coordMat[0] == coordMat[1] + k || coordMat[0] + k == coordMat[1])
        {
          index_sender = coordMat[1];

          MPI_Bcast(&index_sender, 1, MPI_INT, coordMat[1], newcomm);
          MPI_Bcast(under_mat_A, number_of_element_under_mat, MPI_INT, root, newcomm);
          int *tmp_mat_C = matMatMult(under_mat_A, under_mat_B, size_under_mat);
          for (int i = 0; i < size_under_mat; i++)
          {
            for (int j = 0; j < size_under_mat; j++)
            {
              int v = i * size_under_mat + j;
              under_mat_C_General[v + (k * n)] = tmp_mat_C[v];
            }
          }
          free(tmp_mat_C);
        }
        else
        {
          //recoit broadcast de la sous matrice a
          MPI_Bcast(&index_sender, 1, MPI_INT, coordMat[0], newcomm);
          MPI_Bcast(under_mat_A, number_of_element_under_mat, MPI_INT, index_sender, newcomm);
          int *tmp_mat_C = matMatMult(under_mat_A, under_mat_B, size_under_mat);
          for (int i = 0; i < size_under_mat; i++)
          {
            for (int j = 0; j < size_under_mat; j++)
            {
              int v = i * size_under_mat + j;
              under_mat_C_General[v + (k * n)] = tmp_mat_C[v];
            }
          }
          free(tmp_mat_C);
        }
        //définit le numéro du process auquel on va envoyer notre sous matrice B et celui du process qui nous envoye sa sous matrice B
        int next = verif_out_of_band((myRank + size_under_mat), size_under_mat, myRank);
        int befor = verif_out_of_band((myRank - size_under_mat), size_under_mat, myRank);
        //envoie et recoie la sous matrice B
        MPI_Send(under_mat_B, number_of_element_under_mat, MPI_INT, befor, tag, MPI_COMM_WORLD);
        MPI_Recv(under_mat_B, number_of_element_under_mat, MPI_INT, next, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }

      // add toute les matrice qui se trouve dans under_mat_C_General a C
      int *tmp_C = (int *)malloc(NN * sizeof(int));
      for (int i = 0; i < number_of_element_under_mat; i++)
      {
        for (int j = 0; j < size_under_mat; j++)
        {
          int p = i + (j * number_of_element_under_mat);
          under_mat_C[i] += under_mat_C_General[p];
        }
      }
      // récupere toute les sous matrice de c
      MPI_Gather(under_mat_C, number_of_element_under_mat, MPI_INT, tmp_C, number_of_element_under_mat, MPI_INT, root, MPI_COMM_WORLD);
      int indice = 0;
      int *index_mat = generateOneToNN(n);
      int *reorder_mat = (int *)malloc(NN * sizeof(int));
      int num = 0;
      index_coor = 0;
      //creation matrice contenant les indice pour réaranger les sous matrice c dans la matrice C
      for (int i = 0; i < nProc; i++)
      {
        int *under_mat_tmp = (int *)malloc(n * sizeof(int));
        Creat_matrix_local(under_mat_tmp, index_mat, ensemble_coord_Mat[index_coor], ensemble_coord_Mat[index_coor + 1], n, nProc, i);

        for (int j = 0; j < size_under_mat; j++)
        {
          for (int n = 0; n < size_under_mat; n++)
          {
            int k = j * size_under_mat + n;
            reorder_mat[num] = under_mat_tmp[k];
            num += 1;
          }
        }
        index_coor += 2;
      }

      //recréer la matrice c en utilisant les sous matrice de c

      for (int i = 0; i < number_of_element_under_mat; i++)
      {
        for (int j = 0; j < nProc; j++)
        {

          C[reorder_mat[indice]] = tmp_C[indice];

          indice += 1;
        }
      }

      //on récupère toute les sous matrice de c et maintenant reste plus qu'a réorganiser cette matrice
      print_mat(C, n, n);

      free(tmp_C);
      free(coordMat);
      free(under_mat_C_General);
      free(ensemble_coord_Mat);
      free(under_mat_A);
      free(under_mat_B);
      free(under_mat_C);
      MPI_Comm_free(&newcomm);
      return C;
    }
    else
    {
      int size_under_mat = sqrt(n);
      int number_of_element_under_mat = size_under_mat * size_under_mat;
      MPI_Request request = MPI_REQUEST_NULL;
      int *under_mat_A = (int *)malloc(n * sizeof(int));
      int *under_mat_B = (int *)malloc(n * sizeof(int));
      int *under_mat_C = (int *)malloc(n * sizeof(int));
      int index_sender = 0;
      int *coordMat = malloc(size_under_mat * sizeof(int));
      //determine coordonée cartésienne de la sous matrice
      get_Mat_Coord(coordMat, nProc, myRank);
      int tag = 1234;
      //envoie au process root
      MPI_Gather(coordMat, 2, MPI_INT, coordMat, 2, MPI_INT, root, MPI_COMM_WORLD);
      //récupère sous matrice A et B 
      MPI_Recv(under_mat_A, n, MPI_INT, root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(under_mat_B, n, MPI_INT, root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      int *under_mat_C_General = (int *)malloc(size_under_mat * number_of_element_under_mat * sizeof(int));
      //créer un communicator pour chaque ligne de la matrice (cela servira a effectuer le boradcast de la sous matrice A)
      MPI_Comm newcomm;
      MPI_Comm_split(MPI_COMM_WORLD, coordMat[0], myRank, &newcomm);
      int rRank, rSize;
      MPI_Comm_rank(newcomm, &rRank);
      MPI_Comm_size(newcomm, &rSize);
      int test_size = size_under_mat - 1;
      for (int k = 0; k < size_under_mat; k++)
      {
        //algorithme de fox Boradcast sous matrice A et envoie ainsi que récéption des sous matrice B
        if (coordMat[0] == coordMat[1] + k || coordMat[0] + k == coordMat[1])
        {
          index_sender = coordMat[1];
          MPI_Bcast(&index_sender, 1, MPI_INT, coordMat[0], newcomm);
          MPI_Bcast(under_mat_A, number_of_element_under_mat, MPI_INT, coordMat[1], newcomm);
          int *tmp_mat_C = matMatMult(under_mat_A, under_mat_B, size_under_mat);
          for (int i = 0; i < size_under_mat; i++)
          {
            for (int j = 0; j < size_under_mat; j++)
            {
              int v = i * size_under_mat + j;
              under_mat_C_General[v + (k * n)] = tmp_mat_C[v];
            }
          }
          free(tmp_mat_C);
          int next = verif_out_of_band((myRank + size_under_mat), size_under_mat, myRank);
          int befor = verif_out_of_band((myRank - size_under_mat), size_under_mat, myRank);
          MPI_Send(under_mat_B, number_of_element_under_mat, MPI_INT, befor, tag, MPI_COMM_WORLD);
          MPI_Recv(under_mat_B, number_of_element_under_mat, MPI_INT, next, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else
        {
          //algorithme de fox Boradcast sous matrice A et envoie ainsi que récéption des sous matrice B
          MPI_Bcast(&index_sender, 1, MPI_INT, coordMat[0], newcomm);
          MPI_Bcast(under_mat_A, number_of_element_under_mat, MPI_INT, index_sender, newcomm);

          int *tmp_mat_C = matMatMult(under_mat_A, under_mat_B, size_under_mat);
          for (int i = 0; i < size_under_mat; i++)
          {
            for (int j = 0; j < size_under_mat; j++)
            {
              int v = i * size_under_mat + j;
              under_mat_C_General[v + (k * n)] = tmp_mat_C[v];
            }
          }
          free(tmp_mat_C);
          int next = verif_out_of_band((myRank + size_under_mat), size_under_mat, myRank);
          int befor = verif_out_of_band((myRank - size_under_mat), size_under_mat, myRank);
          MPI_Send(under_mat_B, number_of_element_under_mat, MPI_INT, befor, tag, MPI_COMM_WORLD);
          MPI_Recv(under_mat_B, number_of_element_under_mat, MPI_INT, next, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
      }
      for (int i = 0; i < number_of_element_under_mat; i++)
      {
        for (int j = 0; j < size_under_mat; j++)
        {
          int p = i + (j * number_of_element_under_mat);
          under_mat_C[i] += under_mat_C_General[p];
        }
      }
      //envoie la sous matric c de son process au process root 
      MPI_Gather(under_mat_C, number_of_element_under_mat, MPI_INT, under_mat_C, number_of_element_under_mat, MPI_INT, root, MPI_COMM_WORLD);

      free(under_mat_C);
      free(under_mat_B);
      free(coordMat);
      free(under_mat_A);
      free(under_mat_C_General);
      MPI_Comm_free(&newcomm);
      return NULL;
    }
  }
}
