#ifndef MATRIX_H
#define MATRIX_H
#include <stdio.h>
#include <stdlib.h>

typedef struct DoubleMatrix
{
    double ** M;
    int n;
    int m;
}DoubleMatrix;


DoubleMatrix matrixMultiplication(DoubleMatrix M , DoubleMatrix N);


char ** createMap2D() {
  char ** map2D = (char **) malloc(sizeof(char *)*N); // alloue le tableau 1D
  if(map2D == NULL) errorInCreate2D();                // appelle fonction qui envoie message + quit prog si NULL (malloc impossible)
  for(int i = 0; i < N; i++) {
    map2D[i] = (char *) malloc(sizeof(char)*N);       // alloue le tab2D
    if(map2D[i] == NULL) {                            // free proprement + appelle fonction qui envoie message + quit prog si NULL (malloc impossible)
      for(int j = 0; j < i; j++) {
        free(map2D[j]); 
        map2D[j] = NULL;
      }
      free(map2D);
      errorInCreate2D();
    }
  }
  return map2D;
}

DoubleMatrix createMatrix(int n, int m)
{
    DoubleMatrix M;
    M.n = n;
    M.m = m;
    M.M = malloc(n * sizeof * M.M);
    if (M.M != NULL)
    {
        for (int i=0;i<m;i++)
        {
            M.M[i] = malloc(m * sizeof * M.M[i]);
            if (M.M[i] == NULL)
            {
                for (int j=0;j<i;j++)
                {
                    free(M.M[j]);
                }
                free(M.M);
                M.M = NULL;
                return M;
            }
        }
    }
    return M;
}


#endif