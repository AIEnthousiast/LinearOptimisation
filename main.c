#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>


int main()
{
    DoubleMatrix M = createMatrix(2,2);
    DoubleMatrix N = createMatrix(2,2);
    M.M[0][0] = 1;
    M.M[0][1] = 2;
    M.M[1][0] = 3;
    M.M[1][1] = 4;

    N.M[0][0] = 1;
    N.M[0][1] = 0;
    N.M[1][0] = 0;
    N.M[1][1] = 1;



    DoubleMatrix y = createMatrix(2,1);
    y.M[0][0] = 10;
    y.M[1][0] = 5;

    int code;

    DoubleMatrix s = solve(M,y,1e-5,&code);

   // printMatrix(s);

    DoubleMatrix yprime = matrixMultiplication(M,s);

    //printMatrix(yprime);

    //freeMatrix(&y);
    //freeMatrix(&yprime);
    //freeMatrix(&s);
    //freeMatrix(&M);
    //freeMatrix(&N);
    
    //freeMatrix(&AT);
    return EXIT_SUCCESS;
}
