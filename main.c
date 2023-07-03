#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include "simplex.h"

int main()
{

    
    FILE * f = fopen("TestPL2.txt","r");
    DoubleMatrix A;
    DoubleMatrix b;
    DoubleMatrix c;

    int * B;
    int sb = 0;
    int code = 0;
    

    readPLFromFile(&A,&b,&c,&B,&sb,f);
    fclose(f);

    DoubleMatrix x = simplexMethod(A,b,c,B,sb,1e-5,&code);

    printMatrix(x);

    freeMatrix(&A);
    freeMatrix(&b);
    freeMatrix(&c);

    freeMatrix(&x);

    free(B);
    
    return EXIT_SUCCESS;
}
