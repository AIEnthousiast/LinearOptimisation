#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include "simplex.h"
#include "phase1.h"

int main()
{

    
    FILE * f = fopen("Phase1PL.txt","r");
    DoubleMatrix A;
    DoubleMatrix b;
    DoubleMatrix c;

    int * B;
    int sb = 0;
    int code = 0;
    

    readPLFromFile(&A,&b,&c,&B,&sb,f);
    fclose(f);


    int * C = findFeasibleBasis(A,b,1e-5,&code);
    //printTab(C,sb);
    DoubleMatrix x = simplexMethod(A,b,c,C,sb,1e-5,&code);

    printf("Solution:\n");
    printMatrix(x);
    //printTab(B,sb);

    freeMatrix(&A);
    freeMatrix(&b);
    freeMatrix(&c);

    //freeMatrix(&x);

    free(B);
    
    return EXIT_SUCCESS;
}
