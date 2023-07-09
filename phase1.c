#include "phase1.h"
#include "matrix.h"
#include "simplex.h"
#include <math.h>

int * findFeasibleBasis(DoubleMatrix A , DoubleMatrix b,  double eps, int * code)
{
    int sb = A.n;

    int * B = malloc(sb * sizeof* B);

    for (int i = 0; i < sb; i++)
    {
        B[sb-1-i] = A.m -1 - i;
    }

    DoubleMatrix Aprime = createMatrix(A.n,A.m+1);
    DoubleMatrix cprime = createMatrix(A.m+1,1);


    for (int i=0;i<A.n;i++)
    {
        
        Aprime.M[i][0] = -1;
        for (int j=0;j<A.m;j++)
        {
            Aprime.M[i][j+1] = A.M[i][j];
        }
    }

    cprime.M[0][0] = -1;
    for (int i=0;i<cprime.m;i++)
    {
        cprime.M[i+1][0] = 0; 
    }

    //choose a starting basis


    double t = INFINITY;
    int r = -1;

    for (int i = 0;i<sb;i++)
    {
        if (b.M[i][0] < 0 && b.M[i][0] < t)
        {
            t = b.M[i][0];
            r = i;
        }
    }
    

    if (r != -1)
    {
        for (int i = 0; i < r;i++)
        {
            B[i+1] = B[i];
        }
        B[0] = 0;


        for (int i=1;i<sb;i++)
        {
            B[i] += 1;
        }

        //printTab(B,sb);
        //printf("++++++++\n");

        DoubleMatrix x = simplexMethod(Aprime,b,cprime,B,sb,eps,&r);
        //printMatrix(x);
        freeMatrix(&x);

        for (int i=0;i<sb;i++)
        {
            B[i]--;
        }
        return B;
    }
    else
    {
        return B;
    }


}