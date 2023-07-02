#include "matrix.h"
#include <math.h>



DoubleMatrix createMatrix(int n, int m)
{
    DoubleMatrix M;
    M.n = n;
    M.m = m;
    M.M = malloc(n * sizeof * M.M);
    if (M.M != NULL)
    {
        for (int i=0;i<n;i++)
        {
            M.M[i] = calloc(m , sizeof * M.M[i]);
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

DoubleMatrix matrixMultiplication(DoubleMatrix M , DoubleMatrix N)
{
    DoubleMatrix A = createMatrix(M.n,N.m);

    if (A.M != NULL)
    {
        for (int i = 0;i<M.n;i++)
        {
            for (int j = 0; j < N.m;j++)
            {
                double sum = 0.0;
                for (int k = 0; k < M.m;k++)
                {
                    sum += M.M[i][k] * N.M[k][j];
                }
                A.M[i][j] = sum;    
            }
        }
    }
    
    return A;
}

DoubleMatrix transpose(DoubleMatrix M)
{
    DoubleMatrix A = createMatrix(M.m,M.n);

    if (A.M != NULL)
    {
        for (int i=0;i<M.m;i++)
        {
            for (int j=0;j<M.n;j++)
            {
                A.M[i][j] = M.M[j][i];
            }
        }
    }
    return A;
}

void printMatrix(DoubleMatrix M)
{
    if (M.M != NULL)
    {
        for (int i=0;i<M.n;i++)
        {
            for (int j=0;j<M.m;j++)
            {
                printf("%lf ",M.M[i][j]);
            }
            printf("\n");
        }
    }
}
DoubleMatrix copyMatrix(DoubleMatrix M)
{
    DoubleMatrix A = createMatrix(M.n,M.m);

    for (int i=0;i<M.n;i++)
    {
        for (int j=0;j<M.m;j++)
        {
            A.M[i][j] = M.M[i][j];
        }
    }
    return A;
}
DoubleMatrix reducedEchelonForm(DoubleMatrix M,double eps)
{
    DoubleMatrix A = copyMatrix(M);

    if (A.M != NULL)
    {
        for (int i = 0; i < A.n - 1;i++)
        {
            // search for the largest pivot
            int pivotIndex = i;
            for (int j=i+1;j<A.n;j++)
            {
                if (fabs(M.M[pivotIndex][i]) < abs(M.M[j][i]))
                {
                    pivotIndex = j;
                }
            }

            exchangeRow(&A,i,pivotIndex);
            if (A.M[i][i] > eps)
            {
                for (int j=i+1;j<A.n;j++)
                {
                    double a = A.M[j][i] / A.M[i][i];
                    for (int k=i;k<A.m;k++)
                    {

                        A.M[j][k] -= a * A.M[i][k];
                    }
                }
            }
        }
    }

    return A;
}

void exchangeRow(DoubleMatrix * M,int i, int j)
{
    if (i < M->n && j < M->n)
    {
        double * temp = M->M[i];
        M->M[i] = M->M[j];
        M->M[j] = temp;
    }
}


DoubleMatrix hstack(DoubleMatrix M1, DoubleMatrix M2)
{
    DoubleMatrix A = createMatrix(M1.n, M1.m + M2.m);


    if (A.M != NULL)
    {
        for (int i = 0; i < M1.n;i++)
        {
            for (int j = 0; j < M1.m; j++)
            {
                A.M[i][j] = M1.M[i][j]; 
            }
        }
        for (int i = 0; i < M2.n;i++)
        {
            for (int j = 0; j < M2.m; j++)
            {
                A.M[i][M1.m + j] = M2.M[i][j];
            }
        }
    }

    return A;
}

DoubleMatrix solve(DoubleMatrix M, DoubleMatrix b,double eps,int *code)
{
    *code = 0; 
    DoubleMatrix solution = createMatrix(M.n,1);
    printf("Solution:\n");
    printMatrix(solution);
    DoubleMatrix Augmented = hstack(M,b);
    printf("\nAugmented:\n\n");
    printMatrix(Augmented);

    DoubleMatrix ReducedAugmented = reducedEchelonForm(Augmented,eps);
    printf("ReducedAugmented:\n\n");
    printMatrix(ReducedAugmented);
    printf("\n\n");
    
    for (int i = 0; i < ReducedAugmented.n;i++)
    {
        if (fabs(ReducedAugmented.M[i][i]) <= eps)
        {
            *code = 1;
        }
    }
    

    if (*code == 0)
    {

        for (int i = ReducedAugmented.n - 1; i > -1;i--)
        {
            solution.M[i][0] = ReducedAugmented.M[i][M.m];
            for (int j = ReducedAugmented.n-1; j > i;j--)
            {
                solution.M[i][0] -= ReducedAugmented.M[i][j] * solution.M[j][0];
            }
            
            
            solution.M[i][0] = solution.M[i][0] / ReducedAugmented.M[i][i];
        }
    }
    
    freeMatrix(&Augmented);
    freeMatrix(&ReducedAugmented);
    return solution;

}   

void freeMatrix(DoubleMatrix * M)
{
    if (M->M != NULL)
    {
        for (int i=0;i<M->n;i++)
        {
            free(M->M[i]);
        }
        free(M->M);
    }
}