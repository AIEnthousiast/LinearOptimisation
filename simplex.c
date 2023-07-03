#include <stdio.h>
#include <stdlib.h>
#include "simplex.h"
#include <math.h>

IntTabPr createIntTabPr()
{
    IntTabPr T;
    T.nbEl = 0;
    return T;
}
DoubleTabPr createDoubleTabPr()
{
    DoubleTabPr T;
    T.nbEl = 0;
    return T;
}

void printTab(int * T, int n)
{
    for (int i=0;i<n;i++)
    {
        printf("%d ",T[i]);
    }
    printf("\n");
}


int getIndex(int * T,int n,int el)
{
    for (int i = 0;i<n;i++)
    {
        if (T[i] == el)
        {
            return i;
        }
    }
    
    return -1;
}

void readPLFromFile(DoubleMatrix * A,DoubleMatrix * b,DoubleMatrix * c,int ** B, int * sb,FILE * f)
{
    char temp[1000];
    int n;
    int m;
    fscanf(f,"size- %d",&n);
    fscanf(f," %d\n",&m);
    *c = createMatrix(m,1);
    fscanf(f,"c- %lf",&(c->M[0][0]));
    for (int i=1;i<m;i++)
    {
        fscanf(f," %lf",&(c->M[i][0]));
    }
    fscanf(f,"\nA:\n");
    *A = createMatrix(n,m);
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<m;j++)
        {
            fscanf(f,"%lf ",&(A->M[i][j]));
        }
        fscanf(f,"\n");
    }
    fscanf(f,"B-");
    *B = malloc(n * sizeof * B);
    for (int i=0;i<n;i++)
    {
        fscanf(f," %d",&((*B)[i]));
    }
    *b = createMatrix(n,1);
    fscanf(f,"\nb-");
    for (int i=0;i<n;i++)
    {
        fscanf(f," %lf",&(b->M[i][0]));
    }
    *sb = n;
}
void insertIntTabPr(IntTabPr * T,int el,int * code)
{
    if (T->nbEl < MAX_TAB)
    {
        T->tab[T->nbEl++] = el;
        *code = 0;
    }
    *code = 1;
}

void insertDoubleTabPr(DoubleTabPr * T,double el,int * code)
{
    if (T->nbEl < MAX_TAB)
    {
        T->tab[T->nbEl++] = el;
        *code = 0;
    }
    *code = 1;
}


int inside(int * array,int n,int el)
{
    int i = 0;
    while (i < n && array[i] != el)
    {
        i++;
    }
    return (i< n) && (array[i] == el);
}

DoubleMatrix simplexMethod(DoubleMatrix A, DoubleMatrix b, DoubleMatrix c, int * B,int sb,double eps,int * code)
{

    *code = 0;
    DoubleMatrix x = createMatrix(A.m,1);

    DoubleMatrix Ab = extract(A,B,sb,0);
    DoubleMatrix I = inverse(Ab,1e-5);

    if (I.M == NULL)
    {
        *code = 2;
        freeMatrix(&Ab);
        
        freeMatrix(&I);
        return x;
    }

    DoubleMatrix xB = matrixMultiplication(I,b);

    freeMatrix(&I);
    freeMatrix(&Ab);


    int a = 0;
    for (int i=0;i<x.n;i++)
    {
        
        if (inside(B,sb,i))
        {
            x.M[i][0] = xB.M[a][0];

            a++;
        }
        else
        {
            x.M[i][0] = 0;
        }
    }
    


    int * N = malloc((A.m - sb) * sizeof * N);
    //Construct N
    a = 0;
    int i = 0;
    while (a < A.m - sb)
    {
        if (!inside(B,sb,i))
        {
            N[a] = i;
            a++;
        }
        i++;
    }
    
    while (1)
    {
        DoubleMatrix An = extract(A,N,A.m - sb,0);
        DoubleMatrix cn = extract(c,N,A.m - sb,1);

        DoubleMatrix Ab = extract(A,B,sb,0);
        DoubleMatrix cB = extract(c,B,sb,1);
        DoubleMatrix AbT = transpose(Ab);
        
        DoubleMatrix y = solve(AbT,cB,eps,&i);
       
        DoubleMatrix AnT = transpose(An);
        DoubleMatrix R = matrixMultiplication(AnT,y);
        
        DoubleMatrix cnBar = flop(cn,R,-1);

        

        

        IntTabPr iTPr;
        iTPr.nbEl = 0;
        for (int i = 0; i<cnBar.n;i++)
        {
            if (cnBar.M[i][0] > eps)
            {
                insertIntTabPr(&iTPr,i,&i);
            }
        }

        freeMatrix(&AbT);
        freeMatrix(&y);
        freeMatrix(&An);
        freeMatrix(&AnT);
        freeMatrix(&cn);
        freeMatrix(&R);
        freeMatrix(&cnBar);


        if (iTPr.nbEl == 0)
        {
            freeMatrix(&Ab);
            freeMatrix(&cB);
            freeMatrix(&xB);
            free(N);
            return x;
        }
        else{
            int e = iTPr.tab[0];
            DoubleMatrix in = extractColumn(A,N[e]);
            DoubleMatrix d = solve(Ab,in,eps,&i);

            freeMatrix(&in);

            int r = -1;
            double t = INFINITY;
            
            for (int i=0;i<d.n;i++)
            {
                if (d.M[i][0] > 0)
                {
                    double tpot = xB.M[i][0] / d.M[i][0]; 
                    if (tpot < t)
                    {
                    t = tpot;
                    r = i; 
                    }
                }
                
            }

            freeMatrix(&Ab);
            freeMatrix(&cB);

            int a = 0;

            
            if (r != -1)
            {

                
                for (int i=0;i<xB.n;i++)
                {
                    xB.M[i][0] -= d.M[i][0] * t;
                }

                freeMatrix(&d);

                
                xB.M[r][0] = t;

                int a = 0;

                int temp = B[r];
                B[r] = N[e];
                N[e] = temp;

                for (int i=0;i<x.n;i++)
                {
                    if (!inside(B,sb,i))
                    {
                        x.M[i][0] = 0.0;
                    }
                    else
                    {
                        x.M[i][0] = xB.M[getIndex(B,sb,i)][0];
                    }
                }
            }
            else{
                *code = 1;
                freeMatrix(&xB);
                free(N);
                return x;
            }
        }
    }
    
}