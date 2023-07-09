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
    while (i < n && i >= 0 && array[i] != el)
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
    DoubleMatrix xB;

   
    DoubleMatrix I = inverse(Ab,1e-5);
    if (I.M == NULL)
    {
        *code = 2;
        freeMatrix(&Ab);
        
        freeMatrix(&I);
        return x;
    }
    xB = matrixMultiplication(I,b);
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
    int l = 0;
    while (1)
    {
        DoubleMatrix cnBar;
       
        DoubleMatrix An = extract(A,N,A.m - sb,0);
        DoubleMatrix cn = extract(c,N,A.m - sb,1);


        DoubleMatrix Ab = extract(A,B,sb,0);
        DoubleMatrix cB = extract(c,B,sb,1);
        DoubleMatrix AbT = transpose(Ab);
        
        DoubleMatrix y = solve(AbT,cB,eps,&i);  
        DoubleMatrix AnT = transpose(An);
        DoubleMatrix R = matrixMultiplication(AnT,y);

        
        
        
        cnBar = flop(cn,R,-1);



        freeMatrix(&AbT);
        freeMatrix(&y);
        freeMatrix(&An);
        freeMatrix(&AnT);
        freeMatrix(&cn);
        freeMatrix(&R);
        

        

        IntTabPr iTPr;
        iTPr.nbEl = 0;
        
        int coder;
        for (int i = 0; i<cnBar.n;i++)
        {

            if (cnBar.M[i][0] > eps)
            {
                insertIntTabPr(&iTPr,i,&coder);
            }
        }

        //printf("N:\n");
        //printTab(N,A.m - sb);
        //printf("CN:\n");
        //printMatrix(cnBar);
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

            DoubleMatrix comp = hstack(xB,d);
            //printf("Comp:\n");
            //printMatrix(comp);
            //printf(".............................\n\n");
            freeMatrix(&comp);
            
            if (r != -1)
            {

                for (int i=0;i<xB.n;i++)
                {
                    xB.M[i][0] -= d.M[i][0] * t;
                }


                freeMatrix(&d);

                xB.M[r][0] = t;



                int a = 0;

                //printf("Sortant:%d-Entrant:%d\n",B[r],N[e]);
                int temp = B[r];
                B[r] = N[e];
                N[e] = temp;

            


                //sort B and xB

                

                if (r > 0 && B[r] < B[r-1])
                {
                    int i = r;
                    while (i > 0  && i < sb && B[i] < B[i-1])
                    {
                        int temp = B[i-1];
                        float tempx = xB.M[i-1][0];
                        

                        B[i-1] = B[i];
                        xB.M[i-1][0] = xB.M[i][0];

                        B[i] = temp;
                        xB.M[i][0] = tempx;
                        i++;
                    }
                }
                else
                {
                    if (r < sb - 1 && B[r] > B[r+1])
                    {
                        int i = r;
                        
                        while (i < sb - 1 && i >= 0 && B[i] > B[i+1])
                        {
                            int temp = B[i+1];
                            float tempx = xB.M[i+1][0];


                            B[i+1] = B[i];
                            xB.M[i+1][0] = xB.M[i][0];

                            B[i] = temp;
                            xB.M[i][0] = tempx;
                            i++;
                        }  
                    }
                }
                //printTab(B,sb);


                //printMatrix(xB);
                //printf("-------------------\n");
                //printf("*********************\n");
                //printTab(B,sb);printf("----------------------\n");


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
            
                l++;

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