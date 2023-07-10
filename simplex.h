#ifndef SIMPLEX_H
#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"

#define MAX_TAB 10000
typedef struct IntTabPr
{
    int tab[MAX_TAB];
    int nbEl;
}IntTabPr;

typedef struct DoubleTabPr
{
    double tab[MAX_TAB];
    int nbEl;
}DoubleTabPr;

int inside(int * array,int n,int el);
void readPLFromFile(DoubleMatrix * A,DoubleMatrix * b,DoubleMatrix * c,int ** B, int * sb,FILE * f);
IntTabPr createIntTabPr();
DoubleTabPr createDoubleTabPr();
void insertIntTabPr(IntTabPr * T,int el,int * code);
void insertDoubleTabPr(DoubleTabPr * T,double el,int * code);
void printTab(int * T, int n);
DoubleMatrix simplexMethod(DoubleMatrix A, DoubleMatrix b, DoubleMatrix c, int * B,int sb,double eps,int * code);
void transformInequalities(DoubleMatrix * A, DoubleMatrix * b);

DoubleMatrix *  solveLpProblem(char * filename,int * code);


#endif