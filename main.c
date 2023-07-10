#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include "simplex.h"

int main(int argc, char ** argv)
{

    if (argc > 1)
    {
        
        int code = 0;
        DoubleMatrix *x = solveLpProblem(argv[1], &code);

        if (code == 0)
        {
            printf("Solution:\n");
            printMatrix(*x);
        }
        else
        {
            if (code = -1)
            {
                printf("Non réalisable\n");
            }
            else
            {
                printf("Non borné\n");
            }

        }
        freeMatrix(x);

    }
    return EXIT_SUCCESS;
}
