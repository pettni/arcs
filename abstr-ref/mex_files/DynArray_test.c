#include<stdio.h>
#include<stdlib.h>
#include "DynArray.h"

typedef unsigned int uint;

array_declare(NumList, uint);

void main(int argc, char* argv[])
{
    array_init(test, NumList, 10, 10);
    for (int i = 0; i < array_len(test); i++)
    {
        printf("%d\n", i);
    }
}
