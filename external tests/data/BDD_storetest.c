#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "cudd.h"
#include "dddmp.h"


/**
*    TODO:
*    Test of construction of BDD
*    Test of AND OR operations on trees
*    Test of changing truth assignments
*    Test of calculating pre of things
*    Test of adding variables
*    Test of extracting truth assignments
*/

int main (int argc, char *argv[])
{
    if (argc < 2)
    {
        printf("Too few arguments! Format: [BDD file] ([transition file])\n");
        return 1;
    }

    DdManager* manager = Cudd_Init(0,0,CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS,0); /* Global BDD manager. */
    Cudd_AutodynEnable(manager, CUDD_REORDER_SYMM_SIFT);
    int* ids = (int*)malloc(sizeof(int)*35);
    for (int i = 0; i < 35; i++)
    {
        ids[i] = i;
    }

    FILE* save_file = fopen(argv[1], "r");
    DdNode* restored_BDD = Dddmp_cuddBddLoad(manager, DDDMP_VAR_MATCHIDS, NULL, NULL,
                            NULL, DDDMP_MODE_BINARY, argv[1], save_file);
    fclose(save_file);
    Cudd_Ref(restored_BDD);
    // Cudd_ReduceHeap(manager, CUDD_REORDER_SYMM_SIFT, 100);

    int trans_num, a_enc_len, s_enc_len;
    if (argc >= 3)
    {
        char* s1;
        FILE* trans_file = fopen(argv[2], "r");
        fscanf(trans_file, "%*s %*s %*s");
        fscanf(trans_file, "%d %d %d %*d", &trans_num, &a_enc_len, &s_enc_len);
        fscanf(trans_file, "%*s %*s %*s");
        for (int i = 0; i < trans_num; i++)
        {
            char* a_enc = malloc(a_enc_len);
            char *s1_enc = malloc(s_enc_len);
            char *s2_enc = malloc(s_enc_len);
            fscanf(trans_file, "%s %s %s", a_enc, s1_enc, s2_enc);
            //printf("Read %s %s %s\n", a_enc, s1_enc, s2_enc);
            int* input = (int*)malloc(sizeof(int)*(a_enc_len + 2*s_enc_len));
            for (int j = 0; j < a_enc_len; j++)
                input[j] = a_enc[j] == '1';
            for (int j = 0; j < s_enc_len; j++)
            {
                input[j+a_enc_len] = s1_enc[j] == '1';
                input[j+a_enc_len + s_enc_len] = s2_enc[j] == '1';
            }
            if (Cudd_EquivDC(manager, Cudd_Eval(manager, restored_BDD, input), Cudd_ReadZero(manager), Cudd_ReadLogicZero(manager)))
            {
                printf("transition %d existent in BDD!\n", i);
            }
        }
        fclose(trans_file);
    }

    if (argc >= 3)
    {
        printf("%d variables\n", a_enc_len + 2*s_enc_len);
        printf("%d BDD nodes vs %f possible\n", Cudd_DagSize(restored_BDD), pow(2, a_enc_len + 2*s_enc_len));
        printf("%d many nonzero paths vs %d transitions \n", (int)Cudd_CountPathsToNonZero(restored_BDD), trans_num);
        Cudd_PrintSummary(manager, restored_BDD, a_enc_len + 2*s_enc_len, 0);
    }
    else
    {
        printf("%d variables\n", Cudd_ReadSize(manager));
        printf("%d BDD nodes\n", Cudd_DagSize(restored_BDD));
        printf("%d many nonzero paths\n", (int)Cudd_CountPathsToNonZero(restored_BDD));
        Cudd_PrintSummary(manager, restored_BDD, Cudd_ReadSize(manager), 0);
    }
    printf("%d total paths\n", (int)Cudd_CountPath(restored_BDD));
    printf("Density: %f\n", Cudd_Density(manager, restored_BDD, a_enc_len + 2*s_enc_len));


    Cudd_DebugCheck(manager);

    Cudd_RecursiveDeref(manager, restored_BDD);
    Cudd_Quit(manager);
    return 0;
}
// construction test
