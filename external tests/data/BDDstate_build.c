#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include "cudd.h"
#include "dddmp.h"

int BDDcount = 0;

int debugCheck(DdManager *manager, char* id)
{
    if (Cudd_DebugCheck(manager) > 0)
    {
        printf("Error found at %s\n", id);
        return 1;
    }
    return 0;
}

int max(int a, int b)
{
    if (a < b)
        return b;
    else
        return a;
}

int* toBitArray(int num, int len)
{
    int bits = max(ceil(log10(num)/log10(2)), len);
    int* array = (int*)malloc(sizeof(int)*bits);
    printf("Returning array for : %d, %d, %d\n", num, len, bits);
    for (int i = 0; i < bits; i++)
    {
        array[i] = 1 & (num >> i);
        printf("%d", array[i]);
    }
    if (bits < len)
    {
        for (int i = bits; i < len; i++)
        {
            array[i] = 0;
        }
    }
    printf("\n");
    return array;
}

char* toBitString(int num, int len)
{
    int bits = max(ceil(log10(num)/log10(2)), len);
    char* string = (char*)calloc(sizeof(char), bits);
    string[0] = 0;
    for (int i = 0; i < bits; i++)
        string[i] = 1 & (num >> i) ? '1' : '0';
    if (bits < len)
    {
        for (int i = bits; i < len; i++)
        {
            string[i] = '0';
        }
    }
    return string;
}

DdNode* make_state_BDD(DdManager* manager, DdNode** vars, char* encoding, int len)
{
    DdNode *state_BDD = Cudd_ReadOne(manager);
    Cudd_Ref(state_BDD);
    DdNode* tmp;
    for (int i = 0; i < len; i++)
    {
        //printf("Encoded letter: %c\n", encoding[i]);
        if (encoding[i] == '0')
        {
            DdNode* temp_neg = Cudd_Not(vars[i]);
            Cudd_Ref(temp_neg);
            tmp = Cudd_bddAnd(manager, temp_neg, state_BDD);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, state_BDD);
            Cudd_RecursiveDeref(manager, temp_neg);
        }
        else if (encoding[i] == '1')
        {
            tmp = Cudd_bddAnd(manager, vars[i], state_BDD);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, state_BDD);
        }
        else
        {
            printf("Error in encoding (1)! %c %d\n", encoding[i], len);
            return NULL;
        }
        state_BDD = tmp;
    }
    return state_BDD;
}

DdNode* make_state_BDD2(DdManager* manager, DdNode** vars, int* encoding, int len)
{
    DdNode *state_BDD = Cudd_ReadOne(manager);
    Cudd_Ref(state_BDD);
    DdNode* tmp;
    for (int i = 0; i < len; i++)
    {
        //printf("Encoded letter: %c\n", encoding[i]);
        if (encoding[i] == 0)
        {
            DdNode* temp_neg = Cudd_Not(vars[i]);
            Cudd_Ref(temp_neg);
            tmp = Cudd_bddAnd(manager, temp_neg, state_BDD);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, state_BDD);
            Cudd_RecursiveDeref(manager, temp_neg);
        }
        else if (encoding[i] == 1)
        {
            tmp = Cudd_bddAnd(manager, vars[i], state_BDD);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, state_BDD);
        }
        else
        {
            printf("Error in encoding (2)! Read %d %d\n", encoding[i], len);
            return NULL;
        }
        state_BDD = tmp;
    }
    return state_BDD;
}

DdNode* make_trans_BDD(DdManager* manager, DdNode** inVars,
    DdNode** outVars, DdNode** actionVar, char* inEnc, char* outEnc, int s_len, char* actionEnc, int a_len)
{
    DdNode* inState = make_state_BDD(manager, inVars, inEnc, s_len);

    DdNode* outState = make_state_BDD(manager, outVars, outEnc, s_len);

    DdNode* actionState = make_state_BDD(manager, actionVar, actionEnc, a_len);

    DdNode* tmp = Cudd_bddAnd(manager, actionState, outState);
    Cudd_Ref(tmp);
    DdNode* transState = Cudd_bddAnd(manager, inState, tmp);
    Cudd_Ref(transState);
    Cudd_RecursiveDeref(manager, tmp);

    Cudd_RecursiveDeref(manager, inState);
    Cudd_RecursiveDeref(manager, outState);
    Cudd_RecursiveDeref(manager, actionState);

    return transState;
}
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
    if (argc < 3)
    {
        printf("Incorrect input format! [BDD id] [encoding type={split, log}]\n");
        return 1;
    }

    char* enc_type = argv[2];
    int N;
    int vars;
    printf("started reading\n");
    scanf("%*s %d %d\n", &N, &vars);
    printf("Done reading. Read: %d, %d\n", N, vars);
    char **encodings = (char**)malloc(sizeof(char*)*(N+1));

    for (int i = 0; i < N+1; i++)
    {
        int id;
        char *encoding = (char *)calloc(sizeof(char),vars);
        scanf("%d %s\n", &id, encoding);
        encodings[i] = encoding;
        printf("states: %d, %s\n", id, encodings[i]);
    }
    for (int i = 0; i < N+1; i++)
    {
        printf("Made states: %d, %s\n", i+1, encodings[i]);
    }

    printf("done reading states\n");

    DdManager *manager; /* Global BDD manager. */
    manager = Cudd_Init(0,0,CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS,0); /* Initialize a new BDD manager. */
    debugCheck(manager, "start");
    Cudd_AutodynEnable(manager, CUDD_REORDER_SYMM_SIFT);

    int a_states = 0;
    int trans_num = 0;

    scanf("%*s %d %d\n", &a_states, &trans_num);
    printf("Reading actions: %d, %d\n", a_states, trans_num);
    int a_vars = ceil(log10(a_states)/(float)log10(2));
    DdNode **inVars = (DdNode **)malloc(sizeof(DdNode*)*vars);
    DdNode **outVars = (DdNode **)malloc(sizeof(DdNode*)*vars);
    DdNode **actionVars = (DdNode**)malloc(sizeof(DdNode*)*a_vars);
    // char **varNames = (char**)malloc(sizeof(char*)*(2*(vars+1) + a_vars));
    int** actionEncodings = (int**)malloc(sizeof(int*)*a_states);

    for (int i = 0; i < a_vars; i++)
    {
        // char *str = (char*)malloc(sizeof(char)*15);
        // sprintf(str, "a%d", i);
        // varNames[i] = str;
        // printf("ActionVars: %d\n", i);
        actionVars[i] = Cudd_bddNewVar(manager);

    }
    for (int i = 0; i < a_states; i++)
    {
        actionEncodings[i] = toBitArray(i, a_vars);
    }
    for (int i = 0; i < vars; i++)
    {
        printf("Invars: %d\n", i);
        inVars[i] = Cudd_bddNewVar(manager);
        // char numstr[25];
        // numstr[0] = 0;
        // sprintf(numstr, "%d", i);
        // char *result = (char*)malloc(sizeof(char)*(4*(vars+1)));
        // result[0] = 0;
        // strcat(result, "x");
        // strcat(result, numstr);
        // varNames[i+a_vars] = result;
    }
    for (int i = 0; i < vars; i++)
    {
        printf("Outvars: %d\n", i);
        outVars[i] = Cudd_bddNewVar(manager);
        // char numstr[25];
        // numstr[0] = 0;
        // sprintf(numstr, "%d", i);
        // char *result = (char*)malloc(sizeof(char)*(4*(vars+1)));
        // result[0] = 0;
        // strcat(result, "x");
        // strcat(result, numstr);
        // strcat(result, "\'");
        // varNames[i + vars+1 + a_vars] = result;
    }

    printf("Calculations done\n");
    DdNode* trans_system = Cudd_ReadLogicZero(manager);
    Cudd_Ref(trans_system);
    printf("Zero taken\n");
    debugCheck(manager, "var creation");
    int in = 0;
    int out = 0;
    char* a_enc = (char*)malloc(sizeof(char)*(a_vars));
    int counter = 0;
    FILE* trans_file = fopen("trans_file", "w");
    printf("trans a_l s_l\n");
    printf("%d %d %d %d\n", trans_num, a_vars, strlen(encodings[0]));
    printf("a s1 s2\n");
    for (int i = 0; i < trans_num; i++)
    {
        a_enc[0] = 0;
        scanf("%d %s %d\n", &in, a_enc, &out);
        // if ((i % 10000) == 0)
        // {
        //     printf("Rebuilding");
        //     Cudd_ReduceHeap(manager, CUDD_REORDER_ANNEALING, 500);
        // }
        if (1)
        {
            printf("Reading transition %d: %d, %s, %d\n", counter, in, a_enc, out);
            printf("With encodings: %s, %s\n", encodings[in-1], encodings[out-1]);
        }
        counter++;
        //printf("Here!\n");
        DdNode* transition = make_trans_BDD(manager, inVars, outVars, actionVars,
            encodings[in-1], encodings[out-1], vars,  a_enc, a_vars);
        //printf("Here!\n");
        DdNode* tmp = Cudd_bddOr(manager, trans_system, transition);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, trans_system);
        trans_system = tmp;
        // if (debugCheck(manager, "after transitions") > 0)
        //     return 0;
        // fflush(stdout);
        // fprintf(trans_file, "%s %s %s\n", a_enc, encodings[in-1], encodings[out-1]);
    }
    fclose(trans_file);

    // Cudd_DebugCheck(manager);

    int pg_num;
    scanf("%*s %d\n", &pg_num);
    printf("Reading progress groups: %d\n", pg_num);
    int** pg_U = (int**)malloc(sizeof(int*)*pg_num);
    int** pg_G = (int**)malloc(sizeof(int*)*pg_num);
    printf("Allocating BDDs\n");
    DdNode** BDD_pg_U = (DdNode**)malloc(sizeof(DdNode*)*pg_num);
    DdNode** BDD_pg_G = (DdNode**)malloc(sizeof(DdNode*)*pg_num);
    // Cudd_DebugCheck(manager);
    for (int i = 0; i < pg_num; i++)
    {
        int U_num;
        int G_num;
        scanf("%*s %d %d\n", &U_num, &G_num);
        pg_U[i] = calloc(sizeof(int),U_num);
        pg_G[i] = calloc(sizeof(int),G_num);
        BDD_pg_U[i] = Cudd_ReadLogicZero(manager);
        BDD_pg_G[i] = Cudd_ReadLogicZero(manager);
        DdNode* tmp;
        Cudd_Ref(BDD_pg_U[i]);
        Cudd_Ref(BDD_pg_G[i]);
        printf("Group %d: %d %d.\n", i, U_num, G_num);
        // printf("Starting reading of progs\n");
        // Cudd_DebugCheck(manager);
        for (int j = 0; j < U_num; j++)
        {
            scanf("%d", &pg_U[i][j]);
            //  printf("Read action %d\n", pg_U[i][j]);
            // for (int k = 0; k < a_vars; k++)
            // {
            //     printf("%d, %d", pg_U[i][j]-1, actionEncodings[pg_U[i][j]-1][k]);
            // }
            // printf("\n");
            DdNode* group = make_state_BDD2(manager, actionVars, actionEncodings[pg_U[i][j] - 1], a_vars);
            tmp = Cudd_bddOr(manager, group, BDD_pg_U[i]);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, BDD_pg_U[i]);
            Cudd_RecursiveDeref(manager, group);
            BDD_pg_U[i] = tmp;
            // printf("%d\n", Cudd_EquivDC(manager, group, Cudd_ReadOne(manager), Cudd_ReadLogicZero(manager)));
            // Cudd_DebugCheck(manager);
        }
        printf("Read actions.");
        // Cudd_DebugCheck(manager);
        for (int j = 0; j < G_num; j++)
        {
            scanf("%d", &pg_G[i][j]);
            DdNode* group = make_state_BDD(manager, outVars, encodings[pg_G[i][j] - 1], vars);
            tmp = Cudd_bddOr(manager, group, BDD_pg_G[i]);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, BDD_pg_G[i]);
            Cudd_RecursiveDeref(manager, group);
            BDD_pg_G[i] = tmp;

        }
        printf("Read states\n");
    }

    for (int i = 0; i < pg_num; i++)
    {
        printf("%d\n", Cudd_EquivDC(manager, BDD_pg_U[i], Cudd_ReadOne(manager), Cudd_ReadLogicZero(manager)));
        printf("%d\n", Cudd_EquivDC(manager, BDD_pg_G[i], Cudd_ReadOne(manager), Cudd_ReadLogicZero(manager)));
    }

    printf("Progress groups read\n");

    // debugCheck(manager, "Before ADD conv");
    Cudd_ReduceHeap(manager, CUDD_REORDER_SYMM_SIFT, 100);
    // FILE* file = fopen("dump.dot", "w");
    // DdNode* trans_system_ADD = Cudd_BddToAdd(manager, trans_system);
    // Cudd_Ref(trans_system_ADD);
    // debugCheck(manager, "After ADD conv");
    // // Cudd_DumpDot(manager, 1, &trans_system_ADD, (char const * const *)varNames, NULL, file);
    // fclose(file);

    char* name = calloc(sizeof(char),90);
    name[0] = 0;
    strcpy(name, "saved_BDD");
    if (argc >= 2)
    {
        strcat(name, "_");
        strcat(name, argv[1]);
        strcat(name, "_");
        strcat(name, enc_type);
    }

    FILE* save_file = fopen(name, "w");
    int* IDs = (int*)malloc(sizeof(int)*(2*vars+3));
    for (int i = 0; i < 2*vars+3; i++)
        IDs[i] = i;

    Dddmp_cuddBddStore(manager, NULL, trans_system, NULL, IDs,
                       DDDMP_MODE_BINARY, DDDMP_VARIDS, name, save_file);
    //Cudd_PrintInfo(manager, trans_system);
    fclose(save_file);

    save_file = fopen(name, "r");
    DdNode* restored_BDD = Dddmp_cuddBddLoad(manager, DDDMP_VAR_MATCHIDS, NULL, NULL,
                            NULL, DDDMP_MODE_BINARY, name, save_file);
    fclose(save_file);

    if (!Cudd_EquivDC(manager, trans_system, restored_BDD, Cudd_ReadLogicZero(manager)))
        printf("System kaput!\n");



    name[0] = 0;
    strcpy(name, "saved_U_groups");
    if (argc >= 2)
    {
        strcat(name, "_");
        strcat(name, argv[1]);
        strcat(name, "_");
        strcat(name, enc_type);
    }
    save_file = fopen(name, "w");
    Dddmp_cuddBddArrayStore(manager, NULL, pg_num, BDD_pg_U, NULL, NULL, NULL,
                            DDDMP_MODE_BINARY, DDDMP_VARIDS, NULL, save_file);
    fclose(save_file);

    save_file = fopen(name, "r");
    DdNode** U_groups;
    int prog_num = Dddmp_cuddBddArrayLoad(manager, DDDMP_ROOT_MATCHLIST, NULL,
                            DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_BINARY, NULL, save_file, &U_groups);

    fclose(save_file);
    for (int i = 0; i < prog_num; i++)
    {
        Cudd_Ref(U_groups[i]);
        // if (Cudd_EquivDC(manager, U_groups[i], Cudd_ReadOne(manager), Cudd_ReadLogicZero(manager)))
        //     printf("U NOOO! %d\n", i);
        if (!Cudd_EquivDC(manager, U_groups[i], BDD_pg_U[i], Cudd_ReadLogicZero(manager)))
            printf("U NOOOOOOO! %d\n", i);
    }

    name[0] = 0;
    strcat(name, "saved_G_groups");
    if (argc >= 2)
    {
        strcat(name, "_");
        strcat(name, argv[1]);
        strcat(name, "_");
        strcat(name, enc_type);
    }
    save_file = fopen(name, "w");
    Dddmp_cuddBddArrayStore(manager, NULL, pg_num, BDD_pg_G, NULL, NULL, NULL,
                            DDDMP_MODE_BINARY, DDDMP_VARIDS, "saved_G_groups", save_file);
    fclose(save_file);

    save_file = fopen(name, "r");
    DdNode** G_groups;
    printf("Loading array\n");
    prog_num = Dddmp_cuddBddArrayLoad(manager, DDDMP_ROOT_MATCHLIST, NULL,
                            DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_BINARY, NULL, save_file, &G_groups);
    for (int i = 0; i < prog_num; i++)
    {
        Cudd_Ref(U_groups[i]);
        if (Cudd_EquivDC(manager, G_groups[i], Cudd_ReadOne(manager), Cudd_ReadLogicZero(manager)))
            printf("G NOOO! %d\n", i);
        if (!Cudd_EquivDC(manager, G_groups[i], BDD_pg_G[i], Cudd_ReadZero(manager)))
            printf("G NOOOOOOOO! %d\n", i);
    }
    fclose(save_file);

    // save_file = fopen("saved_BDD", "r");
    // DdNode* restored_BDD = Dddmp_cuddBddLoad(manager, DDDMP_VAR_MATCHIDS, NULL, NULL,
    //                         NULL, DDDMP_MODE_BINARY, "saved_BDD", save_file);
    // fclose(save_file);
    // Cudd_Ref(restored_BDD);
    // // check equivalence of created and stored BDD
    // printf("saved == restored: %d\n", Cudd_EquivDC(manager, trans_system, restored_BDD, Cudd_ReadLogicZero(manager)));
    // printf("Saved size: %d\nRestored size: %d\n", Cudd_DagSize(trans_system), Cudd_DagSize(restored_BDD));

    Cudd_RecursiveDeref(manager, trans_system);
    // Cudd_RecursiveDeref(manager, trans_system_ADD);
    // Cudd_RecursiveDeref(manager, restored_BDD);


    printf("Invalid nodes: %d\n", Cudd_CheckZeroRef(manager));
    Cudd_Quit(manager);
    return 0;
}
// construction test
