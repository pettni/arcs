#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
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

#define PRINT_A_IN(manager, system, a_num, var_num, showBin) { \
    int interval[] = {0, a_num + var_num};\
    int delims[] = {a_num};\
    printStatesIn(manager, system, interval, delims, 1, showBin);\
}
#define PRINT_IN(manager, system, a_num, var_num, showBin) {\
    int interval[] = {a_num, a_num + var_num};\
    printStatesIn(manager, system, interval, NULL, 0, showBin);\
}
#define PRINT_OUT(manager, system, a_num, var_num, showBin) {\
    int interval[] = {a_num + var_num, a_num + 2*var_num};\
    printStatesIn(manager, system, interval, NULL, 0, showBin);\
}
#define PRINT_A(manager, system, a_num, var_num, showBin) {\
    int interval[] = {0, a_num};\
    printStatesIn(manager, system, interval, NULL, 0, showBin);\
}

void printState(int* enc, int num, int* start, int i, int* end, int* delims, int delim_ind)
{
    // printf("%d\n", i);
    if (i == *start)
        printf("(");
    if (i == *end-1)
    {
        printf("%d) ", num + 1);
        return;
    }
    else if (delim_ind != -1 && i == delims[delim_ind]-1)
    {
        printf("%d,", num + 1);
        num = 0;
        delim_ind--;
    }

    if (enc[i] == 1 || enc[i] == 0)
        printState(enc, 2*num + enc[i], start, i-1, end, delims, delim_ind);
    else
    {
        printState(enc, 2*num + 1, start, i-1, end, delims, delim_ind);
        printf("(");
        printState(enc, 2*num, start, i-1, end, delims, delim_ind);
    }
}

void printStatesIn(DdManager* manager, DdNode* BDD, int* interval, int* delims, int delim_num, int showBin)
{
    int* cube;
    DdGen* gen;
    CUDD_VALUE_TYPE value;
    Cudd_ForeachCube(manager, BDD, gen, cube, value)
    {
        int start = interval[1]-1;
        int end = interval[0];
        if (showBin)
        {
            if (cube == NULL)
            {
                printf("Invalid system\n");
                return;
            }
            for (int j = interval[0]; j < interval[1]; j++)
            {
                printf("%d", cube[j]);
            }
            printf("=");
        }
        printState(cube, 0, &start, start, &end, delims, delim_num-1);
    }
    printf("\n");
}

int max(int a, int b)
{
    if (a < b)
        return b;
    else
        return a;
}

DdNode* pre(DdManager* manager, DdNode** outVars, DdNode** inVars, int Var_num, DdNode** aVars, int aVar_num,
     DdNode* all_out, DdNode* T, DdNode* X_prime, DdNode* A, char quant1, char quant2)
{
    DdNode* all_in = Cudd_bddSwapVariables(manager, all_out, outVars, inVars, Var_num);
    Cudd_Ref(all_in);
    DdNode* valid_transitions = T;
    DdNode* tmp;

    DdNode* outCube = Cudd_bddComputeCube(manager, outVars, NULL, Var_num);
    Cudd_Ref(outCube);
    DdNode* aCube = Cudd_bddComputeCube(manager, aVars, NULL, aVar_num);
    Cudd_Ref(aCube);

    // Cudd_DebugCheck(manager);
    if (quant1 == 'e' && quant2 == 'e')
    {
        valid_transitions = Cudd_bddAnd(manager, T, X_prime);
        Cudd_Ref(valid_transitions);
        tmp = Cudd_bddExistAbstract(manager, valid_transitions, outCube);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, valid_transitions);
        valid_transitions = tmp;

        tmp = Cudd_bddAnd(manager, valid_transitions, A);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, valid_transitions);
        valid_transitions = tmp;
        tmp = Cudd_bddExistAbstract(manager, valid_transitions, aCube);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, valid_transitions);
        valid_transitions = tmp;
    }
    else if (quant1 == 'e' && quant2 == 'a')
    {
        DdNode* not_X = Cudd_Not(X_prime);
        Cudd_Ref(not_X);
        DdNode* state_a_pairs_not_X = Cudd_bddAndAbstract(manager, T, not_X, outCube);
        Cudd_Ref(state_a_pairs_not_X);
        Cudd_RecursiveDeref(manager, not_X);

        tmp = Cudd_Not(state_a_pairs_not_X);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, state_a_pairs_not_X);
        state_a_pairs_not_X = tmp;

        DdNode* state_a_pairs_X = Cudd_bddExistAbstract(manager, T, outCube);
        Cudd_Ref(state_a_pairs_X);

        valid_transitions = Cudd_bddAnd(manager, state_a_pairs_X, state_a_pairs_not_X);
        Cudd_Ref(valid_transitions);
        Cudd_RecursiveDeref(manager, state_a_pairs_X);
        Cudd_RecursiveDeref(manager, state_a_pairs_not_X);

        tmp = Cudd_bddAndAbstract(manager, valid_transitions, A, aCube);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, valid_transitions);
        valid_transitions = tmp;
    }
    else if(quant1 == 'a' && quant2 == 'a')
    {
        DdNode* not_X = Cudd_Not(X_prime);
        Cudd_Ref(not_X);
        DdNode* state_a_pairs_not_X = Cudd_bddAndAbstract(manager, T, not_X, outCube);
        Cudd_Ref(state_a_pairs_not_X);
        Cudd_RecursiveDeref(manager, not_X);

        tmp = Cudd_Not(state_a_pairs_not_X);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, state_a_pairs_not_X);
        state_a_pairs_not_X = tmp;

        DdNode* state_a_pairs_X = Cudd_bddExistAbstract(manager, T, outCube);
        Cudd_Ref(state_a_pairs_X);

        valid_transitions = Cudd_bddAnd(manager, state_a_pairs_X, state_a_pairs_not_X);
        Cudd_Ref(valid_transitions);
        Cudd_RecursiveDeref(manager, state_a_pairs_X);
        Cudd_RecursiveDeref(manager, state_a_pairs_not_X);

        DdNode* not_A = Cudd_Not(A);
        Cudd_Ref(not_A);

        tmp = Cudd_bddOr(manager, valid_transitions, not_A);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, not_A);
        Cudd_RecursiveDeref(manager, valid_transitions);
        valid_transitions = tmp;

        tmp = Cudd_bddUnivAbstract(manager, valid_transitions, aCube);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, valid_transitions);
        valid_transitions = tmp;
    }
    else if (quant1 == 'a' && quant2 == 'e')
    {
        tmp = Cudd_bddAndAbstract(manager, T, X_prime, outCube);
        Cudd_Ref(tmp);
        valid_transitions = tmp;

        // printf("E(T & X)\n");
        // PRINT_A_IN(manager, valid_transitions, aVar_num, Var_num, 0);

        DdNode* not_A = Cudd_Not(A);
        Cudd_Ref(not_A);
        tmp = Cudd_bddOr(manager, valid_transitions, not_A);
        Cudd_Ref(tmp);

        // printf("-A\n");
        // PRINT_A(manager, not_A, aVar_num, Var_num, 0);

        // printf("-A | E(T & X)\n");
        // PRINT_A_IN(manager, tmp, aVar_num, Var_num, 0);

        Cudd_RecursiveDeref(manager, valid_transitions);
        Cudd_RecursiveDeref(manager, not_A);
        valid_transitions = tmp;

        tmp = Cudd_bddUnivAbstract(manager, valid_transitions, aCube);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, valid_transitions);
        valid_transitions = tmp;
    }
    else
    {
        printf("Invalid quantifiers!\n");
        return NULL;
    }

    Cudd_RecursiveDeref(manager, outCube);
    Cudd_RecursiveDeref(manager, aCube);


    return valid_transitions;
}

int* toBitArray(int num, int len)
{
    int bits = max(ceil(log10(num)/log10(2)), len);
    int* array = (int*)malloc(sizeof(int)*bits);
    //printf("Returning array for : %d, %d\n", num, len);
    for (int i = 0; i < bits; i++)
    {
        array[i] = 1 & (num >> i);
        //printf("%d", array[i]);
    }
    if (bits < len)
    {
        for (int i = bits; i < len; i++)
        {
            array[i] = 0;
        }
    }
    //printf("\n");
    return array;
}

DdNode* makeSet(DdManager* manager, DdNode** vars, int var_num, int* inds, int n, int** encodings)
{
    DdNode* set = Cudd_ReadLogicZero(manager);
    Cudd_Ref(set);
    for (int i = 0; i < n; i++)
    {
        int* enc = encodings[inds[i]];
        DdNode* cube = Cudd_ReadOne(manager);
        Cudd_Ref(cube);
        for (int j = var_num-1; j >= 0; j--)
        {
            DdNode* var = enc[j] ? vars[j] : Cudd_Not(vars[j]);
            if (!enc[j])
                Cudd_Ref(var);
            DdNode* tmp = Cudd_bddAnd(manager, cube, var);
            Cudd_Ref(tmp);
            if (!enc[j])
                Cudd_RecursiveDeref(manager, var);
            Cudd_RecursiveDeref(manager, cube);
            cube = tmp;
        }
        DdNode* tmp = Cudd_bddOr(manager, set, cube);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, cube);
        Cudd_RecursiveDeref(manager, set);
        set = tmp;
    }
    //Cudd_DebugCheck(manager);
    return set;
}

int main (int argc, char *argv[])
{
    if (argc < 3)
    {
        printf("Too few arguments! Format: [BDD file] [BDD text file] [BDD test file]\n");
        return 1;
    }

    DdManager* manager = Cudd_Init(0,0,CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS,0); /* Global BDD manager. */
    Cudd_AutodynEnable(manager, CUDD_REORDER_SYMM_SIFT);
    // int* ids = (int*)malloc(sizeof(int)*35);
    // for (int i = 0; i < 35; i++)
    // {
    //     ids[i] = i;
    // }

    printf("Loading BDD\n");
    FILE* save_file = fopen(argv[1], "r");
    DdNode* restored_BDD = Dddmp_cuddBddLoad(manager, DDDMP_VAR_MATCHIDS, NULL, NULL,
                            NULL, DDDMP_MODE_BINARY, argv[1], save_file);
    fclose(save_file);
    Cudd_Ref(restored_BDD);
    printf("BDD loaded : %d vars\n", Cudd_ReadSize(manager));
    // Cudd_ReduceHeap(manager, CUDD_REORDER_SYMM_SIFT, 100);

    int** encodings;
    int** action_encodings;
    int var_num, a_var_num, state_num, action_num;
    DdNode** state_in_vars;
    DdNode** state_out_vars;
    DdNode** action_vars;
    printf("Loading BDD parameters\n");
    if (argc >= 3)
    {
        FILE* trans_file = fopen(argv[2], "r");
        // printf("Parameter file loaded\n");
        fscanf(trans_file, "%*c %d %d", &state_num, &var_num);
        // printf("Read params: state_num=%d var_num=%d\n", state_num, var_num);
        // printf("Read state and action parameters\n");
        state_in_vars = malloc(sizeof(DdNode*)*var_num);
        state_out_vars = malloc(sizeof(DdNode*)*var_num);
        encodings = (int**)malloc(sizeof(int*)*state_num+1);
        int* visited = calloc(sizeof(int*), pow(2,var_num+1));
        int not_unique = 0;
        // printf("Allocated space\n");
        for (int i = 0; i < state_num+1; i++)
        {
            // printf("Encoding %d loading\n", i);
            char* enc = malloc(sizeof(char)*var_num);
            memset(enc, '0', var_num);
            fscanf(trans_file, "%*d %s", enc);
            int map_id = fromBitCharArray(enc, var_num);
            if (visited[map_id])
                not_unique = 1;
            else
                visited[map_id] = 1;

            // printf("Encoding %d read from file loading\n", i);
            int* enc_num = malloc(sizeof(int)*var_num);
            for (int j = 0; j < var_num; j++)
                enc_num[j] = enc[j] == '1';
            encodings[i] = enc_num;
            printf("Read encoding %d = %s\n", i, enc);
            free(enc);
        }
        if (not_unique)
            printf()

        printf("Encodings loaded\n");
        fscanf(trans_file, "%*s %d %*d", &action_num);

        a_var_num = ceil(log10(action_num)/log10(2));
        action_encodings = malloc(sizeof(int*)*action_num);
        for (int i = 0; i < action_num; i++)
            action_encodings[i] = toBitArray(i, a_var_num);
        action_vars = malloc(sizeof(DdNode*)*a_var_num);

        fclose(trans_file);

        for (int i = 0; i < a_var_num; i++)
        {
            action_vars[i] = Cudd_bddIthVar(manager, i);
        }
        for (int i = 0; i < var_num; i++)
        {
            DdNode* var = Cudd_bddIthVar(manager, i+a_var_num);
            if (var == NULL)
                printf("HELP %d\n", i+a_var_num);
            state_in_vars[i] = var;
        }
        for (int i = 0; i < var_num; i++)
        {
            state_out_vars[i] = Cudd_bddIthVar(manager, i+a_var_num + var_num);
        }
    }
    printf("BDD parameters loaded : %d vars\n", Cudd_ReadSize(manager));

    printf("Creating complete state set\n");
    DdNode* all_out_states;
    DdNode* inCube = Cudd_bddComputeCube(manager, state_in_vars, NULL, var_num);
    Cudd_Ref(inCube);
    // PRINT_IN(manager, inCube, a_var_num, var_num, 1);


    DdNode* actionCube = Cudd_bddComputeCube(manager, action_vars, NULL, a_var_num);
    Cudd_Ref(actionCube);
    DdNode* cube = Cudd_bddAnd(manager, inCube, actionCube);
    Cudd_Ref(cube);

    Cudd_RecursiveDeref(manager, inCube);
    Cudd_RecursiveDeref(manager, actionCube);

    all_out_states = Cudd_bddExistAbstract(manager, restored_BDD, cube);
    Cudd_Ref(all_out_states);
    Cudd_RecursiveDeref(manager, cube);


    //PRINT_OUT(manager, all_out_states, a_var_num, var_num, 1);
    printf("Set created\n");

    int number_of_states;
    int number_of_actions;
    int number_of_tests;
    int** test_states;
    int* test_state_num;
    int** test_actions;
    int* test_action_num;
    char** test_types;
    float* pre_e_e_times;
    int pre_e_e_num = 0;
    float* pre_e_a_times;
    int pre_e_a_num = 0;
    float* pre_a_e_times;
    int pre_a_e_num = 0;
    float* pre_a_a_times;
    int pre_a_a_num = 0;
    printf("Loading test parameters\n");
    FILE* file = fopen(argv[3], "r");
    fscanf(file, "%d %d %d", &number_of_states, &number_of_actions, &number_of_tests);
    printf("Read test parameters: number_of_states=%d, number_of_actions=%d, number_of_tests=%d\n", number_of_states, number_of_actions, number_of_tests);
    test_states = (int**)malloc(sizeof(int*)*number_of_tests);
    test_actions = (int**)malloc(sizeof(int*)*number_of_tests);
    test_types = (char**)malloc(sizeof(char*)*number_of_tests);
    test_state_num = (int*)malloc(sizeof(int)*number_of_tests);
    test_action_num = (int*)malloc(sizeof(int)*number_of_tests);
    for (int i = 0; i < number_of_tests; i++)
    {
        printf("Loading test %d\n", i);
        int states_in_test, actions_in_test;
        test_types[i] = malloc(sizeof(char)*2);
        fscanf(file, "%*s %d %*s %d %c %c", &test_state_num[i], &test_action_num[i], &test_types[i][0], &test_types[i][1]);
        printf("Read %d %d %c %c\n", test_state_num[i], test_action_num[i], test_types[i][0], test_types[i][1]);
        test_states[i] = malloc(sizeof(int)*test_state_num[i]);
        test_actions[i] = malloc(sizeof(int)*test_action_num[i]);
        printf("States: ");
        for (int j = 0; j < test_state_num[i]; j++)
        {
            fscanf(file, "%d", &test_states[i][j]);
            test_states[i][j]--;
        }
        for (int j = 0; j < test_action_num[i]; j++)
        {
            fscanf(file, "%d", &test_actions[i][j]);
            test_actions[i][j]--;
        }

        if (test_types[i][0] == 'e' && test_types[i][1] == 'e')
            pre_e_e_num++;
        if(test_types[i][0] == 'e' && test_types[i][1] == 'a')
            pre_e_a_num++;
        if (test_types[i][0] == 'a' && test_types[i][1] == 'e')
            pre_a_e_num++;
        if (test_types[i][0] == 'a' && test_types[i][1] == 'a')
            pre_a_a_num++;

        printf("Loaded test %d\n", i);
    }
    fclose(file);
    printf("Loaded test parameters\n");

    pre_e_e_times = malloc(sizeof(float)*pre_e_e_num);
    pre_e_a_times = malloc(sizeof(float)*pre_e_a_num);
    pre_a_e_times = malloc(sizeof(float)*pre_a_e_num);
    pre_a_a_times = malloc(sizeof(float)*pre_a_a_num);
    pre_e_e_num = 0;
    pre_a_e_num = 0;
    pre_e_a_num = 0;
    pre_a_a_num = 0;


    clock_t reorder_start = clock();
    Cudd_ReduceHeap(manager, CUDD_REORDER_ANNEALING, 100);
    float reorder_time = (clock() - reorder_time) /(float)CLOCKS_PER_SEC;
    printf("Performing tests\n");
    int test_counter = 0;
    clock_t start_time;
    clock_t end_time;
    float t;
    float tot_time;
    float build_time = 0;
    for (int i = 0; i < number_of_tests; i++)
    {
        if ((100*i / number_of_tests) >= test_counter)
            printf("%%%d\n", test_counter++);
        // printf("Performing test %d\n", i+1);

        start_time = clock();
        DdNode* state_set = makeSet(manager, state_out_vars, var_num,
                                    test_states[i], test_state_num[i], encodings);
        // PRINT_OUT(manager, state_set, a_var_num, var_num, 0);
        DdNode* action_set = makeSet(manager, action_vars, a_var_num,
                                    test_actions[i], test_action_num[i], action_encodings);
        // printf("action_set created\n");
        build_time += (clock() - start_time) /(float)CLOCKS_PER_SEC;

        start_time = clock();

        DdNode* pre_set = pre(manager, state_out_vars, state_in_vars, var_num, action_vars,
                            a_var_num, all_out_states, restored_BDD, state_set, action_set,
                            test_types[i][0], test_types[i][1]);

        end_time = clock();

        t = (float)(end_time - start_time) / (float) CLOCKS_PER_SEC;
        tot_time += t;
        if (test_types[i][0] == 'e' && test_types[i][1] == 'e')
            pre_e_e_times[pre_e_e_num++] = t;
        else if (test_types[i][0] == 'e' && test_types[i][1] == 'a')
            pre_e_a_times[pre_e_a_num++] = t;
        else if (test_types[i][0] == 'a' && test_types[i][1] == 'e')
            pre_a_e_times[pre_a_e_num++] = t;
        else
            pre_a_a_times[pre_a_a_num++] = t;
        // if (Cudd_IsConstant(pre_set))
        // {
        //     printf("Constant pre set = %d\n", Cudd_EquivDC(manager, pre_set, Cudd_ReadLogicZero(manager), Cudd_ReadLogicZero(manager)));
        // }
        Cudd_RecursiveDeref(manager, pre_set);
        Cudd_RecursiveDeref(manager, state_set);
        Cudd_RecursiveDeref(manager, action_set);
    }
    printf("Tests performed\n");

    float tot_e_e_time = 0;
    float tot_e_a_time = 0;
    float tot_a_e_time = 0;
    float tot_a_a_time = 0;
    for (int i = 0; i < pre_e_e_num; i++)
        tot_e_e_time += pre_e_e_times[i];
    for (int i = 0; i < pre_a_e_num; i++)
        tot_a_e_time += pre_a_e_times[i];
    for (int i = 0; i < pre_e_a_num; i++)
        tot_e_a_time += pre_e_a_times[i];
    for (int i = 0; i < pre_a_a_num; i++)
        tot_a_a_time += pre_a_a_times[i];
    printf("Pre times: ee = %f, ea = %f, ae = %f, aa = %f\n",
                        tot_e_e_time, tot_e_a_time, tot_a_e_time, tot_a_a_time);
    printf("Total time: %f\n", tot_time);
    printf("Reordering: %f\n", reorder_time);
    printf("Build time: %f\n", build_time);
    Cudd_DebugCheck(manager);
    Cudd_RecursiveDeref(manager, restored_BDD);
    Cudd_Quit(manager);
    return 0;
}
// construction test
