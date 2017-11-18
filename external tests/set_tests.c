#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include "cudd.h"
#include "dddmp.h"
#include "system_types.h"
#include "set_calcs.h"
#include "DynArray.h"


/**
*    TODO:
*    Test of construction of BDD
*    Test of AND OR operations on trees
*    Test of changing truth assignments
*    Test of calculating pre of things
*    Test of adding variables
*    Test of extracting truth assignments
*/


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
    //fprintf(stderr, "Returning array for : %d, %d\n", num, len);
    for (int i = 0; i < bits; i++)
    {
        array[i] = 1 & (num >> i);
        //fprintf(stderr, "%d", array[i]);
    }
    if (bits < len)
    {
        for (int i = bits; i < len; i++)
        {
            array[i] = 0;
        }
    }
    //fprintf(stderr, "\n");
    return array;
}

// int fromBitCharArray(char* array, int len)
// {
//     int n = 0;
//     for (int i = len-1; i >= 0; i--)
//     {
//         n = (array[i] == '1') | (n << 1);
//     }
//     return n;
// }

// DdNode* makeSet(DdManager* manager, DdNode** vars, int var_num, int* inds, int n, int** encodings)
// {
//   DdNode* set = Cudd_ReadLogicZero(manager);
//   Cudd_Ref(set);
//   // FILE* test = fopen("data.txt", "w");
//   for (int i = 0; i < n; i++)
//   {
//       int* enc = encodings[inds[i]];
//       // printf("Adding encoding for %d: ", inds[i]);
//       // for (int j = 0; j < var_num; j++)
//       // {
//       //     printf("%d ", enc[j]);
//       // }
//       // printf("\n");
//       DdNode* cube = Cudd_bddComputeCube(manager, vars, enc, var_num);
//       Cudd_Ref(cube);
//       DdNode* tmp = Cudd_bddOr(manager, set, cube);
//       Cudd_Ref(tmp);
//       Cudd_RecursiveDeref(manager, cube);
//       Cudd_RecursiveDeref(manager, set);
//       set = tmp;
//   }
//   // fclose(test);
//
//   //Cudd_DebugCheck(manager);
//   return set;
// }

int fromBitCharArray(char* array, int len)
{
    int n = 0;
    for (int i = len-1; i >= 0; i--)
    {
        n = (array[i] == '1') | (n << 1);
    }
    return n;
}

DdNode* makeSet(DdManager* manager, DdNode** vars, int var_num, int* inds, int n, int** encodings)
{
    DdNode* set = Cudd_ReadLogicZero(manager);
    Cudd_Ref(set);
    // FILE* test = fopen("data.txt", "w");
    for (int i = 0; i < n; i++)
    {
        int* enc = encodings[inds[i]];
        // printf("Adding encoding for %d: ", inds[i]);
        // for (int j = 0; j < var_num; j++)
        // {
        //     printf("%d ", enc[j]);
        // }
        // printf("\n");
        DdNode* cube = Cudd_bddComputeCube(manager, vars, enc, var_num);
        Cudd_Ref(cube);
        DdNode* tmp = Cudd_bddOr(manager, set, cube);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, cube);
        Cudd_RecursiveDeref(manager, set);
        set = tmp;
    }
    // fclose(test);

    //Cudd_DebugCheck(manager);
    return set;
}

DdNode* make_state_BDD(DdManager* manager, DdNode** vars, int* encoding, int len)
{
    DdNode *state_BDD = Cudd_ReadOne(manager);
    Cudd_Ref(state_BDD);
    DdNode* tmp;
    for (int i = 0; i < len; i++)
    {
        //printf("Encoded letter: %c\n", encoding[i]);
        if (!encoding[i])
        {
            DdNode* temp_neg = Cudd_Not(vars[i]);
            Cudd_Ref(temp_neg);
            tmp = Cudd_bddAnd(manager, temp_neg, state_BDD);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, state_BDD);
            Cudd_RecursiveDeref(manager, temp_neg);
        }
        else if (encoding[i])
        {
            tmp = Cudd_bddAnd(manager, vars[i], state_BDD);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, state_BDD);
        }
        else
        {
            printf("Error in encoding (1)! %d %d\n", encoding[i], len);
            return NULL;
        }
        state_BDD = tmp;
    }
    return state_BDD;
}

int main (int argc, char *argv[])
{
    if (argc < 5)
    {
        fprintf(stderr, "Too few arguments! Format: [BDD structure file] [BDD transition file] [U groups] [G groups] [test file] [mode={test, time}] [reorder bool] [Number of tests]\n");
        return 1;
    }

    char* BDD_text = argv[2];
    char* save_file_name = argv[1];
    char* test_file = argv[5];
    char* U_file = argv[3];
    char* G_file = argv[4];
    char* exec_mode = argv[6];
    char* reorder_mode = argv[7];
    int tests = atoi(argv[8]);


    // load manager
    DdManager* manager = Cudd_Init(0,0,CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS,0); /* Global BDD manager. */
    Cudd_AutodynEnable(manager, CUDD_REORDER_SYMM_SIFT);

    fprintf(stderr, "Loading BDD\n");
    FILE* save_file = fopen(save_file_name, "r");
    // fprintf(stderr, "file name: %s\n", name);
    DdNode* restored_BDD;
    restored_BDD = Dddmp_cuddBddLoad(manager, DDDMP_VAR_MATCHIDS, NULL, NULL,
                            NULL, DDDMP_MODE_BINARY, save_file_name, save_file);
    fclose(save_file);
    Cudd_Ref(restored_BDD);
    fprintf(stderr, "BDD loaded : %d vars\n", Cudd_ReadSize(manager));
    // Cudd_ReduceHeap(manager, CUDD_REORDER_SYMM_SIFT, 100);

    int** encodings;
    int* enc_state_map;
    int** action_encodings;
    int var_num, a_var_num, state_num, action_num;
    DdNode** state_in_vars;
    DdNode** state_out_vars;
    DdNode** action_vars;
    fprintf(stderr, "Loading BDD parameters\n");
    FILE* trans_file = fopen(BDD_text, "r");
    // fprintf(stderr, "Parameter file loaded\n");
    fscanf(trans_file, "%*c %d %d\n", &state_num, &var_num);
    state_num++;

    enc_state_map = malloc(sizeof(int)*pow(2, var_num));
    memset(enc_state_map, '2', pow(2, var_num));
    encodings = malloc(sizeof(int*)*state_num);
    int not_unique = 0;
    int* visited = calloc(sizeof(int), pow(2, var_num));
    for (int i = 0; i < state_num; i++)
    {
        int id;
        char* enc = malloc(sizeof(char)*var_num);
        encodings[i] = malloc(sizeof(int*)*var_num);
        fscanf(trans_file, "%d %s\n", &id, enc);
        // fprintf(stderr, "Read %d %s\n", id, enc);
        int map_id = fromBitCharArray(enc, var_num);
        if (visited[map_id])
        {
            fprintf(stderr, "Duplicate encodings: %d=%d\n", enc_state_map[map_id]+1, i);
            not_unique = 1;
        }
        else
            visited[map_id] = 1;
        enc_state_map[map_id] = id-1;
        for (int k = 0; k < var_num; k++)
        {
            encodings[i][k] = enc[k] == '1';
        }
        free(enc);
    }
    if (not_unique)
    {
        fprintf(stderr, "Encoding not unique!\n");
        return 0;
    }
    // fprintf(stderr, "Read params: state_num=%d var_num=%d\n", state_num, var_num);
    // fprintf(stderr, "Read state and action parameters\n");
    state_in_vars = malloc(sizeof(DdNode*)*var_num);
    state_out_vars = malloc(sizeof(DdNode*)*var_num);
    // fprintf(stderr, "Allocated space\n");

    fprintf(stderr, "Encodings loaded\n");

    int trans_num;
    fscanf(trans_file, "%*s %d %d\n", &action_num, &trans_num);

    a_var_num = ceil(log10(action_num)/log10(2));
    action_encodings = malloc(sizeof(int*)*action_num);
    for (int i = 0; i < action_num; i++)
        action_encodings[i] = toBitArray(i, a_var_num);
    action_vars = malloc(sizeof(DdNode*)*a_var_num);

    for (int i = 0; i < a_var_num; i++)
    {
        action_vars[i] = Cudd_bddIthVar(manager, i);
    }
    for (int i = 0; i < var_num; i++)
    {
        state_in_vars[i] = Cudd_bddIthVar(manager, i+a_var_num);
    }
    for (int i = 0; i < var_num; i++)
    {
        state_out_vars[i] = Cudd_bddIthVar(manager, i+a_var_num + var_num);
    }

        fprintf(stderr, "BDD parameters loaded : %d vars\n", Cudd_ReadSize(manager));

    // restored_BDD = Cudd_ReadLogicZero(manager);
    // Cudd_Ref(restored_BDD);
    for (int i = 0; i < trans_num; i++)
    {
        int start;
        int end;
        char* a_enc;
        // int* int_enc = malloc(sizeof(int)*a_var_num);
        fscanf(trans_file, "%*d %*s %d\n", &end);
        // for (int j = 0; j < a_var_num; j++)
        // {
        //     int_enc[i] = a_enc[j] == '1';
        // }
        // DdNode* start_state = make_state_BDD(manager, state_in_vars, encodings[start-1], var_num);
        // DdNode* end_state = make_state_BDD(manager, state_out_vars, encodings[end-1], var_num);
        // DdNode* action_state = make_state_BDD(manager, action_vars, int_enc, a_var_num);
        //
        // DdNode* trans_state = Cudd_bddAnd(manager, start_state, action_state);
        // Cudd_Ref(trans_state);
        // Cudd_RecursiveDeref(manager, start_state);
        // Cudd_RecursiveDeref(manager, end_state);
        // DdNode* tmp = Cudd_bddAnd(manager, trans_state, end_state);
        // Cudd_Ref(tmp);
        // Cudd_RecursiveDeref(manager, trans_state);
        // Cudd_RecursiveDeref(manager, end_state);
        // trans_state = tmp;
        //
        // tmp = Cudd_bddOr(manager, restored_BDD, trans_state);
        // Cudd_Ref(tmp);
        // Cudd_RecursiveDeref(manager, restored_BDD);
        // Cudd_RecursiveDeref(manager, trans_state);
        // restored_BDD = tmp;
    }

    fprintf(stderr, "Transitions skipped\n");
    // fprintf("BDD created from file\n");

    // int pg_num;
    // fscanf(trans_file, "%*s %d\n", &pg_num);
    // int** pg_U = malloc(sizeof(int*)*pg_num);
    // int** pg_G = malloc(sizeof(int*)*pg_num);
    //
    // for (int i = 0; i < pg_num; i++)
    // {
    //     int U_num;
    //     int G_num;
    //     fscanf(trans_file, "%*s %d %d\n", &U_num, &G_num);
    //
    //     pg_U[i] = malloc(sizeof(int)*U_num);
    //     for (int j = 0; j < U_num; j++)
    //     {
    //         fscanf(trans_file, "%d", &pg_U[i][j]-1);
    //     }
    //     pg_G[i] = malloc(sizeof(int)*G_num);
    //     for (int j = 0; j < G_num; j++)
    //     {
    //         fscanf(trans_file, "%d", &pg_G[i][j]-1);
    //     }
    // }
    // fclose(trans_file);




    fprintf(stderr, "Loading progress groups: %s\n", U_file);
    save_file = fopen(U_file, "r");
    DdNode** U_groups;
    int prog_num = Dddmp_cuddBddArrayLoad(manager, DDDMP_ROOT_MATCHLIST, NULL,
                            DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_BINARY, NULL, save_file, &U_groups);

    // for (int i = 0; i < prog_num; i++)
    //     PRINT_A(manager, U_groups[i], a_var_num, var_num, 1);

    fprintf(stderr, "Loading progress groups: %s\n", G_file);
    save_file = fopen(G_file, "r");
    DdNode** G_groups;

    prog_num = Dddmp_cuddBddArrayLoad(manager, DDDMP_ROOT_MATCHLIST, NULL,
                            DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_BINARY, NULL, save_file, &G_groups);

    fprintf(stderr, "Creating complete state set\n");
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

    int* states = calloc(sizeof(int), state_num);
    for (int i = 0; i < state_num; i++)
        states[i] = i;
    all_out_states = makeSet(manager, state_out_vars, var_num, states, state_num, encodings);
    free(states);
    states = calloc(sizeof(int), action_num);
    for (int i = 0; i < action_num; i++)
        states[i] = i;
    DdNode* all_actions = makeSet(manager, action_vars, a_var_num, states, action_num, action_encodings);
    fprintf(stderr, "Set created\n");

    BDDSys system;
    system.manager = manager;
    system.trans_sys = restored_BDD;
    array_init(s_in_vars, BDDlist, var_num, 1);
    for (int i = 0; i < var_num; i++)
        array_set(s_in_vars, i, state_in_vars[i]);
    system.s_in_vars = s_in_vars;
    array_init(s_out_vars, BDDlist, var_num, 1);
    for (int i = 0; i < var_num; i++)
        array_set(s_out_vars, i, state_out_vars[i]);
    system.s_out_vars = s_out_vars;
    array_init(a_vars, BDDlist, a_var_num, 1);
    for (int i = 0; i < a_var_num; i++)
        array_set(a_vars, i, action_vars[i]);
    system.a_vars = a_vars;
    system.s_var_num = var_num;
    system.a_var_num = a_var_num;
    system.all_states = all_out_states;
    system.all_actions = all_actions;
    // system.enc_state_map = enc_state_map;
    array_init(pg_U, BDDlist, prog_num, 1);
    for (int i = 0; i < prog_num; i++)
        array_set(pg_U, i, U_groups[i]);
    system.pg_U = pg_U;
    array_init(pg_G, BDDlist, prog_num, 1);
    for (int i = 0; i < prog_num; i++)
        array_set(pg_G, i, G_groups[i]);
    system.pg_G = pg_G;
    system.pg_num = prog_num;

    char input[100];
    if (strcmp(exec_mode, "test") == 0)
    {
        // test section : win_primal
        FILE* file = fopen(test_file, "r");
        input[0] = 0;
        fscanf(file, "%s", input);
        fprintf(stderr, "Read start: %s\n", input);
        int counter = 1;
        while (strcmp(input, "test") == 0)
        {
            int A_vec_len = 0;
            fscanf(file, "%*s %d", &A_vec_len);
            int* A_vec = malloc(sizeof(int)*A_vec_len);
            for (int i = 0; i < A_vec_len; i++)
            {
                fscanf(file, "%d", &A_vec[i]);
                A_vec[i]--;
            }
            DdNode* A = makeSet(manager, state_out_vars, var_num, A_vec, A_vec_len, encodings);
            free(A_vec);
            // read B vector
            int B_vec_len = 0;
            fscanf(file, "%*s %d", &B_vec_len);
            int* B_vec = malloc(sizeof(int)*B_vec_len);
            for (int i = 0; i < B_vec_len; i++)
            {
                fscanf(file, "%d", &B_vec[i]);
                B_vec[i]--;
            }
            DdNode* B = makeSet(manager, state_out_vars, var_num, B_vec, B_vec_len, encodings);
            free(B_vec);

            int C_list_len = 0;
            fscanf(file, "%*s %d\n", &C_list_len);
            int** C_vec = malloc(sizeof(int*)*C_list_len);
            DdNode** C_list = malloc(sizeof(DdNode*)*C_list_len);
            for (int i = 0; i < C_list_len; i++)
            {
                int C_vec_len = 0;
                fscanf(file, "%d", &C_vec_len);
                C_vec[i] = malloc(sizeof(int)*C_vec_len);
                for (int j = 0; j < C_vec_len; j++)
                {
                    fscanf(file, "%d", &C_vec[i][j]);
                    C_vec[i][j]--;
                }
                C_list[i] = makeSet(manager, state_out_vars, var_num, C_vec[i], C_vec_len, encodings);
                free(C_vec[i]);
            }
            free(C_vec);

            char quant1;
            fscanf(file, "%*s %c", &quant1);
            char quant2;
            fscanf(file, "%*s %c", &quant2);

            // read answer
            int ans_vec_len = 0;
            fscanf(file, "%*s %d", &ans_vec_len);
            int* ans_vec = malloc(sizeof(int)*ans_vec_len);
            for (int i = 0; i < ans_vec_len; i++)
            {
                fscanf(file, "%d", &ans_vec[i]);
                ans_vec[i]--;
            }
            DdNode* ans = makeSet(manager, state_out_vars, var_num, ans_vec, ans_vec_len, encodings);
            free(ans_vec);

            int mode;
            fscanf(file, "%*s %d\n", &mode);

            DdNode* cand;
            if (mode >= WIN_CANDIDATE_SET)
            {
                // read candidate
                int cand_vec_len = 0;
                fscanf(file, "%*s %d", &cand_vec_len);
                int* cand_vec = malloc(sizeof(int)*cand_vec_len);
                for (int i = 0; i < cand_vec_len; i++)
                {
                    fscanf(file, "%d", &cand_vec[i]);
                    cand_vec[i]--;
                }
                cand = makeSet(manager, state_out_vars, var_num, cand_vec, cand_vec_len, encodings);
                free(cand_vec);
            }

            // calculate inv
            DdNode* res;
            DdNode* cand_res;
            DdNode** in = (DdNode**)win_primal(&system, A, B, C_list, C_list_len,
                              quant1, quant2, NULL, mode);
            if (mode >= WIN_SET)
            {
                res = in[0];
                if (mode >= WIN_CANDIDATE_SET)
                    cand_res = in[1];
            }
            free(in);

            // compare with ans
            fprintf(stderr, "win_primal test %d ", counter);
            if (Cudd_EquivDC(manager, ans, res, Cudd_ReadLogicZero(manager))
                && (mode == WIN_SET || Cudd_EquivDC(manager, cand, cand_res, Cudd_ReadLogicZero(manager))))
                fprintf(stderr, "passed\n");
            else
            {
                fprintf(stderr, "failed\n");
                if (!Cudd_EquivDC(manager, ans, res, Cudd_ReadLogicZero(manager)))
                    fprintf(stderr, "answer did not match\n");
                if (mode >= WIN_CANDIDATE_SET && !Cudd_EquivDC(manager, cand, cand_res, Cudd_ReadLogicZero(manager)))
                {
                    fprintf(stderr, "candidate did not match\n");
                    fprintf(stderr, "Candidate:\n");
                    // PRINT_OUT(manager, cand, a_var_num, var_num, enc_state_map, 0, NULL);
                    fprintf(stderr, "Result\n");
                    // PRINT_OUT(manager, cand_res, a_var_num, var_num, enc_state_map, 0, NULL);
                }
                return 0;
            }
            counter++;
            input[0] = 0;
            int num = fscanf(file, "%s", input);
        }
        fclose(file);
    }

    if (strcmp(exec_mode, "time") == 0)
    {
        // test section : win_primal
        FILE* file = fopen(test_file, "r");
        input[0] = '\0';
        fscanf(file, "%s", input);
        fprintf(stderr, "Read start: %s\n", input);
        int counter = 1;
        while (strcmp(input, "test") == 0)
        {
            // read A vector
            int A_vec_len = 0;
            fscanf(file, "%*s %d", &A_vec_len);
            int* A_vec = malloc(sizeof(int)*A_vec_len);
            for (int i = 0; i < A_vec_len; i++)
            {
                fscanf(file, "%d", &A_vec[i]);
                A_vec[i]--;
            }
            DdNode* A = makeSet(manager, state_out_vars, var_num, A_vec, A_vec_len, encodings);
            free(A_vec);

            // read B vector
            int B_vec_len = 0;
            fscanf(file, "%*s %d", &B_vec_len);
            int* B_vec = malloc(sizeof(int)*B_vec_len);
            for (int i = 0; i < B_vec_len; i++)
            {
                fscanf(file, "%d", &B_vec[i]);
                B_vec[i]--;
            }
            DdNode* B = makeSet(manager, state_out_vars, var_num, B_vec, B_vec_len, encodings);
            free(B_vec);

            int C_list_len = 0;
            fscanf(file, "%*s %d\n", &C_list_len);
            int** C_vec = malloc(sizeof(int*)*C_list_len);
            DdNode** C_list = malloc(sizeof(DdNode*)*C_list_len);
            for (int i = 0; i < C_list_len; i++)
            {
                int C_vec_len = 0;
                fscanf(file, "%d", &C_vec_len);
                C_vec[i] = malloc(sizeof(int)*C_vec_len);
                for (int j = 0; j < C_vec_len; j++)
                {
                    fscanf(file, "%d", &C_vec[i][j]);
                    C_vec[i][j]--;
                }
                C_list[i] = makeSet(manager, state_out_vars, var_num, C_vec[i], C_vec_len, encodings);
                free(C_vec[i]);
            }
            free(C_vec);

            char quant1;
            fscanf(file, "%*s %c", &quant1);
            char quant2;
            fscanf(file, "%*s %c", &quant2);

            // read answer
            int ans_vec_len = 0;
            fscanf(file, "%*s %d", &ans_vec_len);
            int* ans_vec = malloc(sizeof(int)*ans_vec_len);
            for (int i = 0; i < ans_vec_len; i++)
            {
                fscanf(file, "%d", &ans_vec[i]);
                ans_vec[i]--;
            }
            DdNode* ans = makeSet(manager, state_out_vars, var_num, ans_vec, ans_vec_len, encodings);
            free(ans_vec);

            int mode;
            fscanf(file, "%*s %d\n", &mode);

            DdNode* cand;
            if (mode >= WIN_CANDIDATE_SET)
            {
                // read candidate
                int cand_vec_len = 0;
                fscanf(file, "%*s %d", &cand_vec_len);
                int* cand_vec = malloc(sizeof(int)*cand_vec_len);
                for (int i = 0; i < cand_vec_len; i++)
                {
                    fscanf(file, "%d", &cand_vec[i]);
                    cand_vec[i]--;
                }
                cand = makeSet(manager, state_out_vars, var_num, cand_vec, cand_vec_len, encodings);
                free(cand_vec);
            }

            //reorder option
            if (*reorder_mode == '1')
            {
                clock_t reorder_time = clock();
                Cudd_ReduceHeap(manager, CUDD_REORDER_ANNEALING, 500);
                double reorder_secs = (double)(clock() - reorder_time) / (double)CLOCKS_PER_SEC;
                fprintf(stderr, "Reorder time: %f seconds\n", reorder_secs);
            }

            // calculate inv
            DdNode* res;
            DdNode* cand_res;
            clock_t time = 0;
            clock_t start = clock();
            for (int k = 0; k < tests; k++)
            {
                //start = clock();
                DdNode** in = (DdNode**)win_primal(&system, A, B, C_list, C_list_len,
                                  quant1, quant2, NULL, mode);
                for (int j = 0; j < mode; j++)
                {
                    Cudd_RecursiveDeref(manager, in[j]);
                }
            }
            double seconds = (double)(clock() - start)/(double)(CLOCKS_PER_SEC*tests);
            fprintf(stderr, "Win_primal %d time: %f seconds\n", counter, seconds);
            printf("%.10f", seconds);

            Cudd_RecursiveDeref(manager, A);
            Cudd_RecursiveDeref(manager, B);
            for (int i = 0; i < C_list_len; i++)
                Cudd_RecursiveDeref(manager, C_list[i]);
            free(C_list);

            counter++;
            input[0] = 0;
            int num = fscanf(file, "%s", input);
        }
        fclose(file);
    }

    Cudd_RecursiveDeref(manager, restored_BDD);
    Cudd_Quit(manager);
    return 0;
}
// construction test
