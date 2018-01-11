#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include "cudd.h"
#include "dddmp.h"
#include "system_types.h"
#include "set_calcs.h"
#include "DynArray.h"

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

    for (int i = 0; i < n; i++)
    {
        int* enc = encodings[inds[i]];
        DdNode* cube = Cudd_bddComputeCube(manager, vars, enc, var_num);
        Cudd_Ref(cube);
        DdNode* tmp = Cudd_bddOr(manager, set, cube);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, cube);
        Cudd_RecursiveDeref(manager, set);
        set = tmp;
    }
    return set;
}

/*
    This script is capable of performing comparison tests and performance tests
    of the winning set calculation functions for BDDs

    The test loads pre-computed BDDs for transition systems and progress groups
    and loads the encodings of states other critical information to the test
    from an accompanying transition file on the format

    S [number of states, excluding "outside state"] [number of state variables]
    [state 1] [encoding 1]
    [state 2] [encoding 2]
    ...
    ...
    [last state] [last encoding]
    T [number of actions] [number of transitions]
    [start state 1] [action encoding] [end state 1]
    [start state 2] [action encoding] [end state 2]
    ...
    ...
    [last start state] [action encoding] [last end state]
    (any test after this is irrelevant)

    The BDD files should be files generated using dddmp functions for saving
    BDDs

    It then takes a file containing test information in the format

    test
    A [number of states] [state 1] [state 2] ... [last state in A]
    B [number of states] [state 1] [state 2] ... [last state in B]
    C [number of sets]
    [size of set 1] [state 1] [state 2] ... [last state in set 1]
    ...
    [size of last set] [state 1] [state 2] ... [last state in last set]
    q1 [quantifier 1]
    q2 [quantifier 2]
    ans [length of answer winning set] [state 1] ... [last state in winning set]
    mode [mode of test (winning=1 or candidate set=2)]
    (if mode is 2)
    cand [length of answer candidate set] [state 1] ... [last state in candidate set]
    (repeat if more than one test is present)

    Output if time test:
    [size of BDD transition system (node count)] [average run time of test]
*/
int main (int argc, char *argv[])
{
    if (argc < 5)
    {
        fprintf(stderr, "Too few arguments! Format: [BDD structure file] [BDD transition file] [U groups] [G groups] [test file] [mode={test, time}] [reorder bool] [Number of tests]\n");
        return 1;
    }

    // Load test arguments
    char* BDD_text = argv[2];
    char* save_file_name = argv[1];
    char* test_file = argv[5];
    char* U_file = argv[3];
    char* G_file = argv[4];
    char* exec_mode = argv[6];
    char* reorder_mode = argv[7];
    int tests = atoi(argv[8]);

    // load CUDD BDD node manager
    DdManager* manager = Cudd_Init(0,0,CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS,0); /* Global BDD manager. */
    Cudd_AutodynEnable(manager, CUDD_REORDER_SYMM_SIFT);

    // Load pre-computed transition for test
    FILE* save_file = fopen(save_file_name, "r");
    DdNode* restored_BDD;
    // load BDD for transition system
    restored_BDD = Dddmp_cuddBddLoad(manager, DDDMP_VAR_MATCHIDS, NULL, NULL,
                            NULL, DDDMP_MODE_BINARY, save_file_name, save_file);
    fclose(save_file);
    Cudd_Ref(restored_BDD);
    fprintf(stderr, "Transition BDD loaded, variable count: %d vars, node size: %d\n", Cudd_ReadSize(manager), Cudd_DagSize(restored_BDD));
    // give size of BDD transition system as output, unless reordering is performed later
    if (*reorder_mode != '1')
        printf("%d ", Cudd_DagSize(restored_BDD));

    // initialize variables for encodings, inverse encoding map and variable
    // parameters
    int** encodings;
    int* enc_state_map;
    int** action_encodings;
    int var_num, a_var_num, state_num, action_num;
    DdNode** state_in_vars;
    DdNode** state_out_vars;
    DdNode** action_vars;

    // Open transition file, to load encodings and variable parameters
    FILE* trans_file = fopen(BDD_text, "r");
    fscanf(trans_file, "%*c %d %d\n", &state_num, &var_num);
    state_num++;

    // Initialize inverse encoding map
    enc_state_map = malloc(sizeof(int)*pow(2, var_num));
    memset(enc_state_map, '2', pow(2, var_num));
    encodings = malloc(sizeof(int*)*state_num);
    int not_unique = 0;
    // list of mapped states, to check for injective inverse map
    int* visited = calloc(sizeof(int), pow(2, var_num));

    // Load state encodings and make inverse map from transition file
    for (int i = 0; i < state_num; i++)
    {
        int id;
        char* enc = malloc(sizeof(char)*var_num);
        encodings[i] = malloc(sizeof(int*)*var_num);
        fscanf(trans_file, "%d %s\n", &id, enc); // read encoding for state "id"

        int map_id = fromBitCharArray(enc, var_num);
        if (visited[map_id]) // map not injective
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
        // map was found not to be injective, terminates test
        fprintf(stderr, "Encoding not unique!\n");
        return 0;
    }


    // Load action and transition parameters
    int trans_num;
    fscanf(trans_file, "%*s %d %d\n", &action_num, &trans_num);
    a_var_num = ceil(log10(action_num)/log10(2));

    // Create action encoding, assumes the logarithmic encoding for actions
    action_encodings = malloc(sizeof(int*)*action_num);
    for (int i = 0; i < action_num; i++)
        action_encodings[i] = toBitArray(i, a_var_num);

    fprintf(stderr, "Encodings loaded\n");

    // Initialize variable BDDs
    action_vars = malloc(sizeof(DdNode*)*a_var_num);
    state_in_vars = malloc(sizeof(DdNode*)*var_num);
    state_out_vars = malloc(sizeof(DdNode*)*var_num);

    // Create variable BDDs
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

    // Skip the transition part of the file
    // Not needed as tests load pre-constructed BDDs made from the same
    // transition file
    for (int i = 0; i < trans_num; i++)
    {
        int start;
        int end;
        char* a_enc;
        fscanf(trans_file, "%*d %*s %d\n", &end);
    }

    // Load pre-constructed progress group BDDs
    // action groups
    save_file = fopen(U_file, "r");
    DdNode** U_groups;
    int prog_num = Dddmp_cuddBddArrayLoad(manager, DDDMP_ROOT_MATCHLIST, NULL,
                            DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_BINARY, NULL, save_file, &U_groups);

    // state groups
    save_file = fopen(G_file, "r");
    DdNode** G_groups;
    prog_num = Dddmp_cuddBddArrayLoad(manager, DDDMP_ROOT_MATCHLIST, NULL,
                            DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_BINARY, NULL, save_file, &G_groups);

    fprintf(stderr, "Progress group BDDs loaded\n");

    // Construct BDD for the state set (as end states of transitions)
    DdNode* all_out_states;
    int* states = calloc(sizeof(int), state_num);
    for (int i = 0; i < state_num; i++)
        states[i] = i;
    all_out_states = makeSet(manager, state_out_vars, var_num, states, state_num, encodings);
    free(states);
    // Construct BDD for the action set
    states = calloc(sizeof(int), action_num);
    for (int i = 0; i < action_num; i++)
        states[i] = i;
    DdNode* all_actions = makeSet(manager, action_vars, a_var_num, states, action_num, action_encodings);
    free(states);

    fprintf(stderr, "State and action sets created\n");

    // Enclose all relevant variables into a BDDSys struct and dynamic arrays
    // to pass to winning set calculation functions
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
    array_init(pg_U, BDDlist, prog_num, 1);
    for (int i = 0; i < prog_num; i++)
        array_set(pg_U, i, U_groups[i]);
    system.pg_U = pg_U;
    array_init(pg_G, BDDlist, prog_num, 1);
    for (int i = 0; i < prog_num; i++)
        array_set(pg_G, i, G_groups[i]);
    system.pg_G = pg_G;
    system.pg_num = prog_num;

    char input[100]; // test file input buffer
    if (strcmp(exec_mode, "test") == 0) // Perform a result comparison test
    {
        // Load test file parameters
        FILE* file = fopen(test_file, "r");
        input[0] = 0;
        fscanf(file, "%s", input); // read start line of test
        int counter = 1;
        while (strcmp(input, "test") == 0) // As long as there is a test to run
        {
            // Load the A set and make BDD for win_primal test arguments
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

            // Load the B set and make BDD for win_primal test arguments
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

            // Load the C sets and make list BDDs for win_primal test arguments
            int C_list_len = 0;
            fscanf(file, "%*s %d\n", &C_list_len);
            int** C_vec = malloc(sizeof(int*)*C_list_len);
            DdNode** C_list = malloc(sizeof(DdNode*)*C_list_len);
            for (int i = 0; i < C_list_len; i++) // for each C set
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

            // Load test quantifiers
            char quant1;
            fscanf(file, "%*s %c", &quant1);
            char quant2;
            fscanf(file, "%*s %c", &quant2);

            // Load answer winning set
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

            // Read win_primal mode (winning set only or win and candidate set)
            int mode;
            fscanf(file, "%*s %d\n", &mode);

            DdNode* cand;
            if (mode >= WIN_CANDIDATE_SET) // If test includes a candidate set
            {
                // Load candidate set and create corresponding BDD
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

            // Calculate the winning (and candidate set)
            DdNode* res;
            DdNode* cand_res;
            DdNode** in = (DdNode**)win_primal(&system, A, B, C_list, C_list_len,
                              quant1, quant2, NULL, mode);
            if (mode >= WIN_SET) // load result
            {
                res = in[0];
                if (mode >= WIN_CANDIDATE_SET)
                    cand_res = in[1];
            }
            free(in);

            // Compare with answer
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
                }
                return 0;
            }
            counter++;
            input[0] = 0;
            int num = fscanf(file, "%s", input);
        }
        fclose(file);
    }

    if (strcmp(exec_mode, "time") == 0)  // Perform timing test
    {
        // Load test file parameters
        FILE* file = fopen(test_file, "r");
        input[0] = '\0';
        fscanf(file, "%s", input);
        int counter = 1;
        while (strcmp(input, "test") == 0) // read start line of test
        {
            // Load the A set and make BDD for win_primal test arguments
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

            // Load the B set and make BDD for win_primal test arguments
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

            // Load the C sets and make list BDDs for win_primal test arguments
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

            // Load test quantifiers
            char quant1;
            fscanf(file, "%*s %c", &quant1);
            char quant2;
            fscanf(file, "%*s %c", &quant2);

            // Read answer winning set (not needed for test, but to move stream pointer)
            int ans_vec_len = 0;
            fscanf(file, "%*s %d", &ans_vec_len);
            int* ans_vec = malloc(sizeof(int)*ans_vec_len);
            for (int i = 0; i < ans_vec_len; i++)
            {
                fscanf(file, "%d", &ans_vec[i]);
            }
            free(ans_vec);

            // Read win_primal mode (winning set only or win and candidate set)
            int mode;
            fscanf(file, "%*s %d\n", &mode);

            // Check the reorder option of the test
            if (*reorder_mode == '1') // perform reordering of BDDs
            {
                clock_t reorder_time = clock();
                // Reorder with annealing
                Cudd_ReduceHeap(manager, CUDD_REORDER_ANNEALING, 500);
                double reorder_secs = (double)(clock() - reorder_time) / (double)CLOCKS_PER_SEC;
                fprintf(stderr, "Reordering performed, time: %f seconds\n", reorder_secs);
                printf("%d ", Cudd_DagSize(restored_BDD));
            }

            // Performance test
            fprintf(stderr, "Performing test\n");
            clock_t time = 0;
            clock_t start = clock();
            // Perform many tests and take average (needed since clock lacks precision)
            for (int k = 0; k < tests; k++)
            {
                DdNode** in = (DdNode**)win_primal(&system, A, B, C_list, C_list_len,
                                  quant1, quant2, NULL, mode);
                // have to release node memory here to not cause other performance issues
                // will impact the performance test, but minimally
                for (int j = 0; j < mode; j++)
                {
                    Cudd_RecursiveDeref(manager, in[j]);
                }
                free(in);
            }
            double seconds = (double)(clock() - start)/(double)(CLOCKS_PER_SEC*tests);
            fprintf(stderr, "Win_primal %d time: %f seconds\n", counter, seconds);
            printf("%.10f", seconds); // output time to script

            // Release node memory for this test
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
