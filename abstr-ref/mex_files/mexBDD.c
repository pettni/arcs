#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdint.h>
#include "cudd.h"
#include "mex.h"
#include "matrix.h"
#include "system_types.h"
#include "set_calcs.h"
#include "DynArray.h"

#define MAX(x1, x2) x1 < x2 ? x2 : x1

void free_all();
void free_system(BDDSys* system_to_free);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
BDDSys* initializeBDD(uint var_num, EncList* _s_encs, uint a_var_num, EncList* _a_encs);
// pads the uint list with the amount of zeros specified by spaces
void zeropad(NumList* list, uint pos, uint pad_len);
void add_action(BDDSys* sys, uint index, NumList* enc, uint new_var_num);
void add_state(BDDSys* sys, uint index, NumList* enc, uint new_var_num);
void set_state_enc(BDDSys* sys, uint ind, NumList* enc);
void add_trans(BDDSys* sys, NumList* in_state, NumList* action, NumList* out_state);
DdNode* get_trans_with_s(BDDSys* sys, NumList* states);
void rm_trans_with_s(BDDSys* sys, NumList* state);
void add_progress_group(BDDSys* sys, NumList* U, NumList* G);
void rm_progress_group(BDDSys* sys, uint index);
int has_superior_pg(BDDSys* sys, NumList* U, NumList* G);
void add_to_pg(BDDSys* sys, NumList* indices, NumList* states);
NumList* is_member_of_pg(BDDSys* sys, NumList* states);
uint get_num_of_trans(BDDSys* sys, DdNode* bdd);
uint** read_transitions(BDDSys* sys, DdNode* T);
uint** read_actions(BDDSys* sys, DdNode* T);
uint** read_out_states(BDDSys* sys, DdNode* T);
uint** read_in_states(BDDSys* sys, DdNode* T);
void read_in(NumList** intervals, uint interv_num, NumList* cube_inds, int* cube, uint cube_pos, uint** outputs, uint* output_num);
uint fromBitArray(uint* arr, uint len);
uint fromBitArray_int(int* arr, uint len);
DdNode* makeSet(DdManager* manager, DdNode** vars, int var_num, int* inds, int n, EncList* encodings);
void write_BDD_to_mxarray(BDDSys* sys, DdNode* bdd, mxArray** array, char setting);
DdNode* read_BDD_from_mxarray(BDDSys* sys, const mxArray* array, char setting);

BDDSysList* allocated_BDDs;
const uint BDD_alloc_buffer = 10;
const uint var_buffer = 10;
const uint s_buffer = 500;
const uint pg_buffer = 10;
int initialized = 0;
Cudd_ReorderingType reorder_alg = CUDD_REORDER_ANNEALING;

// DdManager*  manager;
// DdNode*     trans_sys;
// DdNode**    s_in_vars;
// DdNode**    s_out_vars;
// DdNode**    a_vars;
// uint        s_var_num;
// uint        a_var_num;
// DdNode*     all_states;
// DdNode*     all_actions;
// uint**      s_encs;
// uint**      a_encs;
// uint        s_encs_len;
// uint        a_encs_len;
// uint*       enc_state_map;
// DdNode**    pg_U;
// DdNode**    pg_G;
// uint        pg_num;

// int initialized = 0;


/*
     BDD commands:
     initialize       -- Initialize BDD and needed properties
     add_a            -- add action to system
     add_s            -- add state to system
     add_trans        -- add transition to system
     get_trans_with_s -- get all transitions from or to state
     rm_trans_with_s  -- remove all transition from or to state
     rm_pg            -- remove progress group
     add_pg           -- add progress

     pre              -- get pre of given action and state set
     win_primal       -- calculate winning/candidate set

     General conventions:
     All BDDs given as return arguments are references before returned
     and dereferenced outside function
*/

/*
     Destructor for when the mex file is cleared or MATLAB shuts down.
     Releases all memory allocated for created BDD systems.
*/
void free_all()
{
     mexPrintf("Entering free\n");
     for (int i = 0; i < array_len(allocated_BDDs); i++)
     {
          mexPrintf("Clearing system %d\n", i);
          BDDSys* sys = array_get(allocated_BDDs, i);
          if (sys != NULL)
               free_system(sys);
          mexPrintf("System cleared\n");
     }
     mexPrintf("Freeing system list\n");
     array_free(allocated_BDDs);
     mexPrintf("Memory freed\n");
}
void free_system(BDDSys* system_to_free)
{
     if (initialized && system_to_free != NULL)
     {
          mexPrintf("Clearing transition system\n");
          Cudd_RecursiveDeref(system_to_free->manager, system_to_free->trans_sys);
          mexPrintf("Clearing state BDD\n");
          Cudd_RecursiveDeref(system_to_free->manager, system_to_free->all_states);
          mexPrintf("Clearing action BDD\n");
          Cudd_RecursiveDeref(system_to_free->manager, system_to_free->all_actions);
          mexPrintf("Clearing s_in_vars\n");
          array_free(system_to_free->s_in_vars);
          array_free(system_to_free->s_in_inds);
          mexPrintf("Clearing s_out_vars\n");
          array_free(system_to_free->s_out_vars);
          array_free(system_to_free->s_out_inds);
          mexPrintf("Clearing a_vars\n");
          array_free(system_to_free->a_vars);
          array_free(system_to_free->a_inds);
          mexPrintf("Clearing s_encodings\n");
          for (int i = 0; i < array_len(system_to_free->s_encs); i++)
          {
               mexPrintf("Clearing s_encoding %d\n", i);
               mexEvalString("drawnow;");
               if (array_get(system_to_free->s_encs, i) != NULL)
                    array_free(array_get(system_to_free->s_encs, i));
          }
          mexPrintf("Clearing s_encoding list\n");
          array_free(system_to_free->s_encs);

          mexPrintf("Clearing a_encodings\n");
          for (int i = 0; i < array_len(system_to_free->a_encs); i++)
          {
               mexPrintf("Clearing a_encoding %d\n", i);
               if (array_get(system_to_free->a_encs, i) != NULL)
                    array_free(array_get(system_to_free->a_encs, i));
          }
          mexPrintf("Clearing a_encoding list\n");
          array_free(system_to_free->a_encs);
          mexPrintf("Clearing enc_state_map\n");
          array_free(system_to_free->enc_state_map);
          mexPrintf("Clearing pg_U\n");
          for (int i = 0; i < array_len(system_to_free->pg_U); i++)
          {
               mexPrintf("Clearing U group %d\n", i);
               if (array_get(system_to_free->pg_U, i) != NULL)
                    Cudd_RecursiveDeref(system_to_free->manager, array_get(system_to_free->pg_U, i));
          }
          array_free(system_to_free->pg_U);
          mexPrintf("Clearing pg_G\n");
          for (int i = 0; i < array_len(system_to_free->pg_G); i++)
          {
               mexPrintf("Clearing G group %d\n", i);
               if (array_get(system_to_free->pg_G, i) != NULL)
                    Cudd_RecursiveDeref(system_to_free->manager, array_get(system_to_free->pg_G, i));
          }
          array_free(system_to_free->pg_G);
          mexPrintf("Checking non-zero references: %d\n", Cudd_CheckZeroRef(system_to_free->manager));
          mexPrintf("Quitting CUDD\n");
          Cudd_Quit(system_to_free->manager);
     }
}


/*
     Gateway function for MATLAB
     Index convention: All MATLAB indices start at 1, all C indices start at 0
     All indices are decremented here before passed to outer functions
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
     if (!initialized)
     {
          array_init(list, BDDSysList, 0, BDD_alloc_buffer);
          allocated_BDDs = list;
          array_make_persist(allocated_BDDs);
          for (int i = 0; i < array_len(allocated_BDDs); i++)
               array_set(allocated_BDDs, i, NULL);
          mexAtExit(&free_all);
          initialized = 1;
     }
     if (nrhs == 0)
     {
          // TODO: argument error
          return;
     }

     // interpret command
     char* command = mxArrayToString(prhs[0]);
     // mexPrintf("MexFunction called with %s\n", command);

     BDDSys* given_sys = NULL;

     if (strcmp(command, "initialize") == 0)
     {
          // read arguments:
          // number of state vars, number of states, state encoding,
          // number of action vars, number of actions, action encodings
          uint state_var_num = (uint)mxGetScalar(prhs[1]);
          uint state_num = (uint)mxGetScalar(prhs[2]);
          mexPrintf("Initializing with s_vars=%d, s_num=%d\n", state_var_num, state_num);
          const mxArray* s_enc_array = prhs[3];

          array_init(s_encodings, EncList, state_num, s_buffer);
          // read contents of encoding array (cell array)
          for (int i = 0; i < state_num; i++)
          {
               array_init(new_enc, NumList, state_var_num, var_buffer);
               mxArray* cell = mxGetCell(s_enc_array, i);
               int size = mxGetNumberOfElements(cell);
               double* list_ptr = mxGetPr(cell);
               for (int k = 0; k < size; k++)
                    array_set(new_enc, k, (uint)list_ptr[k]);
               int len_diff = state_var_num - size;
               // pad with zeros
               if (len_diff > 0)
                    zeropad(new_enc, size, len_diff);
               mexPrintf("Read ");
               for (int k = 0; k < array_len(new_enc); k++)
                    mexPrintf("%d", array_get(new_enc, k));
               mexPrintf("\n");
               array_set(s_encodings, i, new_enc);
          }
          uint action_var_num = (uint)mxGetScalar(prhs[4]);
          uint action_num = (uint)mxGetScalar(prhs[5]);
          mexPrintf("Read a_vars=%d, a_num=%d\n", action_var_num, action_num);
          const mxArray* a_enc_array = prhs[6];
          array_init(a_encodings, EncList, action_num, s_buffer);
          // read contents of encoding array (cell array)
          for (int i = 0; i < action_num; i++)
          {
               mxArray* cell = mxGetCell(a_enc_array, i);
               int size = mxGetNumberOfElements(cell);
               array_init(list, NumList, action_var_num, var_buffer);
               double* list_ptr = (double*)mxGetData(cell);
               for (int k = 0; k < array_len(list); k++)
                    array_set(list, k, (uint)list_ptr[k]);

               int len_diff = action_var_num - size;
               // pad with zeros
               if (len_diff > 0)
                    zeropad(list, size, len_diff);
               mexPrintf("Read ");
               for (int k = 0; k < array_len(list); k++)
                    mexPrintf("%d", array_get(list, k));
               mexPrintf("\n");
               array_set(a_encodings, i, list);
          }

          BDDSys* sys = initializeBDD(state_var_num, s_encodings, action_var_num, a_encodings);
          // add new system to list of allocated BDD systems
          mexPrintf("%d\n", array_len(allocated_BDDs));
          array_pushback(allocated_BDDs, sys);
          mexPrintf("%d\n", array_len(allocated_BDDs));

          // return BDD system index
          plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
          double* output_ptr = mxGetPr(plhs[0]);
          *output_ptr = array_len(allocated_BDDs);
          array_free(a_encodings);
          array_free(s_encodings);
     }
     else
     {
          // if not given a "initialize" command, then system properties will
          // be changed or destroyed.
          // Read given system
          uint system_ID = (uint)mxGetScalar(prhs[1]) - 1;
          if (system_ID >= array_len(allocated_BDDs) || array_get(allocated_BDDs, system_ID) == NULL)
          {
               //TODO: Prompt user of desposed or nonexistent system
               return;
          }

          given_sys = array_get(allocated_BDDs, system_ID);
     }

     if (strcmp(command, "add_a") == 0)
     {
          /*
          read arguments:
          action ind, new number of action vars, encoding
          TODO: Error handling, size of encoding can cause segfault, read new_var_num from enc?
          */
          // mexPrintf("Reading input\n");
          uint action = (uint)mxGetScalar(prhs[2]) - 1;
          uint new_var_num = (uint)mxGetScalar(prhs[3]);
          array_init(encoding, NumList, new_var_num, var_buffer);
          uint enc_len = MAX(mxGetN(prhs[4]), mxGetM(prhs[4]));
          array_cpy(encoding, 0, (uint*)mxGetData(prhs[4]), sizeof(uint)*enc_len);
          uint len_diff = new_var_num - enc_len;
          // mexPrintf("Padding\n");
          if (len_diff > 0)
               zeropad(encoding, enc_len, len_diff);
          // mexPrintf("Adding action to system\n");
          add_action(given_sys, action, encoding, new_var_num);
     }
     else if (strcmp(command, "add_s") == 0)
     {
          /*
          read arguments:
          state ind, new number of state variables, encoding
          TODO: Error handling, size of encoding can cause segfault, read new_var_num from enc?
          */
          // mexPrintf("Reading input\n");
          uint state = (uint)mxGetScalar(prhs[2]) - 1;
          uint new_var_num = (uint)mxGetScalar(prhs[3]);
          uint enc_len = mxGetNumberOfElements(prhs[4]);
          array_init(encoding, NumList, enc_len, var_buffer);
          uint* ptr = (uint*)mxGetData(prhs[4]);
          // mexPrintf("Read encoding: ");
          for (int i = 0; i < enc_len; i++)
          {
               array_set(encoding, i, (uint)ptr[i]);
               // mexPrintf("%d", array_get(encoding, i));
          }
          // mexPrintf("\n");
          uint len_diff = new_var_num - enc_len;
          // mexPrintf("Padding by %d\n", len_diff);
          if (len_diff > 0)
               zeropad(encoding, enc_len, len_diff);
          // mexPrintf("Adding state\n");
          // mexEvalString("drawnow;");
          add_state(given_sys, state, encoding, new_var_num);
          // mexPrintf("Added state\n");
          // mexPrintf("Encoding at %d: ", state-1);
          // NumList* test_enc = array_get(given_sys->s_encs, state-1);
          // for (int j = 0; j < array_len(test_enc); j++)
          //      mexPrintf("%d", array_get(test_enc, j));
          // mexPrintf("\n");
          // mexPrintf("And at %d: ", state);
          // test_enc = array_get(given_sys->s_encs, state);
          // for (int j = 0; j < array_len(test_enc); j++)
          //      mexPrintf("%d", array_get(test_enc, j));
          // mexPrintf("\n");

     }
     else if(strcmp(command, "set_state_enc") == 0)
     {
          /*
          Read arguments:
          old state index, new state index
          */
          uint ind = (uint)mxGetScalar(prhs[2]) - 1;
          uint enc_len = mxGetNumberOfElements(prhs[3]);
          array_init(encoding, NumList, enc_len, var_buffer);
          uint* ptr = (uint*)mxGetData(prhs[3]);
          for (int i = 0; i < enc_len; i++)
               array_set(encoding, i, ptr[i]);
          // array_cpy(encoding, 0, (uint*)mxGetData(prhs[3]), sizeof(uint)*enc_len);
          uint len_diff = given_sys->s_var_num - enc_len;
          if (len_diff > 0)
               zeropad(encoding, enc_len, len_diff);
          // mexPrintf("Setting enc of length %d ", enc_len);
          // for (int i = 0; i < array_len(encoding); i++)
          //      mexPrintf("%d", array_get(encoding, i));
          // mexPrintf(" to %d\n", ind+1);
          set_state_enc(given_sys, ind, encoding);
          // for (int i = 0; i < array_len(given_sys->enc_state_map); i++)
          // {
          //      mexPrintf("key %d to %d\n", i, array_get(given_sys->enc_state_map, i));
          // }
     }
     else if (strcmp(command, "add_trans") == 0)
     {
          /*
          Read arguments:
          in state(s), action(s), out state(s)
          where states, and action lists have same dimensions
          TODO: Error handling, unequal number of input arguments,
                states or action not in system references
          */
          uint trans_num = mxGetNumberOfElements(prhs[2]);
          double* in_states_ptr = mxGetPr(prhs[2]);
          array_init(in_states, NumList, trans_num, 1);
          double* actions_ptr = mxGetPr(prhs[3]);
          array_init(actions, NumList, trans_num, 1);
          double* out_states_ptr = mxGetPr(prhs[4]);
          array_init(out_states, NumList, trans_num, 1);
          for (int i = 0; i < trans_num; i++)
          {
               array_set(in_states, i, (uint)in_states_ptr[i] - 1);
               array_set(actions, i, (uint)actions_ptr[i] - 1);
               array_set(out_states, i, (uint)out_states_ptr[i] - 1);
               // mexPrintf("C Adding trans: (%d, %d, %d)\n",
               //           array_get(in_states, i),
               //           array_get(actions, i),
               //           array_get(out_states, i));
               // mexEvalString("drawnow;");
          }

          // mexPrintf("Adding %d transitions\n", trans_num);
               // mexPrintf("Adding a transition (%d, %d, %d)\n", (uint)in_states[i], (uint)actions[i], (uint)out_states[i]);
          add_trans(given_sys, in_states, actions, out_states);
     }
     else if (strcmp(command, "get_trans_with_s") == 0)
     {
          /*
          Read arguments:
          state
          */
          // mexPrintf("Gets in\n");
          double* input = mxGetPr(prhs[2]);
          uint state_num = mxGetNumberOfElements(prhs[2]);
          array_init(states, NumList, state_num, s_buffer);
          for (int i = 0; i < state_num; i++)
          {
               // mexPrintf("Read %d\n", (uint)input[i]);
               array_set(states, i, (uint)input[i] - 1);
          }
          // mexPrintf("Getting transitions\n");
          DdNode* transitions = get_trans_with_s(given_sys, states);
          array_free(states);

          // mexPrintf("Reading transitions");
          // mexEvalString("drawnow;");
          uint** transition_parts = read_transitions(given_sys, transitions);
          Cudd_RecursiveDeref(given_sys->manager, transitions);
          uint trans_num = (*transition_parts)[0];


          if (nlhs != 3)
          {
               // TODO: Output argument errors
          }
          // mexPrintf("Giving output\n");
          // create outputs
          for (int i = 0; i < 3; i++)
          {
               plhs[i] = mxCreateDoubleMatrix(trans_num, 1, mxREAL);
               double* ptr = (double*)mxGetData(plhs[i]);
               for (int j = 0; j < trans_num; j++)
               {
                    // mexPrintf("output %d: %d\n", i, transition_parts[i+1][j]);
                    ptr[j] = (double)(transition_parts[i+1][j]+1);
               }
               mxFree(transition_parts[i+1]);
          }
          mxFree(transition_parts);
     }
     else if (strcmp(command, "rm_trans_with_s") == 0)
     {
          /*
          Read arguments:
          state
          */
          uint state_num = mxGetNumberOfElements(prhs[2]);
          double* input = mxGetPr(prhs[2]);
          array_init(states, NumList, state_num, 1);
          for (int i = 0; i < state_num; i++)
               array_set(states, i, (uint)input[i]-1);
          rm_trans_with_s(given_sys, states);
     }
     else if (strcmp(command, "add_pg") == 0)
     {
          /*
          Read arguments:
          action group, state group
          */
          // mexPrintf("Reading inputs\n");
          uint U_len = mxGetNumberOfElements(prhs[2]);
          uint G_len = mxGetNumberOfElements(prhs[3]);
          array_init(U, NumList, U_len, 1);
          array_init(G, NumList, G_len, 1);
          uint* output_U = (uint*)mxGetPr(prhs[2]);
          double* output_G = mxGetPr(prhs[3]);
          // mexPrintf("Setting groups\n");
          for (int i = 0; i < U_len; i++)
               array_set(U, i, (uint)output_U[i] - 1);
          for (int i = 0; i < G_len; i++)
               array_set(G, i, (uint)output_G[i] - 1);
          // mexPrintf("Adding groups\n");
          add_progress_group(given_sys, U, G);
          array_free(U);
          array_free(G);
     }
     else if (strcmp(command, "read_pg") == 0)
     {
          uint ind = (uint)mxGetScalar(prhs[2]) - 1;
          // mexPrintf("Checking group %d\n", ind);
          if (array_len(given_sys->pg_U) < ind
               || array_len(given_sys->pg_G) < ind
               || array_get(given_sys->pg_U, ind) == NULL
               || array_get(given_sys->pg_G, ind) == NULL)
          {
               plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
               plhs[1] = mxCreateDoubleMatrix(0, 0, mxREAL);
               return;
          }

          uint U_len = array_len(given_sys->pg_U);
          // mexPrintf("U_group_len: %d\n", U_len);
          uint G_len = array_len(given_sys->pg_G);
          // mexPrintf("G_group_len: %d\n", G_len);

          // mexPrintf("Reading groups %d\n", ind);
          uint** transition_parts = read_actions(given_sys, array_get(given_sys->pg_U, ind));
          uint U_num = *(transition_parts[0]);
          // mexPrintf("writing read parts\n");
          plhs[0] = mxCreateNumericMatrix(U_num, 1, mxUINT32_CLASS, mxREAL);
          uint* ptr = (uint*)mxGetData(plhs[0]);
          // mexPrintf("U_len: %d\n", U_num);
          for (int i = 0; i < U_num; i++)
          {
               // mexPrintf("Read U=%d\n", transition_parts[1][i] + 1);
               ptr[i] = (uint)(transition_parts[1][i]+1);
          }
          // mexPrintf("Freeing middle parts\n");
          mxFree(transition_parts[1]);
          mxFree(transition_parts);

          transition_parts = read_out_states(given_sys, array_get(given_sys->pg_G, ind));
          uint G_num = *(transition_parts[0]);
          plhs[1] = mxCreateNumericMatrix(G_num, 1, mxUINT32_CLASS, mxREAL);
          uint* G_ptr = (uint*)mxGetData(plhs[1]);
          for (int i = 0; i < G_num; i++)
          {
               // mexPrintf("Read G=%d\n", transition_parts[1][i]);
               G_ptr[i] = transition_parts[1][i]+(uint)1;
          }
          mxFree(transition_parts[1]);
          mxFree(transition_parts);
     }
     else if(strcmp(command, "check_sup_pg") == 0)
     {
          /*
          Read arguments
          action group, state group
          */

          uint U_len = mxGetNumberOfElements(prhs[2]);
          uint G_len = mxGetNumberOfElements(prhs[3]);
          // mexPrintf("Reading groups of length: U=%d,G=%d\n", U_len, G_len);
          uint* input_U = (uint*)mxGetData(prhs[2]);
          double* input_G = mxGetPr(prhs[3]);
          array_init(U, NumList, U_len, 1);
          array_init(G, NumList, G_len, 1);
          // mexPrintf("U group: ");
          for (int i = 0; i < U_len; i++)
          {
               // mexPrintf("%d ", (uint)input_U[i]);
               array_set(U, i, (uint)input_U[i] - 1);
          }
          // mexPrintf("\n");
          // mexPrintf("G group: ");
          for (int i = 0; i < G_len; i++)
          {
               // mexPrintf("%d ", (uint)input_G[i]);
               array_set(G, i, (uint)input_G[i] - 1);
          }
          // mexPrintf("\n");
          int has_sup = has_superior_pg(given_sys, U, G);
          array_free(U);
          array_free(G);
          // mexPrintf("has_sup=%d\n", has_sup);
          // output: true if U and G has superior group
          plhs[0] = mxCreateLogicalScalar(has_sup);
     }
     else if (strcmp(command, "rm_pg") == 0)
     {
          /*
          Read arguments:
          index of progress group to remove
          */
          uint ind = (uint)mxGetScalar(prhs[2]) - 1;
          rm_progress_group(given_sys, ind);
     }
     else if (strcmp(command, "add_to_pg") == 0)
     {
          /*
          Read arguments:
          indices of progress groups, state elements to add
          */
          double* input_inds = mxGetPr(prhs[2]);
          uint ind_len = mxGetNumberOfElements(prhs[2]);
          if (ind_len == 0)
          {
               command = "empty";
               goto end_of_mex;
          }
          array_init(inds, NumList, ind_len, 1);
          for (int i = 0; i < ind_len; i++)
          {
               // mexPrintf("Read index %d\n", (uint)input_inds[i]);
               array_set(inds, i, (uint)input_inds[i]-1);
          }
          double* input_states = mxGetPr(prhs[3]);
          uint state_num = mxGetNumberOfElements(prhs[3]);
          if (state_num == 0)
          {
               command = "empty";
               array_free(inds);
               goto end_of_mex;
          }
          array_init(states, NumList, state_num, 1);
          for (int i = 0; i < state_num; i++)
          {
               // mexPrintf("Read state %d\n", (uint)input_states[i]);
               array_set(states, i, (uint)input_states[i] - 1);
          }
          add_to_pg(given_sys, inds, states);

          array_free(inds);
          array_free(states);
     }
     else if (strcmp(command, "check_mem_in_G_pg") == 0)
     {
          /*
          Read arguments:
          members to check if inferior pg group
          */
          uint* input_states = (uint*)mxGetData(prhs[2]);
          uint state_num = mxGetNumberOfElements(prhs[2]);
          array_init(states, NumList, state_num, 1);
          for (int i = 0; i < state_num; i++)
          {
               // mexPrintf("%d\n", input_states[i]-1);
               array_set(states, i, input_states[i] - 1);
          }
          // mexPrintf("G group: ");
          // for (int i = 0; i < state_num; i++)
          //      mexPrintf("%d ", array_get(states, i));
          // mexPrintf("\n");
          // returns list (length, elements)
          NumList* group_indices = is_member_of_pg(given_sys, states);
          //output: indices of groups it is part of
          array_free(states);
          uint len = array_len(group_indices);

          // mexPrintf("Making output\n");
          if (array_len(group_indices) > 0)
          {
               // mexPrintf("Writing to output\n");
               // mexEvalString("drawnow;");
               plhs[0] = mxCreateDoubleMatrix(1, array_len(group_indices), mxREAL);
               double* storage = mxGetPr(plhs[0]);

               for (int i = 0; i < array_len(group_indices); i++)
                    storage[i] = (double)array_get(group_indices, i) + 1;
          }
          else
               plhs[0] = mxCreateDoubleMatrix(0,0, mxREAL);
          array_free(group_indices);
     }
     else if (strcmp(command, "delete") == 0)
     {
          /*
               Read arguments:
               ID of system to free
          */
          uint ID = (uint)mxGetScalar(prhs[1]) - 1;
          if (ID < array_len(allocated_BDDs))
          {
               free_system(array_get(allocated_BDDs, ID));
               array_set(allocated_BDDs, ID, NULL);
          }
          else
          {
               //TODO: ERROR BAD INDEX
          }

     }
     else if (strcmp(command, "pre") == 0)
     {
          /*
          Read arguments:
          post-states, actions, quant1, quant2
          */
          // reading states
          DdNode* state_BDD = read_BDD_from_mxarray(given_sys, prhs[2], 's');

          // reading actions
          uint action_num = mxGetNumberOfElements(prhs[3]);
          DdNode* action_BDD = read_BDD_from_mxarray(given_sys, prhs[3], 'a');
          if (action_num == 0)
          {
               Cudd_RecursiveDeref(given_sys->manager, action_BDD);
               action_BDD = Cudd_ReadOne(given_sys->manager);
               Cudd_Ref(action_BDD);
          }

          // reading quantifiers
          char quant1 = ((uint)(mxGetScalar(prhs[4]))) ? 'e' : 'a';
          char quant2 = ((uint)(mxGetScalar(prhs[5]))) ? 'e' : 'a';

          mexPrintf("Calling with quants (%c, %c)\n", quant1, quant2);
          // calculate set
          DdNode* pre_set = pre(given_sys, state_BDD, action_BDD, quant1, quant2);
          uint** pre_states = read_in_states(given_sys, pre_set);
          Cudd_RecursiveDeref(given_sys->manager, pre_set);
          Cudd_RecursiveDeref(given_sys->manager, state_BDD);
          Cudd_RecursiveDeref(given_sys->manager, action_BDD);
          // return output
          uint pre_num = (uint)*(pre_states[0]);
          plhs[0] = mxCreateDoubleMatrix(pre_num, 1, mxREAL);
          double* ptr = mxGetPr(plhs[0]);
          mexPrintf("We have %d elements\n", pre_num);
          for (int i = 0; i < pre_num; i++)
               ptr[i] = (double)pre_states[1][i] + 1;
     }
     else if (strcmp(command, "pre_pg") == 0)
     {
          // read V
          DdNode* V = read_BDD_from_mxarray(given_sys, prhs[2], 's');

          // read B
          DdNode* B = read_BDD_from_mxarray(given_sys, prhs[3], 's');

          // read quant
          char quant = (uint)mxGetScalar(prhs[4]) ? 'e' : 'a';
          int mode = nlhs;

          DdNode** PG_sets = PGpre(given_sys, V, B, quant, mode);
          Cudd_RecursiveDeref(given_sys->manager, V);
          Cudd_RecursiveDeref(given_sys->manager, B);

          if (mode >= WIN_SET)
          {
               DdNode* PG_set_W = PG_sets[0];
               write_BDD_to_mxarray(given_sys, PG_set_W, &plhs[0], 's');
               Cudd_RecursiveDeref(given_sys->manager, PG_set_W);
               if (mode >= WIN_CANDIDATE_SET)
               {
                    DdNode* PG_set_Cw = PG_sets[1];
                    write_BDD_to_mxarray(given_sys, PG_set_Cw, &plhs[1], 's');
                    Cudd_RecursiveDeref(given_sys->manager, PG_set_Cw);
               }
          }
          mxFree(PG_sets);
     }
     else if (strcmp(command, "pg_inv") == 0)
     {
          // read U
          DdNode* U = read_BDD_from_mxarray(given_sys, prhs[2], 'a');

          // read G
          DdNode* G = read_BDD_from_mxarray(given_sys, prhs[3], 's');

          // read Z
          DdNode* Z = read_BDD_from_mxarray(given_sys, prhs[4], 's');

          // read B
          DdNode* B = read_BDD_from_mxarray(given_sys, prhs[5], 's');

          // read quant1
          char quant1 = (uint)mxGetScalar(prhs[6]) ? 'e' : 'a';

          // set mode
          int mode = nlhs;
          mexPrintf("Calling pg_inv with quant %c and mode %d\n", quant1, mode);
          DdNode** pg_inv_set = inv(given_sys, Z, B, U, G, quant1, mode);
          Cudd_RecursiveDeref(given_sys->manager, U);
          Cudd_RecursiveDeref(given_sys->manager, G);
          Cudd_RecursiveDeref(given_sys->manager, Z);
          Cudd_RecursiveDeref(given_sys->manager, B);

          if (mode >= WIN_SET)
          {
               DdNode* pg_inv_set_W = pg_inv_set[0];
               write_BDD_to_mxarray(given_sys, pg_inv_set_W, &plhs[0], 's');
               Cudd_RecursiveDeref(given_sys->manager, pg_inv_set_W);
               if (mode >= WIN_CANDIDATE_SET)
               {
                    DdNode* pg_inv_set_Cw = pg_inv_set[0];
                    write_BDD_to_mxarray(given_sys, pg_inv_set_Cw, &plhs[1], 's');
                    Cudd_RecursiveDeref(given_sys->manager, pg_inv_set_Cw);
               }
          }
          mxFree(pg_inv_set);
     }
     else if (strcmp(command, "win_until") == 0)
     {
          // read B
          DdNode* B = read_BDD_from_mxarray(given_sys, prhs[2], 's');

          // read P
          DdNode* P = read_BDD_from_mxarray(given_sys, prhs[3], 's');

          // read quant1
          char quant = (uint)mxGetScalar(prhs[4]) ? 'e' : 'a';
          int mode = nlhs;

          DdNode** win_sets = win_until(given_sys, P, B, quant, mode);
          Cudd_RecursiveDeref(given_sys->manager, B);
          Cudd_RecursiveDeref(given_sys->manager, P);

          if (mode >= WIN_SET)
          {
               write_BDD_to_mxarray(given_sys, win_sets[0], &plhs[0], 's');
               Cudd_RecursiveDeref(given_sys->manager, win_sets[0]);
               if (mode >= WIN_CANDIDATE_SET)
               {
                    write_BDD_to_mxarray(given_sys, win_sets[1], &plhs[1], 's');
                    Cudd_RecursiveDeref(given_sys->manager, win_sets[1]);
               }
          }
          mxFree(win_sets);
     }
     else if (strcmp(command, "win_until_and_always") == 0)
     {
          // read A
          DdNode* A = read_BDD_from_mxarray(given_sys, prhs[2], 's');

          // read B
          DdNode* B = read_BDD_from_mxarray(given_sys, prhs[3], 's');

          // read P
          DdNode* P = read_BDD_from_mxarray(given_sys, prhs[4], 's');

          // quant1
          char quant = (uint)mxGetScalar(prhs[5]) ? 'e' : 'a';
          int mode = nlhs;

          DdNode** win_sets = win_until_and_always(given_sys, A, B, P, quant, mode);
          Cudd_RecursiveDeref(given_sys->manager, A);
          Cudd_RecursiveDeref(given_sys->manager, B);
          Cudd_RecursiveDeref(given_sys->manager, P);

          if (win_sets == NULL)
               return;
          if (mode >= WIN_SET)
          {
               write_BDD_to_mxarray(given_sys, win_sets[0], &plhs[0], 's');
               Cudd_RecursiveDeref(given_sys->manager, win_sets[0]);
               if (mode >= WIN_CANDIDATE_SET)
               {
                    write_BDD_to_mxarray(given_sys, win_sets[1], &plhs[1], 's');
                    Cudd_RecursiveDeref(given_sys->manager, win_sets[1]);
               }
          }
          mxFree(win_sets);
     }
     else if (strcmp(command, "win_intermediate") == 0)
     {
          // A
          DdNode* A = read_BDD_from_mxarray(given_sys, prhs[2], 's');

          // B
          DdNode* B = read_BDD_from_mxarray(given_sys, prhs[3], 's');

          // P
          DdNode* P = read_BDD_from_mxarray(given_sys, prhs[4], 's');

          // C_list
          uint cell_num = mxGetNumberOfElements(prhs[5]);
          DdNode** C_list = mxMalloc(sizeof(DdNode*)*cell_num);
          for (int i = 0; i < cell_num; i++)
          {
               const mxArray* cell = mxGetCell(prhs[5], i);
               C_list[i] = read_BDD_from_mxarray(given_sys, cell, 's');
          }

          // quant
          char quant = (uint)mxGetScalar(prhs[6]) ? 'e' : 'a';
          int mode = nlhs;

          DdNode** win_sets = win_intermediate(given_sys, A, B, P, C_list, cell_num, quant, mode);
          Cudd_RecursiveDeref(given_sys->manager, A);
          Cudd_RecursiveDeref(given_sys->manager, B);
          Cudd_RecursiveDeref(given_sys->manager, P);
          for (int i = 0; i < cell_num; i++)
               Cudd_RecursiveDeref(given_sys->manager, C_list[i]);

          if (mode >= WIN_SET)
          {
               write_BDD_to_mxarray(given_sys, win_sets[0], &plhs[0], 's');
               Cudd_RecursiveDeref(given_sys->manager, win_sets[0]);
               if (mode >= WIN_CANDIDATE_SET)
               {
                    write_BDD_to_mxarray(given_sys, win_sets[1], &plhs[1], 's');
                    Cudd_RecursiveDeref(given_sys->manager, win_sets[1]);
               }
          }
          mxFree(win_sets);
          mxFree(C_list);
     }
     else if (strcmp(command, "win_primal") == 0)
     {
          // A
          DdNode* A = read_BDD_from_mxarray(given_sys, prhs[2], 's');
          // B
          DdNode* B = read_BDD_from_mxarray(given_sys, prhs[3], 's');
          // C_list
          uint cell_num = mxGetNumberOfElements(prhs[4]);
          DdNode** C_list = mxMalloc(sizeof(DdNode*)*cell_num);
          for (int i = 0; i < cell_num; i++)
          {
               const mxArray* cell = mxGetCell(prhs[4], i);
               C_list[i] = read_BDD_from_mxarray(given_sys, cell, 's');
          }

          // quants
          char quant1 = (uint)mxGetScalar(prhs[5]) ? 'e' : 'a';
          char quant2 = (uint)mxGetScalar(prhs[6]) ? 'e' : 'a';

          // head_start
          DdNode* V = read_BDD_from_mxarray(given_sys, prhs[7], 's');

          // mode
          int mode = nlhs;

          DdNode** win_sets = win_primal(given_sys, A, B, C_list, cell_num, quant1, quant2, V, mode);
          Cudd_RecursiveDeref(given_sys->manager, A);
          Cudd_RecursiveDeref(given_sys->manager, B);
          Cudd_RecursiveDeref(given_sys->manager, V);
          for (int i = 0; i < cell_num; i++)
               Cudd_RecursiveDeref(given_sys->manager, C_list[i]);

          if (mode >= WIN_SET)
          {
               write_BDD_to_mxarray(given_sys, win_sets[0], &plhs[0], 's');
               Cudd_RecursiveDeref(given_sys->manager, win_sets[0]);
               if (mode >= WIN_CANDIDATE_SET)
               {
                    write_BDD_to_mxarray(given_sys, win_sets[1], &plhs[1], 's');
                    Cudd_RecursiveDeref(given_sys->manager, win_sets[1]);
               }
          }
          mxFree(win_sets);
          mxFree(C_list);
     }
     else if (strcmp(command, "print_enc_map") == 0)
     {
          mexPrintf("Enc state map:\n");
          for (int i = 0; i < array_len(given_sys->enc_state_map); i++)
               mexPrintf("%d to %d\n", i, array_get(given_sys->enc_state_map, i));
     }
     else if (strcmp(command, "print_all_states") == 0)
     {
          mexPrintf("All states:\n");
          uint** parts = read_out_states(given_sys, given_sys->all_states);
          uint len = *(parts[0]);
          for (int i = 0; i < len; i++)
               mexPrintf("%d\n", parts[1][i]);
          mxFree(parts[1]);
          mxFree(parts);
     }
     else if (strcmp(command, "print_trans") == 0)
     {
          mexPrintf("All transitions:\n");
          uint** parts = read_transitions(given_sys, given_sys->trans_sys);
          uint len = *(parts[0]);
          for (int i = 0; i < len; i++)
          {
               mexPrintf("(%d, %d, %d)\n", parts[1][i], parts[2][i], parts[3][i]);
          }
          mxFree(parts[1]);
          mxFree(parts[2]);
          mxFree(parts[3]);
          mxFree(parts);
     }
     else if (strcmp(command, "num_trans") == 0)
     {
          plhs[0] = mxCreateNumericMatrix(1,1, mxUINT32_CLASS, mxREAL);
          uint* ptr = (uint*)mxGetData(plhs[0]);
          uint trans_num = get_num_of_trans(given_sys, given_sys->trans_sys);
          *ptr = trans_num;
     }
     else if (strcmp(command, "print_encs") == 0)
     {
          mexPrintf("Encodings:\n");
          for (int i = 0; i < array_len(given_sys->s_encs); i++)
          {
               mexPrintf("%d ", i);
               NumList* enc = array_get(given_sys->s_encs, i);
               for (int j = 0; j < array_len(enc); j++)
               {
                    mexPrintf("%d", array_get(enc, j));
               }
               mexPrintf("\n");
          }
     }
     else if (strcmp(command, "debug") == 0)
     {
          if (Cudd_DebugCheck(given_sys->manager))
          {
               mexPrintf("Something is wrong!\n");
          }
     }
     else if (strcmp(command, "count_nodes") == 0)
     {
          plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
          double* ptr = mxGetPr(plhs[0]);
          *ptr = (double)Cudd_ReadNodeCount(given_sys->manager);
     }
     else if (strcmp(command, "count_system_nodes") == 0)
     {
          plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
          double* ptr = mxGetPr(plhs[0]);
          *ptr = (double)Cudd_DagSize(given_sys->trans_sys);
     }
     else if (strcmp(command, "count_dead_nodes") == 0)
     {
          plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
          double* ptr = mxGetPr(plhs[0]);
          *ptr = (double)Cudd_ReadDead(given_sys->manager);
     }
     else if (strcmp(command, "reorder") == 0)
     {
          uint bound = (uint)mxGetScalar(prhs[2]);
          mexPrintf("Reorder with bound %d\n", bound);
          // mexEvalString("drawnow;");
          Cudd_ReduceHeap(given_sys->manager, reorder_alg, bound);
     }
     else if (strcmp(command, "toggle_reorder") == 0)
     {
          int is_on = *mxGetLogicals(prhs[2]);
          if (is_on)
               Cudd_AutodynEnable(given_sys->manager, reorder_alg);
          else
               Cudd_AutodynDisable(given_sys->manager);
     }
     else
     {
          // unknown command
          // TODO: Error
     }

     end_of_mex:
     if (strcmp(command, "empty") == 0)
     {
          for (int i = 0; i < nlhs; i++)
          {
               plhs[i] = mxCreateDoubleMatrix(0,0,mxREAL);
          }
     }
}

/*
     Initialize the BDD and needed fields
     Make order of variables be (a_vars, in_states, out_states)
*/
BDDSys* initializeBDD(uint var_num, EncList* _s_encs, uint a_var_num, EncList* _a_encs)
{
     BDDSys sys;
     mexPrintf("Making system with s=%d, a=%d\n", var_num, a_var_num);
     sys.s_var_num = var_num;
     sys.a_var_num = a_var_num;
     /* Start the BDD manager
     2*var_num + a_var_num initial BDD variables
     0 ZDD variables
     default number of unique slots (cudd figures it out)
     default cache size (cudd figures it out)
     max memory 0 = unlimited */

     mexPrintf("Making manager with %u\n", 2*var_num + a_var_num);
     sys.manager = Cudd_Init((int)(2*var_num + a_var_num),0,CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS,0);
     // Cudd_AutodynEnable(sys.manager, reorder_alg);
     // Gather variables and allocate variable lists
     mexPrintf("Creating action variables\n");
     array_init(a_vars, BDDlist, sys.a_var_num, var_buffer);
     array_init(a_inds, NumList, sys.a_var_num, var_buffer);
     sys.a_vars = a_vars;
     sys.a_inds = a_inds;
     array_make_persist(sys.a_vars);
     array_make_persist(sys.a_inds);
     for (int i = 0; i < sys.a_var_num; i++)
     {
          array_set(sys.a_vars, i, Cudd_bddIthVar(sys.manager, i));
          array_set(sys.a_inds, i, i);
     }
     mexPrintf("Creating state variables\n");
     array_init(s_in_vars, BDDlist, sys.s_var_num, var_buffer);
     array_init(s_in_inds, NumList, sys.s_var_num, var_buffer);
     sys.s_in_vars = s_in_vars;
     sys.s_in_inds = s_in_inds;
     array_make_persist(sys.s_in_vars);
     array_make_persist(sys.s_in_inds);
     for (int i = 0; i < sys.s_var_num; i++)
     {
          array_set(sys.s_in_vars, i, Cudd_bddIthVar(sys.manager, i + sys.a_var_num));
          array_set(sys.s_in_inds, i, i+sys.a_var_num);
     }
     array_init(s_out_vars, BDDlist, sys.s_var_num, var_buffer);
     array_init(s_out_inds, NumList, sys.s_var_num, var_buffer);
     sys.s_out_vars = s_out_vars;
     sys.s_out_inds = s_out_inds;
     array_make_persist(sys.s_out_vars);
     array_make_persist(sys.s_out_inds);
     for (int i = 0; i < sys.s_var_num; i++)
     {
          array_set(sys.s_out_vars, i, Cudd_bddIthVar(sys.manager, i + sys.a_var_num + sys.s_var_num));
          array_set(sys.s_out_inds, i, i + sys.a_var_num + sys.s_var_num);
     }
     // save and allocate encodings and add states to sets
     mexPrintf("Allocating encoding lists\n");
     array_init(s_encs, EncList, array_len(_s_encs), s_buffer);
     array_init(enc_map, NumList, pow(2, sys.s_var_num), pow(2, sys.s_var_num));

     sys.enc_state_map = enc_map;
     sys.all_states = Cudd_ReadLogicZero(sys.manager);
     Cudd_Ref(sys.all_states);
     array_make_persist(sys.enc_state_map);

     mexPrintf("Adding %d encodings\n", array_len(_s_encs));
     for (int i = 0; i < array_len(_s_encs); i++)
     {
          NumList* enc = array_get(_s_encs, i);
          array_set(s_encs, i, enc);
          array_make_persist(array_get(s_encs, i));
          uint key = fromBitArray(array_list(enc), sys.s_var_num);
          array_set(sys.enc_state_map, key, (uint)i);

          mexPrintf("Adding enc key : %d\n", key);
          // Add to set
          DdNode* state_cube = Cudd_bddComputeCube(sys.manager, array_list(sys.s_out_vars), array_list(enc), sys.s_var_num);
          Cudd_Ref(state_cube);
          DdNode* tmp = Cudd_bddOr(sys.manager, sys.all_states, state_cube);
          Cudd_Ref(tmp);
          Cudd_RecursiveDeref(sys.manager, sys.all_states);
          Cudd_RecursiveDeref(sys.manager, state_cube);
          sys.all_states = tmp;
     }
     sys.s_encs = s_encs;
     array_make_persist(sys.s_encs);

     array_init(a_encs, EncList, array_len(_a_encs), s_buffer);
     sys.all_actions = Cudd_ReadLogicZero(sys.manager);
     Cudd_Ref(sys.all_actions);
     for (int i = 0; i < array_len(_a_encs); i++)
     {
          NumList* enc = array_get(_a_encs, i);
          array_set(a_encs, i, enc);
          array_make_persist(array_get(a_encs, i));
          printf("Adding action encoding %d: ", i);
          for (int i = 0; i < array_len(enc); i++)
               printf("%d", array_get(enc, i));
          printf("\n");

          // Add to set
          DdNode* action_cube = Cudd_bddComputeCube(sys.manager, array_list(sys.a_vars), array_list(enc), sys.a_var_num);
          Cudd_Ref(action_cube);
          DdNode* tmp = Cudd_bddOr(sys.manager, sys.all_actions, action_cube);
          Cudd_Ref(tmp);
          Cudd_RecursiveDeref(sys.manager, action_cube);
          Cudd_RecursiveDeref(sys.manager, sys.all_actions);
          sys.all_actions = tmp;
     }
     sys.a_encs = a_encs;
     array_make_persist(sys.a_encs);
     // allocate pg groups
     array_init(pg_U, BDDlist, 0, pg_buffer);
     array_init(pg_G, BDDlist, 0, pg_buffer);
     sys.pg_U = pg_U;
     array_make_persist(sys.pg_U);
     sys.pg_G = pg_G;
     array_make_persist(sys.pg_G);

     sys.trans_sys = Cudd_ReadLogicZero(sys.manager);
     Cudd_Ref(sys.trans_sys);

     BDDSys* ptr = mxMalloc(sizeof(BDDSys));
     *ptr = sys;
     mexMakeMemoryPersistent(ptr);
     return ptr;
}

// pads the uint list with the amount of zeros specified by spaces
void zeropad(NumList* list, uint pos, uint pad_len)
{
     uint zero = (uint)0;
     for (int i = 0; i < pad_len; i++)
     {
          if (pos+i < array_len(list))
          {
               array_set(list, pos+i, zero);
          }
          else
               array_pushback(list, zero);
     }
}

void add_action(BDDSys* sys, uint index, NumList* enc, uint new_var_num)
{
     // mexPrintf("Calculating variable difference\n");
     // mexEvalString("drawnow;");
     int var_diff = new_var_num - sys->a_var_num;
     // mexPrintf("Number of new variables: %d\n", var_diff);
     // mexEvalString("drawnow");
     // Add new variables
     if (new_var_num > sys->a_var_num)
     {
          // zero pad encodings
          // mexPrintf("Padding encodings: %d\n", array_len(sys->a_encs));
          for (int i = 0; i < array_len(sys->a_encs); i++)
          {
               // mexPrintf("Padding enc %d\n", i);
               zeropad(array_get(sys->a_encs, i), sys->a_var_num, var_diff);

          }
          // create action variables and update BDD
          for (int i = 0; i < var_diff; i++)
          {
               // create new var
               DdNode* new_var = Cudd_bddNewVarAtLevel(sys->manager, sys->a_var_num);
               array_pushback(sys->a_vars, new_var);
               array_pushback(sys->a_inds, Cudd_ReadSize(sys->manager)-1);

               // update BDD. T & !v
               DdNode* tmp = Cudd_bddAnd(sys->manager, sys->trans_sys, Cudd_Not(new_var));
               Cudd_Ref(tmp);
               Cudd_RecursiveDeref(sys->manager, sys->trans_sys);
               sys->trans_sys = tmp;

               // update all actions. A & !v
               tmp = Cudd_bddAnd(sys->manager, sys->all_actions, Cudd_Not(new_var));
               Cudd_Ref(tmp);
               Cudd_RecursiveDeref(sys->manager, sys->all_actions);
               sys->all_actions = tmp;

               //update progress groups
               for (int i = 0; i < array_len(sys->pg_U); i++)
               {
                    if (array_get(sys->pg_U, i) != NULL)
                    {
                         tmp = Cudd_bddAnd(sys->manager, array_get(sys->pg_U, i), Cudd_Not(new_var));
                         Cudd_Ref(tmp);
                         Cudd_RecursiveDeref(sys->manager, array_get(sys->pg_U, i));
                         array_set(sys->pg_U, i, tmp);
                    }
               }
          }
          sys->a_var_num = new_var_num;
     }

     // could need insertion
     // mexPrintf("Inserting\n");
     array_insert(sys->a_encs, index, enc);
     array_make_persist(enc);

     // add action to all_actions set
     // mexPrintf("Adding to set\n");
     DdNode* action_cube = Cudd_bddComputeCube(sys->manager, array_list(sys->a_vars), array_list(enc), new_var_num);
     Cudd_Ref(action_cube);
     // mexPrintf("Created cube\n");
     DdNode* tmp = Cudd_bddOr(sys->manager, sys->all_actions, action_cube);
     Cudd_Ref(tmp);
     // mexPrintf("Created union\n");
     Cudd_RecursiveDeref(sys->manager, action_cube);
     Cudd_RecursiveDeref(sys->manager, sys->all_actions);
     sys->all_actions = tmp;

     // debug info
     // mexPrintf("Result encodings: %d\n", array_len(sys->a_encs));
     // for (int i = 0; i < array_len(sys->a_encs); i++)
     // {
     //      // mexPrintf("Writing encoding of len %d: ", array_len(array_get(sys->a_encs, i)));
     //      for (int j = 0; j < array_len(array_get(sys->a_encs, i)); j++)
     //      {
     //           // mexPrintf("%d", array_get(array_get(sys->a_encs, i), j));
     //      }
     //      // mexPrintf(" ");
     // }
     // // mexPrintf("\n");
     // uint** parts = read_actions(sys, sys->all_actions);
     // uint len = *parts[0];
     // for (int i = 0; i < len; i++)
     //      mexPrintf("%d ", parts[1][i]);
     // mexPrintf("\n");
     // mxFree(parts[1]);
     // mxFree(parts);
}

void add_state(BDDSys* sys, uint index, NumList* enc, uint new_var_num)
{
     int var_diff = new_var_num - sys->s_var_num;
     // mexPrintf("New variables: %d\n", var_diff);
     // mexEvalString("drawnow;");
     // Add new variables
     if (new_var_num > sys->s_var_num)
     {
          // zero pad encodings
          // mexPrintf("Padding all encodings\n");
          // mexEvalString("drawnow;");
          for (int i = 0; i < array_len(sys->s_encs); i++)
               zeropad(array_get(sys->s_encs, i), sys->s_var_num, var_diff);

          // create state variables and update BDD
          for (int i = 0; i < var_diff; i++)
          {
               // mexPrintf("Adding variable: %d\n", i);
               // mexEvalString("drawnow;");
               // create new var
               DdNode* new_in_var = Cudd_bddNewVarAtLevel(sys->manager, sys->a_var_num + sys->s_var_num);
               array_pushback(sys->s_in_vars, new_in_var);
               array_pushback(sys->s_in_inds, Cudd_ReadSize(sys->manager)-1);
               DdNode* new_out_var = Cudd_bddNewVar(sys->manager);
               array_pushback(sys->s_out_vars, new_out_var);
               array_pushback(sys->s_out_inds, Cudd_ReadSize(sys->manager)-1);

               // update BDD. T & !v_in & !v_out = T & !(v_in | v_out)
               DdNode* var_bdd = Cudd_bddOr(sys->manager, new_in_var, new_out_var);
               Cudd_Ref(var_bdd);
               DdNode* tmp = Cudd_bddAnd(sys->manager, sys->trans_sys, Cudd_Not(var_bdd));
               Cudd_Ref(tmp);
               Cudd_RecursiveDeref(sys->manager, var_bdd);
               Cudd_RecursiveDeref(sys->manager, sys->trans_sys);
               sys->trans_sys = tmp;

               // update all states. S & !(v_out)
               tmp = Cudd_bddAnd(sys->manager, sys->all_states, Cudd_Not(new_out_var));
               Cudd_Ref(tmp);
               Cudd_RecursiveDeref(sys->manager, sys->all_states);
               sys->all_states = tmp;

               //update progress groups
               for (int i = 0; i < array_len(sys->pg_G); i++)
               {
                    if (array_get(sys->pg_G, i) != NULL)
                    {
                         tmp = Cudd_bddAnd(sys->manager, array_get(sys->pg_G, i), Cudd_Not(new_out_var));
                         Cudd_Ref(tmp);
                         Cudd_RecursiveDeref(sys->manager, array_get(sys->pg_G, i));
                         array_set(sys->pg_G, i, tmp);
                    }
               }
          }
          sys->s_var_num = new_var_num;
     }

     // mexEvalString("memory");
     // mexPrintf("Inserting state\n");
     // mexEvalString("drawnow;");
     array_insert(sys->s_encs, index, enc);
     array_make_persist(enc);
     // mexPrintf("Encoding ");
     // mexEvalString("drawnow;");
     // for (int i = 0; i < array_len(enc); i++)
     //      mexPrintf("%d", array_get(enc, i));
     // mexEvalString("drawnow;");
     // mexPrintf(" inserted at %d\n", index);
     // mexEvalString("drawnow;");
     if (array_len(sys->s_encs) > array_len(sys->enc_state_map))
     {
          array_growby(sys->enc_state_map, array_len(sys->enc_state_map));
     }
     uint key = fromBitArray(array_list(enc), array_len(enc));
     while (key >= array_len(sys->enc_state_map))
          array_pushback(sys->enc_state_map, 0);
     array_set(sys->enc_state_map, key, index);
     // mexEvalString("memory");

     // add state to all_states (out) set
     // mexPrintf("Adding state to union\n");
     // mexEvalString("drawnow;");
     DdNode* state_cube = Cudd_bddComputeCube(sys->manager, array_list(sys->s_out_vars), array_list(enc), array_len(enc));
     Cudd_Ref(state_cube);
     DdNode* tmp = Cudd_bddOr(sys->manager, sys->all_states, state_cube);
     Cudd_Ref(tmp);
     Cudd_RecursiveDeref(sys->manager, state_cube);
     Cudd_RecursiveDeref(sys->manager, sys->all_states);
     sys->all_states = tmp;
}

void set_state_enc(BDDSys* sys, uint ind, NumList* enc)
{
     // mexPrintf("Setting state %d/%d\n", ind, array_len(sys->s_encs));
     while (ind >= array_alloc_size(sys->s_encs))
     {
          // mexPrintf("Growing\n");
          array_grow(sys->s_encs);
     }
     array_set(sys->s_encs, ind, enc);
     array_make_persist(enc);
     uint key = fromBitArray(array_list(enc), array_len(enc));
     // mexPrintf("Setting map %d to %d\n", key, ind);
     array_set(sys->enc_state_map, key, ind);
}

void add_trans(BDDSys* sys, NumList* in_states, NumList* actions, NumList* out_states)
{
     DdNode* trans_BDD = Cudd_ReadLogicZero(sys->manager);
     Cudd_Ref(trans_BDD);
     DdNode* tmp;
     for (int t = 0; t < array_len(in_states); t++)
     {
          uint in_state = array_get(in_states, t);
          uint action = array_get(actions, t);
          uint out_state = array_get(out_states, t);
          // make transition cube
          // mexPrintf("In cube for %d: ", in_state);
          // for (int i = 0; i < array_len(array_get(sys->s_encs, in_state)); i++)
          // {
          //      NumList* enc = array_get(sys->s_encs, in_state);
               // mexPrintf("%d", array_get(enc, i));
          // }
          // mexPrintf("\n");
          DdNode* in_cube = Cudd_bddComputeCube(sys->manager, array_list(sys->s_in_vars), array_list(array_get(sys->s_encs, in_state)), sys->s_var_num);
          // mexPrintf("Made in cube\n");
          // mexPrintf("Out cube for %d: ", out_state);
          // for (int i = 0; i < array_len(array_get(sys->s_encs, out_state)); i++)
          // {
          //      NumList* enc = array_get(sys->s_encs, out_state);
          //      mexPrintf("%d", array_get(enc, i));
          // }
          // mexPrintf("\n");
          DdNode* out_cube = Cudd_bddComputeCube(sys->manager, array_list(sys->s_out_vars), array_list(array_get(sys->s_encs, out_state)), sys->s_var_num);
          // mexPrintf("Made out cube\n");
          // mexPrintf("action cube for %d: ", action);
          // for (int i = 0; i < array_len(array_get(sys->a_encs, action)); i++)
          // {
          //      NumList* enc = array_get(sys->a_encs, action);
               // mexPrintf("%d", array_get(enc, i));
          // }
          // mexPrintf("\n");
          DdNode* action_cube = Cudd_bddComputeCube(sys->manager, array_list(sys->a_vars), array_list(array_get(sys->a_encs, action)), sys->a_var_num);
          // mexPrintf("Made actions cube\n");
          Cudd_Ref(in_cube);
          Cudd_Ref(out_cube);
          Cudd_Ref(action_cube);
          // mexPrintf("Variable cubes computed\n");

          // mexPrintf("Making transition cube\n");
          DdNode* trans_cube = Cudd_bddAnd(sys->manager, action_cube, in_cube);
          Cudd_Ref(trans_cube);
          tmp = Cudd_bddAnd(sys->manager, trans_cube, out_cube);
          Cudd_Ref(tmp);
          Cudd_RecursiveDeref(sys->manager, trans_cube);
          trans_cube = tmp;
          Cudd_RecursiveDeref(sys->manager, in_cube);
          Cudd_RecursiveDeref(sys->manager, out_cube);
          Cudd_RecursiveDeref(sys->manager, action_cube);
          // mexPrintf("Transition cube made\n");

          tmp = Cudd_bddOr(sys->manager, trans_BDD, trans_cube);
          Cudd_Ref(tmp);
          Cudd_RecursiveDeref(sys->manager, trans_cube);
          Cudd_RecursiveDeref(sys->manager, trans_BDD);
          trans_BDD = tmp;
     }
     // append system with transition
     // T | trans
     // mexPrintf("Adding transition to union\n");
     tmp = Cudd_bddOr(sys->manager, sys->trans_sys, trans_BDD);
     Cudd_Ref(tmp);
     Cudd_RecursiveDeref(sys->manager, sys->trans_sys);
     sys->trans_sys = tmp;
     Cudd_RecursiveDeref(sys->manager, trans_BDD);
}

DdNode* get_trans_with_s(BDDSys* sys, NumList* states)
{
     DdNode* tmp;
     // make BDD set X_out
     DdNode* stateBDD = Cudd_ReadLogicZero(sys->manager);
     Cudd_Ref(stateBDD);
     for (int i = 0; i < array_len(states); i++)
     {
          DdNode* state_cube = Cudd_bddComputeCube(sys->manager, array_list(sys->s_out_vars), array_list(array_get(sys->s_encs, array_get(states, i))), sys->s_var_num);
          Cudd_Ref(state_cube);
          tmp = Cudd_bddOr(sys->manager, state_cube, stateBDD);
          Cudd_Ref(tmp);
          Cudd_RecursiveDeref(sys->manager, state_cube);
          Cudd_RecursiveDeref(sys->manager, stateBDD);
          stateBDD = tmp;
     }

     // compute (in, a, state) pairs, T & X_out
     DdNode* in_transitions = Cudd_bddAnd(sys->manager, sys->trans_sys, stateBDD);
     Cudd_Ref(in_transitions);

     // compute (state, a, out) pairs, T & X_in
     tmp = Cudd_bddSwapVariables(sys->manager, stateBDD, array_list(sys->s_out_vars), array_list(sys->s_in_vars), sys->s_var_num);
     Cudd_Ref(tmp);
     Cudd_RecursiveDeref(sys->manager, stateBDD);
     stateBDD = tmp;
     DdNode* out_transitions = Cudd_bddAnd(sys->manager, sys->trans_sys, stateBDD);
     Cudd_Ref(out_transitions);
     Cudd_RecursiveDeref(sys->manager, stateBDD);

     // Make union of all pairs
     DdNode* all_pairs = Cudd_bddOr(sys->manager, in_transitions, out_transitions);
     Cudd_Ref(all_pairs);
     Cudd_RecursiveDeref(sys->manager, in_transitions);
     Cudd_RecursiveDeref(sys->manager, out_transitions);

     return all_pairs;
}

void rm_trans_with_s(BDDSys* sys, NumList* states)
{
     DdNode* tmp;
     DdNode* stateBDD = Cudd_ReadLogicZero(sys->manager);
     Cudd_Ref(stateBDD);
     // make (out)state BDD
     for (int i = 0; i < array_len(states); i++)
     {
          DdNode* state_cube = Cudd_bddComputeCube(sys->manager, array_list(sys->s_out_vars), array_list(array_get(sys->s_encs, array_get(states, i))), sys->s_var_num);
          Cudd_Ref(state_cube);
          tmp = Cudd_bddOr(sys->manager, state_cube, stateBDD);
          Cudd_Ref(tmp);
          Cudd_RecursiveDeref(sys->manager, state_cube);
          Cudd_RecursiveDeref(sys->manager, stateBDD);
          stateBDD = tmp;
     }

     // remove X_out, T = T & !X_out
     tmp = Cudd_bddAnd(sys->manager, sys->trans_sys, Cudd_Not(stateBDD));
     Cudd_Ref(tmp);
     Cudd_RecursiveDeref(sys->manager, sys->trans_sys);
     sys->trans_sys = tmp;

     // make (in)state BDD
     tmp = Cudd_bddSwapVariables(sys->manager, stateBDD, array_list(sys->s_out_vars), array_list(sys->s_in_vars), sys->s_var_num);
     Cudd_Ref(tmp);
     Cudd_RecursiveDeref(sys->manager, stateBDD);
     stateBDD = tmp;

     // remove X_in, T = T & !X_in
     tmp = Cudd_bddAnd(sys->manager, sys->trans_sys, Cudd_Not(stateBDD));
     Cudd_Ref(tmp);
     Cudd_RecursiveDeref(sys->manager, sys->trans_sys);
     Cudd_RecursiveDeref(sys->manager, stateBDD);
     sys->trans_sys = tmp;
}

void add_progress_group(BDDSys* sys, NumList* U, NumList* G)
{
     mexPrintf("Read\n");
     for (int i = 0; i < array_len(U); i++)
          mexPrintf("%d ", array_get(U, i));
     mexPrintf("\n");
     for (int i = 0; i < array_len(G); i++)
          mexPrintf("%d ", array_get(G, i));
     mexPrintf("\n");
     DdNode* U_bdd = makeSet(sys->manager, array_list(sys->a_vars), sys->a_var_num, array_list(U), array_len(U), sys->a_encs);
     DdNode* G_bdd = makeSet(sys->manager, array_list(sys->s_out_vars), sys->s_var_num, array_list(G), array_len(G), sys->s_encs);
     mexPrintf("Made BDDs of groups\n");
     if (U_bdd == NULL)
          mexPrintf("U NULL\n");
     if (G_bdd == NULL)
          mexPrintf("G NULL\n");

     // remove inferior progress groups
     // mexPrintf("Checking inferior groups\n");
     for (int i = 0; i < array_len(sys->pg_U); i++)
     {
          mexPrintf("Checking group %d\n", i);
          if (array_get(sys->pg_U, i) == NULL || array_get(sys->pg_G, i) == NULL)
          {
               mexPrintf("Group %d empty\n", i);
               continue;
          }
          if (Cudd_EquivDC(sys->manager, U_bdd, array_get(sys->pg_U, i), Cudd_Not(array_get(sys->pg_U, i)))
               && Cudd_EquivDC(sys->manager, G_bdd, array_get(sys->pg_G, i), Cudd_Not(array_get(sys->pg_G, i))))
          {
               Cudd_RecursiveDeref(sys->manager, array_get(sys->pg_U, i));
               Cudd_RecursiveDeref(sys->manager, array_get(sys->pg_G, i));
               array_set(sys->pg_U, i, NULL);
               array_set(sys->pg_G, i, NULL);
          }
     }

     // add the progress group
     // mexPrintf("Adding groups\n");
     array_pushback(sys->pg_G, G_bdd);
     mexPrintf("Added G group\n");
     array_pushback(sys->pg_U, U_bdd);
     mexPrintf("Added U group\n");
}

void rm_progress_group(BDDSys* sys, uint index)
{
     if (index < array_len(sys->pg_U) && array_get(sys->pg_U, index) != NULL)
     {
          Cudd_RecursiveDeref(sys->manager, array_get(sys->pg_U, index));
          Cudd_RecursiveDeref(sys->manager, array_get(sys->pg_G, index));
          array_set(sys->pg_U, index, NULL);
          array_set(sys->pg_G, index, NULL);
     }
}

int has_superior_pg(BDDSys* sys, NumList* U, NumList* G)
{
     DdNode* U_bdd = (DdNode*)makeSet(sys->manager, array_list(sys->a_vars), sys->a_var_num, array_list(U), array_len(U), sys->a_encs);
     DdNode* G_bdd = (DdNode*)makeSet(sys->manager, array_list(sys->s_out_vars), sys->s_var_num, array_list(G), array_len(G), sys->s_encs);

     for (int i = 0; i < array_len(sys->pg_U); i++)
     {
          if (array_len(sys->pg_U) == 0
               || array_len(sys->pg_G) == 0
               || array_get(sys->pg_U, i) == NULL
               || array_get(sys->pg_G, i) == NULL)
               continue;

          if (Cudd_EquivDC(sys->manager, array_get(sys->pg_U, i), U_bdd, Cudd_Not(U_bdd))
               && Cudd_EquivDC(sys->manager, array_get(sys->pg_G, i), G_bdd, Cudd_Not(G_bdd)))
          {
               Cudd_RecursiveDeref(sys->manager, U_bdd);
               Cudd_RecursiveDeref(sys->manager, G_bdd);
               return 1;
          }
     }

     Cudd_RecursiveDeref(sys->manager, U_bdd);
     Cudd_RecursiveDeref(sys->manager, G_bdd);
     return 0;
}

void add_to_pg(BDDSys* sys, NumList* indices, NumList* states)
{
     DdNode* G_add = makeSet(sys->manager, array_list(sys->s_out_vars), sys->s_var_num, array_list(states), array_len(states), sys->s_encs);

     for (int i = 0; i < array_len(indices); i++)
     {
          uint j = array_get(indices, i);
          // mexPrintf("Editing group %d\n", j);
          if (array_get(sys->pg_G, j) == NULL)
          {
               // TODO: Do something or error?
               continue;
          }
          DdNode* tmp = Cudd_bddOr(sys->manager, array_get(sys->pg_G, j), G_add);
          Cudd_Ref(tmp);
          Cudd_RecursiveDeref(sys->manager, array_get(sys->pg_G, j));
          array_set(sys->pg_G, j, tmp);
     }
     Cudd_RecursiveDeref(sys->manager, G_add);
}

NumList* is_member_of_pg(BDDSys* sys, NumList* states)
{
     uint ind_num = 0;
     // mexPrintf("Making the index list\n");
     // mexEvalString("drawnow;");
     array_init(indices, NumList, 0, MAX(array_len(sys->pg_G), 1));

     // mexPrintf("Making the state BDD\n");
     // mexEvalString("drawnow;");
     DdNode* state_bdd = makeSet(sys->manager, array_list(sys->s_out_vars), sys->s_var_num, array_list(states), array_len(states), sys->s_encs);
     uint** parts = read_out_states(sys, state_bdd);
     uint len = *(parts[0]);
     // mexPrintf("Checking inputs\n");
     // for (int i = 0; i < len; i++)
     //      mexPrintf("%d ", parts[1][i]);
     // mexPrintf("\n");

     // mexPrintf("Checking progress groups\n");
     // mexEvalString("drawnow;");
     for (int i = 0; i < array_len(sys->pg_U); i++)
     {
          // mexPrintf("Checking inside\n");
          // mexEvalString("drawnow;");
          if (array_get(sys->pg_G, i) != NULL
               && Cudd_EquivDC(sys->manager, array_get(sys->pg_G, i), state_bdd, Cudd_Not(state_bdd)))
          {
               // mexPrintf("Found something\n");
               // mexEvalString("drawnow;");
               array_pushback(indices, i);
          }
     }
     mxFree(parts[1]);
     mxFree(parts);
     Cudd_RecursiveDeref(sys->manager, state_bdd);
     // mexPrintf("Done\n");
     return indices;
}

uint get_num_of_trans(BDDSys* sys, DdNode* bdd)
{
     return Cudd_CountMinterm(sys->manager, bdd, 2*sys->s_var_num + sys->a_var_num);
}

uint** read_transitions(BDDSys* sys, DdNode* T)
{
     uint max_trans_num = get_num_of_trans(sys, T);

     uint* in_states = mxMalloc(sizeof(uint)*max_trans_num);
     uint* actions = mxMalloc(sizeof(uint)*max_trans_num);
     uint* out_states = mxMalloc(sizeof(uint)*max_trans_num);
     uint trans_num = 0;

     uint** outputs = mxMalloc(sizeof(uint*)*4);
     outputs[0] = &trans_num;
     outputs[1] = actions;
     outputs[2] = in_states;
     outputs[3] = out_states;

     NumList** intervals = mxMalloc(sizeof(NumList*)*3);
     intervals[0] = sys->s_in_inds;
     intervals[1] = sys->a_inds;
     intervals[2] = sys->s_out_inds;
     array_init(cube_inds, NumList, 0, sys->a_var_num + 2*sys->s_var_num);
     for (int i = 0; i < sys->a_var_num; i++)
     {
          array_pushback(cube_inds, array_get(intervals[1], i));
     }
     for (int i = 0; i < sys->s_var_num; i++)
     {
          array_pushback(cube_inds, array_get(intervals[0], i));
     }
     for (int i = 0; i < sys->s_var_num; i++)
     {
          array_pushback(cube_inds, array_get(intervals[2], i));
     }

     DdGen* gen;
     int* cube;
     CUDD_VALUE_TYPE value;
     Cudd_ForeachCube(sys->manager, T, gen, cube, value)
     {
          read_in(intervals, 3, cube_inds, cube, 0, &(outputs[1]), &trans_num);
     }
     // convert keys to state
     for (int i = 0; i < trans_num; i++)
     {
          outputs[1][i] = array_get(sys->enc_state_map, outputs[1][i]);
          outputs[3][i] = array_get(sys->enc_state_map, outputs[3][i]);
     }

     mxFree(intervals);
     array_free(cube_inds);

     return outputs;
}

uint** read_in_states(BDDSys* sys, DdNode* T)
{
     uint max_trans_num = array_len(sys->s_encs);
     uint state_num = 0;

     uint* in_states = mxMalloc(sizeof(uint)*max_trans_num);
     uint** outputs = mxMalloc(sizeof(uint*)*2);
     outputs[0] = &state_num;
     outputs[1] = in_states;

     DdGen* gen;
     int* cube;
     CUDD_VALUE_TYPE value;
     Cudd_ForeachCube(sys->manager, T, gen, cube, value)
     {
          read_in(&sys->s_in_inds, 1, sys->s_in_inds, cube, 0, &(outputs[1]), &state_num);
     }
     // convert keys to state
     for (int i = 0; i < state_num; i++)
          outputs[1][i] = array_get(sys->enc_state_map, outputs[1][i]);

     return outputs;
}

uint** read_out_states(BDDSys* sys, DdNode* T)
{
     uint max_state_num = array_len(sys->s_encs);
     uint state_num = 0;

     uint* out_states = mxMalloc(sizeof(uint)*max_state_num);
     uint** outputs = mxMalloc(sizeof(uint*)*2);
     outputs[0] = &state_num;
     outputs[1] = out_states;

     DdGen* gen;
     int* cube;
     CUDD_VALUE_TYPE value;
     Cudd_ForeachCube(sys->manager, T, gen, cube, value)
     {
          read_in(&sys->s_out_inds, 1, sys->s_out_inds, cube, 0, &out_states, &state_num);
     }
     // convert keys to state
     for (int i = 0; i < state_num; i++)
     {
          // mexPrintf("Converting %d to %d\n", outputs[1][i], array_get(sys->enc_state_map, outputs[1][i]));
          outputs[1][i] = array_get(sys->enc_state_map, outputs[1][i]);
     }
     return outputs;
}

uint** read_actions(BDDSys* sys, DdNode* T)
{
     uint max_trans_num = array_len(sys->a_encs);
     uint action_num = 0;

     uint* a_states = mxMalloc(sizeof(uint)*max_trans_num);
     uint** outputs = mxMalloc(sizeof(uint*)*2);

     DdGen* gen;
     int* cube;
     CUDD_VALUE_TYPE value;
     int cube_count = 0;
     Cudd_ForeachCube(sys->manager, T, gen, cube, value)
     {
          read_in(&sys->a_inds, 1, sys->a_inds, cube, 0, &a_states, &action_num);
     }
     outputs[0] = &action_num;
     outputs[1] = a_states;

     return outputs;
}

void read_in(NumList** intervals, uint interv_num, NumList* cube_inds, int* cube, uint cube_pos, uint** outputs, uint* output_num)
{
     if (cube_pos >= array_len(cube_inds))
     {
          // mexPrintf("\n");
          // read action, in_state, out_state
          for (int i = 0; i < interv_num; i++)
          {
               NumList* interv = intervals[i];

               array_init(interval_cube, NumList, array_len(interv), 1);
               // mexPrintf("Cube: ");
               for (int j = 0; j < array_len(interv); j++)
               {
                    // mexPrintf("%d", cube[array_get(interv, j)]);
                    array_set(interval_cube, j, (uint)cube[array_get(interv, j)]);
               }
               outputs[i][*output_num] = fromBitArray(array_list(interval_cube), array_len(interval_cube));
               // mexPrintf("=%d\n", outputs[i][*output_num]);
               array_free(interval_cube);
          }
          (*output_num)++;
          return;
     }

     uint index = array_get(cube_inds, cube_pos);
     // mexPrintf("(%d,%d)", array_get(cube_inds, cube_pos), cube[array_get(cube_inds, cube_pos)]);

     if (cube[index] == 2)
     {
          cube[index] = 1;
          read_in(intervals, interv_num, cube_inds, cube, cube_pos+1, outputs, output_num);
          cube[index] = 0;
          read_in(intervals, interv_num, cube_inds, cube, cube_pos+1, outputs, output_num);
          cube[index] = 2;
     }
     else
          read_in(intervals, interv_num, cube_inds, cube, cube_pos+1, outputs, output_num);
}

uint fromBitArray(uint* arr, uint len)
{
     uint val = 0;
     for (int i = len-1; i >= 0; i--)
     {
          val = (val << 1) | arr[i];
     }
     return val;
}

uint fromBitArray_int(int* arr, uint len)
{
     uint val = 0;
     for (int i = len-1; i >= 0; i--)
     {
          val = (val << 1) | arr[i];
     }
     return val;
}

DdNode* makeSet(DdManager* manager, DdNode** vars, int var_num, int* inds, int n, EncList* encodings)
{
    DdNode* set = Cudd_ReadLogicZero(manager);
    Cudd_Ref(set);

    for (int i = 0; i < n; i++)
    {
        int* enc = (int*)array_list(array_get(encodings, inds[i]));

        DdNode* cube = Cudd_bddComputeCube(manager, vars, enc, var_num);
        Cudd_Ref(cube);
        DdNode* tmp = Cudd_bddOr(manager, set, cube);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, set);
        Cudd_RecursiveDeref(manager, cube);
        set = tmp;
    }

    return set;
}

DdNode* read_BDD_from_mxarray(BDDSys* sys, const mxArray* array, char setting)
{
     uint len = mxGetNumberOfElements(array);
     double* input = mxGetPr(array);
     uint* inds = mxMalloc(sizeof(uint)*len);

     for (int i = 0; i < len; i++)
     {
          inds[i] = (uint)input[i] - 1;
     }
     DdNode* bdd;
     if (setting == 's')
     {
          bdd = makeSet(sys->manager, array_list(sys->s_out_vars), sys->s_var_num, inds, len, sys->s_encs);
     }
     else if (setting == 'a')
          bdd = makeSet(sys->manager, array_list(sys->a_vars), sys->a_var_num, inds, len, sys->a_encs);
     else
     {
          mxFree(inds);
          return NULL;
     }
     mxFree(inds);

     return bdd;
}

void write_BDD_to_mxarray(BDDSys* sys, DdNode* bdd, mxArray** array, char setting)
{
     uint** set_parts;
     if (setting == 's')
          set_parts = read_out_states(sys, bdd);
     else if (setting == 'a')
          set_parts = read_actions(sys, bdd);
     uint set_num = *(set_parts[0]);
     *array = mxCreateDoubleMatrix(set_num, 1, mxREAL);
     double* ptr = mxGetPr(*array);
     for (int i = 0; i < set_num; i++)
          ptr[i] = (double)set_parts[1][i] + 1;
     mxFree(set_parts[1]);
     mxFree(set_parts);
}
