#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include "mex.h"
#include "cudd.h"
#include "system_types.h"

#define deref(system) Cudd_RecursiveDeref(get_mgr(sys), system);

#define WIN_SET 1
#define WIN_CANDIDATE_SET 2
#define WIN_CANDIDATE_CONT 3

DdNode** win_until_and_always(BDDSys* sys, DdNode* A, DdNode* B, DdNode* Z,  char quant, int mode);
DdNode** win_intermediate(BDDSys* sys, DdNode* A, DdNode* Z, DdNode* B,
                          DdNode** C, int C_num, char quant, int mode);



//  uint fromBitArray_int_i(int* arr, uint len)
//  {
//       uint val = 0;
//       for (int i = len-1; i >= 0; i--)
//       {
//            val = (val << 1) | arr[i];
//       }
//       return val;
//  }
//
//  void read_in_i(NumList** intervals, uint interv_num, NumList* cube_inds, int* cube, uint cube_pos, uint** outputs, uint* output_num)
//  {
//       if (cube_pos >= array_len(cube_inds))
//      {
//           // read action, in_state, out_state
//           for (int i = 0; i < interv_num; i++)
//           {
//                NumList* interv = intervals[i];
//                array_init(interval_cube, NumList, array_len(interv), 1);
//                for (int j = 0; j < array_len(interv); j++)
//                {
//                     array_set(interval_cube, j, cube[array_get(interv, j)]);
//                }
//                outputs[i][*output_num] = fromBitArray_int_i(array_list(interval_cube), array_len(interval_cube));
//                array_free(interval_cube);
//           }
//           (*output_num)++;
//           return;
//      }
//
//      uint index = array_get(cube_inds, cube_pos);
//
//      if (cube[index] == 2)
//      {
//           cube[index] = 1;
//           read_in_i(intervals, interv_num, cube_inds, cube, cube_pos+1, outputs, output_num);
//           cube[index] = 0;
//           read_in_i(intervals, interv_num, cube_inds, cube, cube_pos+1, outputs, output_num);
//           cube[index] = 2;
//      }
//      else
//           read_in_i(intervals, interv_num, cube_inds, cube, cube_pos+1, outputs, output_num);
// }
//
//  uint** read_in_states_i(BDDSys* sys, DdNode* T)
//  {
//       uint max_trans_num = array_len(sys->s_encs);
//      uint state_num = 0;
//
//      uint* in_states = mxMalloc(sizeof(uint)*max_trans_num);
//      uint** outputs = mxMalloc(sizeof(uint*)*2);
//      outputs[0] = &state_num;
//      outputs[1] = in_states;
//
//      DdGen* gen;
//      int* cube;
//      CUDD_VALUE_TYPE value;
//      Cudd_ForeachCube(sys->manager, T, gen, cube, value)
//      {
//           read_in_i(&sys->s_in_inds, 1, sys->s_in_inds, cube, 0, &(outputs[1]), &state_num);
//      }
//      // convert keys to state
//      for (int i = 0; i < state_num; i++)
//           outputs[1][i] = array_get(sys->enc_state_map, outputs[1][i]);
//
//      return outputs;
//  }
//
//  uint** read_out_states_i(BDDSys* sys, DdNode* T)
//  {
//       uint max_state_num = array_len(sys->s_encs);
//      uint state_num = 0;
//
//      uint* out_states = mxMalloc(sizeof(uint)*max_state_num);
//      uint** outputs = mxMalloc(sizeof(uint*)*2);
//      outputs[0] = &state_num;
//      outputs[1] = out_states;
//
//      DdGen* gen;
//      int* cube;
//      CUDD_VALUE_TYPE value;
//      Cudd_ForeachCube(sys->manager, T, gen, cube, value)
//      {
//           read_in_i(&sys->s_out_inds, 1, sys->s_out_inds, cube, 0, &(outputs[1]), &state_num);
//      }
//      // convert keys to state
//      for (int i = 0; i < state_num; i++)
//           outputs[1][i] = array_get(sys->enc_state_map, outputs[1][i]);
//
//      return outputs;
//  }
//
//  uint** read_actions_i(BDDSys* sys, DdNode* T)
//  {
//       uint max_trans_num = array_len(sys->a_encs);
//      uint action_num = 0;
//
//      uint* a_states = mxMalloc(sizeof(uint)*max_trans_num);
//      uint** outputs = mxMalloc(sizeof(uint*)*2);
//      outputs[0] = &action_num;
//      outputs[1] = a_states;
//
//      DdGen* gen;
//      int* cube;
//      CUDD_VALUE_TYPE value;
//      int cube_count = 0;
//      Cudd_ForeachCube(sys->manager, T, gen, cube, value)
//      {
//           read_in_i(&sys->a_inds, 1, sys->a_inds, cube, 0, &a_states, &action_num);
//      }
//
//      return outputs;
//  }


DdNode* pre(BDDSys* sys, DdNode* X_prime, DdNode* A, char quant1, char quant2)
{
    DdManager* manager = get_mgr(sys);
    DdNode* T = sys->trans_sys;
    DdNode* valid_transitions;
    DdNode* tmp;

    DdNode* outCube = Cudd_bddComputeCube(manager, array_list(sys->mgr->s_out_vars), NULL, sys->mgr->s_var_num);
    Cudd_Ref(outCube);
    DdNode* aCube = Cudd_bddComputeCube(manager, array_list(sys->mgr->a_vars), NULL, sys->mgr->a_var_num);
    Cudd_Ref(aCube);

    // Cudd_DebugCheck(sys->manager);
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
        DdNode* state_a_pairs_not_X = Cudd_bddAndAbstract(manager, T, Cudd_Not(X_prime), outCube);
        Cudd_Ref(state_a_pairs_not_X);

        DdNode* state_a_pairs_X = Cudd_bddExistAbstract(manager, T, outCube);
        Cudd_Ref(state_a_pairs_X);

        valid_transitions = Cudd_bddAnd(manager, state_a_pairs_X, Cudd_Not(state_a_pairs_not_X));
        Cudd_Ref(valid_transitions);
        Cudd_RecursiveDeref(manager, state_a_pairs_X);
        Cudd_RecursiveDeref(manager, state_a_pairs_not_X);

        tmp = Cudd_bddAndAbstract(manager, valid_transitions, A, aCube);
        Cudd_Ref(tmp);
        deref(valid_transitions);
        valid_transitions = tmp;
    }
    else if(quant1 == 'a' && quant2 == 'a')
    {
        DdNode* state_a_pairs_not_X = Cudd_bddAndAbstract(manager, T, Cudd_Not(X_prime), outCube);
        Cudd_Ref(state_a_pairs_not_X);

        DdNode* state_a_pairs_X = Cudd_bddExistAbstract(manager, T, outCube);
        Cudd_Ref(state_a_pairs_X);

        valid_transitions = Cudd_bddAnd(manager, state_a_pairs_X, Cudd_Not(state_a_pairs_not_X));
        Cudd_Ref(valid_transitions);
        Cudd_RecursiveDeref(manager, state_a_pairs_X);
        Cudd_RecursiveDeref(manager, state_a_pairs_not_X);

        tmp = Cudd_bddOr(manager, valid_transitions, Cudd_Not(A));
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, valid_transitions);
        valid_transitions = tmp;

        tmp = Cudd_bddUnivAbstract(manager, valid_transitions, aCube);
        Cudd_Ref(tmp);
        deref(valid_transitions);
        valid_transitions = tmp;
    }
    else if (quant1 == 'a' && quant2 == 'e')
    {
        tmp = Cudd_bddAndAbstract(manager, T, X_prime, outCube);
        Cudd_Ref(tmp);
        valid_transitions = tmp;

        tmp = Cudd_bddOr(manager, valid_transitions, Cudd_Not(A));
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, valid_transitions);
        valid_transitions = tmp;

        tmp = Cudd_bddUnivAbstract(manager, valid_transitions, aCube);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, valid_transitions);
        valid_transitions = tmp;
    }
    else
    {
        printf("Invalid quantifiers! %c %c\n", quant1, quant2);
        return NULL;
    }

    Cudd_RecursiveDeref(manager, outCube);
    Cudd_RecursiveDeref(manager, aCube);

    return valid_transitions;
}

DdNode** inv(BDDSys* sys, DdNode* Z, DdNode* B, DdNode* U, DdNode* G,
            char quant, int mode)
{
    DdManager* manager = get_mgr(sys);
    // calculate Y_0 = (G & B) \ Z
    DdNode* tmp;
    DdNode* Y = Cudd_bddAnd(manager, G, B);
    Cudd_Ref(Y);
    tmp = Cudd_bddAnd(manager, Y, Cudd_Not(Z));
    Cudd_Ref(tmp);
    deref(Y);
    Y = tmp;
    DdNode* Cw = NULL;

    DdNode* quick_test = pre(sys, Z, U, quant, 'e');
    tmp = Cudd_bddSwapVariables(manager, quick_test, array_list(sys->mgr->s_out_vars), array_list(sys->mgr->s_in_vars), sys->mgr->s_var_num);
    Cudd_Ref(tmp);
    deref(quick_test);
    quick_test = tmp;
    tmp = Cudd_bddAnd(manager, quick_test, Y);
    Cudd_Ref(tmp);
    deref(quick_test);
    quick_test = tmp;

    // no reach to Z
    if (Cudd_EquivDC(manager, quick_test, Cudd_ReadLogicZero(manager), Cudd_ReadLogicZero(manager)))
    {
        deref(quick_test);
        deref(Y);
        Y = Cudd_ReadLogicZero(manager);
        Cudd_Ref(Y);
        if (mode >= WIN_CANDIDATE_SET)
        {
            Cw = Cudd_ReadLogicZero(manager);
            Cudd_Ref(Cw);
        }
        DdNode** out = mxMalloc(sizeof(DdNode*)*mode);
        if (mode >= WIN_SET)
        {
            out[0] = Y;
            if (mode >= WIN_CANDIDATE_SET)
                out[1] = Cw;
        }

        return out;
    }
    deref(quick_test);

    if (mode >= WIN_CANDIDATE_SET) // get candidate space
    {
        Cw = Cudd_bddAnd(manager, Y, Cudd_Not(Z));
        Cudd_Ref(Cw);
    }
    // printf("Y:\n");
    // PRINT_OUT(manager, Y, sys->a_var_num, sys->s_var_num, map, 0);

    // iterate with Y_k+1 = Y_k & Pre_ea(Y_k & Z)
    DdNode* old_Y = Cudd_ReadLogicZero(manager);
    Cudd_Ref(old_Y);
    while (!Cudd_EquivDC(manager, Y, old_Y, Cudd_ReadLogicZero(manager)))
    {
        deref(old_Y);
        old_Y = Y;

        tmp = Cudd_bddOr(manager, old_Y, Z);
        Cudd_Ref(tmp);
        DdNode* pre_set = pre(sys, tmp, U, quant, 'a');
        deref(tmp);

        tmp = Cudd_bddSwapVariables(manager, pre_set, array_list(sys->mgr->s_out_vars), array_list(sys->mgr->s_in_vars), sys->mgr->s_var_num);
        Cudd_Ref(tmp);
        deref(pre_set);
        pre_set = tmp;
        Y = Cudd_bddAnd(manager, old_Y, pre_set);
        Cudd_Ref(Y);
        deref(pre_set);
    }
    deref(old_Y);

    DdNode** out = mxMalloc(sizeof(DdNode*)*mode);
    if (mode >= WIN_SET)
    {
        out[0] = Y;
        if (mode >= WIN_CANDIDATE_SET)
            out[1] = Cw;
    }

    return out;
}

DdNode** PGpre(BDDSys* sys, DdNode* Z, DdNode* B, char quant, int mode)
{
    DdManager* manager = get_mgr(sys);
    DdNode* U;
    DdNode* G;
    DdNode* PG_set = Cudd_bddAnd(manager, Z, Cudd_ReadOne(manager));
    Cudd_Ref(PG_set);
    DdNode* inv_set;
    DdNode* tmp;
    DdNode* Cw;
    if (mode >= WIN_CANDIDATE_SET)
    {
        Cw = Cudd_ReadLogicZero(manager);
        Cudd_Ref(Cw);
    }
    for (int i = 0; i < array_len(sys->pg_U); i++)
    {
        if (array_get(sys->pg_U, i) == NULL)
            continue;
        if (quant == 'a' && !Cudd_EquivDC(manager, array_get(sys->pg_U, i), sys->all_actions, Cudd_Not(sys->all_actions)))
        {
            mexPrintf("This happens! Progress group skipped\n");
            continue;
        }

        U = array_get(sys->pg_U, i);
        G = array_get(sys->pg_G, i);
        DdNode** in = inv(sys, PG_set, B, U, G, quant, mode);
        if (mode >= WIN_SET)
        {
            inv_set = in[0];
            if (mode >= WIN_CANDIDATE_SET)
            {
                deref(Cw);
                Cw = in[1];
            }
        }
        mxFree(in);
        tmp = Cudd_bddOr(manager, inv_set, PG_set);
        Cudd_Ref(tmp);
        deref(inv_set);
        deref(PG_set);
        PG_set = tmp;
    }

    DdNode** out = mxMalloc(sizeof(DdNode*)*mode);

    if (mode >= WIN_SET)
    {
        out[0] = PG_set;
        if (mode >= WIN_CANDIDATE_SET)
            out[1] = Cw;
    }

    return out;
}

DdNode** win_until(BDDSys* sys, DdNode* Z, DdNode* B, char quant, int mode)
{
    DdManager* manager = get_mgr(sys);
    // calculates Win_q,a(B U Z)
    DdNode* PGpre_set;
    DdNode* pre_set;
    DdNode* tmp;
    DdNode* Cw;
    if (mode >= WIN_CANDIDATE_SET)
    {
        Cw = Cudd_ReadLogicZero(manager);
        Cudd_Ref(Cw);
    }
    DdNode* one = Cudd_ReadOne(manager);
    Cudd_Ref(one);
    DdNode* X = Cudd_ReadLogicZero(manager);
    Cudd_Ref(X);
    DdNode* old_X = Cudd_ReadOne(manager);
    Cudd_Ref(old_X);
    while(!Cudd_EquivDC(manager, X, old_X, Cudd_ReadLogicZero(manager)))
    {
        deref(old_X);
        old_X = X;
        DdNode** in = PGpre(sys, X, B, quant, mode);
        if (mode >= WIN_SET)
        {
            PGpre_set = in[0];
            if (mode >= WIN_CANDIDATE_SET)
            {
                deref(Cw);
                Cw = in[1];
            }
        }
        mxFree(in);
        pre_set = pre(sys, X, one, quant, 'a');

        tmp = Cudd_bddSwapVariables(manager, pre_set, array_list(sys->mgr->s_in_vars), array_list(sys->mgr->s_out_vars), sys->mgr->s_var_num);
        Cudd_Ref(tmp);
        deref(pre_set);
        pre_set = tmp;
        tmp = Cudd_bddAnd(manager, B, pre_set);
        Cudd_Ref(tmp);
        X = tmp;
        tmp = Cudd_bddOr(manager, PGpre_set, X);
        Cudd_Ref(tmp);
        deref(X);
        X = tmp;
        tmp = Cudd_bddOr(manager, Z, X);
        Cudd_Ref(tmp);
        deref(X);
        X = tmp;

        deref(PGpre_set);
        deref(pre_set);
    }

    if (mode >= WIN_CANDIDATE_SET)
    {
        DdNode* XPre = pre(sys, X, one, quant, 'a');
        tmp = Cudd_bddSwapVariables(manager, XPre, array_list(sys->mgr->s_out_vars), array_list(sys->mgr->s_in_vars), sys->mgr->s_var_num);
        Cudd_Ref(tmp);
        deref(XPre);
        XPre = tmp;
        tmp = Cudd_bddAnd(manager, XPre, Cudd_Not(X));
        Cudd_Ref(tmp);
        deref(XPre);
        XPre = tmp;
        tmp = Cudd_bddOr(manager, Cw, XPre);
        Cudd_Ref(tmp);
        deref(Cw);
        deref(XPre);
        Cw = tmp;
    }

    deref(one);
    deref(old_X);

    DdNode** out = mxMalloc(sizeof(DdNode*)*mode);
    if (mode >= WIN_SET)
    {
        out[0] = X;
        if (mode >= WIN_CANDIDATE_SET)
            out[1] = Cw;
    }
    return out;
}

DdNode** win_until_and_always(BDDSys* sys, DdNode* A, DdNode* B, DdNode* Z, char quant, int mode)
{
    DdManager* manager = get_mgr(sys);
    DdNode* always_A = Cudd_ReadOne(manager);
    Cudd_Ref(always_A);
    DdNode* win_until_set;
    DdNode* B_always;
    DdNode* Z_always;
    DdNode* Cw;
    if (mode >= WIN_CANDIDATE_SET)
    {
        Cw = Cudd_ReadLogicZero(manager);
        Cudd_Ref(Cw);
    }
    DdNode* tmp;

    if (!Cudd_EquivDC(manager, A, sys->all_states, Cudd_Not(sys->all_states)))
    {
        DdNode* zero = Cudd_ReadLogicZero(manager);
        Cudd_Ref(zero);
        DdNode* one = Cudd_ReadOne(manager);
        Cudd_Ref(one);
        deref(always_A);

        DdNode** in = win_intermediate(sys, sys->all_states, A, zero, &one, 1, quant, mode);
        Cudd_RecursiveDeref(manager, zero);
        Cudd_RecursiveDeref(manager, one);
        if (mode >= WIN_SET)
        {
            always_A = in[0];
            if (mode >= WIN_CANDIDATE_SET)
            {
                deref(Cw);
                Cw = in[1];
            }
        }
        mxFree(in);
    }

    B_always = Cudd_bddAnd(manager, always_A, B);
    Cudd_Ref(B_always);
    Z_always = Cudd_bddAnd(manager, always_A, Z);
    Cudd_Ref(Z_always);
    deref(always_A);

    DdNode** in = win_until(sys, Z_always, B_always, quant, mode);
    Cudd_RecursiveDeref(manager, B_always);
    Cudd_RecursiveDeref(manager, Z_always);
    if (mode >= WIN_SET)
    {
        win_until_set = in[0];
        if (mode >= WIN_CANDIDATE_SET)
        {
            tmp = Cudd_bddOr(manager, Cw, in[1]);
            Cudd_Ref(tmp);
            deref(Cw);
            deref(in[1]);
            Cw = tmp;
        }
    }
    mxFree(in);

    DdNode** out = mxMalloc(sizeof(DdNode*)*mode);
    if (mode >= WIN_SET)
    {
        out[0] = win_until_set;
        if (mode >= WIN_CANDIDATE_SET)
            out[1] = Cw;
    }
    return out;
}

DdNode** win_intermediate(BDDSys* sys, DdNode* A, DdNode* B, DdNode* Z,
                         DdNode** C, int C_num, char quant, int mode)
{
    DdManager* manager = get_mgr(sys);
    DdNode* W;
    DdNode* old_W;
    DdNode** Z_iter = (DdNode**)mxMalloc(sizeof(DdNode*)*C_num);
    DdNode* tmp;
    DdNode* pre_set;
    DdNode* one = Cudd_ReadOne(manager);
    Cudd_Ref(one);
    DdNode* win_until_set;
    DdNode* Cw;
    if (mode >= WIN_CANDIDATE_SET)
    {
        Cw = Cudd_ReadLogicZero(manager);
        Cudd_Ref(Cw);
    }
    DdNode* W_1;

    W = Cudd_ReadOne(manager);
    Cudd_Ref(W);
    old_W = Cudd_ReadLogicZero(manager);
    Cudd_Ref(old_W);
    int first_set = 1;
    while (!Cudd_EquivDC(manager, W, old_W, Cudd_ReadLogicZero(manager)))
    {
        deref(old_W);
        old_W = W;
        pre_set = pre(sys, W, one, quant, 'a');
        tmp = Cudd_bddSwapVariables(manager, pre_set, array_list(sys->mgr->s_in_vars), array_list(sys->mgr->s_out_vars), sys->mgr->s_var_num);
        Cudd_Ref(tmp);
        deref(pre_set);
        pre_set = tmp;

        DdNode* intermediate = Cudd_bddAnd(manager, B, pre_set);
        Cudd_Ref(intermediate);
        deref(pre_set);

        for (int i = 0; i < C_num; i++)
        {
            tmp = Cudd_bddAnd(manager, C[i], intermediate);
            Cudd_Ref(tmp);
            Z_iter[i] = tmp;
            tmp = Cudd_bddOr(manager, Z_iter[i], Z);
            Cudd_Ref(tmp);
            deref(Z_iter[i]);
            Z_iter[i] = tmp;

            DdNode** in = win_until_and_always(sys, A, B, Z_iter[i], quant, first_set ? mode : WIN_SET);
            if (mode >= WIN_SET)
            {
                win_until_set = in[0];
                if (first_set && mode >= WIN_CANDIDATE_SET)
                {
                    tmp = Cudd_bddOr(manager, Cw, in[1]);
                    Cudd_Ref(tmp);
                    deref(Cw);
                    deref(in[1]);
                    Cw = tmp;
                }
            }
            mxFree(in);
            tmp = Cudd_bddAnd(manager, W, win_until_set);
            Cudd_Ref(tmp);
            if (i > 0)
                deref(W);
            W = tmp;
            deref(win_until_set);
            deref(Z_iter[i]);
        }
        deref(intermediate);

        if (mode >= WIN_CANDIDATE_SET && first_set)
        {
            W_1 = Cudd_bddAnd(manager, W, Cudd_ReadOne(manager));
            Cudd_Ref(W_1);
            first_set = 0;
        }

    }
    mxFree(Z_iter);
    deref(one);
    deref(old_W);

    DdNode** out = mxMalloc(sizeof(DdNode*)*mode);
    if (mode >= WIN_SET)
    {
        out[0] = W;
        if (mode >= WIN_CANDIDATE_SET)
        {
            DdNode* W_set = Cudd_bddAnd(manager, W_1, Cudd_Not(W));
            Cudd_Ref(W_set);
            deref(W_1);
            tmp = Cudd_bddOr(manager, Cw, W_set);
            Cudd_Ref(tmp);
            deref(Cw);
            deref(W_set);
            Cw = tmp;
            out[1] = Cw;
        }
    }

    return out;
}

DdNode** win_primal(BDDSys* sys, DdNode* A, DdNode* B, DdNode** C, int C_num,
                   char quant1, char quant2, DdNode* head_start, int mode)
{
    DdManager* manager = get_mgr(sys);
    int C_created = 0;
    if (C_num <= 0 || C == NULL)
    {
        C_num = 1;
        C = mxMalloc(sizeof(DdNode*));
        C[0] = Cudd_ReadOne(manager);
        Cudd_Ref(C[0]);
        C_created = 1;
    }

    // TODO: dualization

    DdNode* V;
    if (head_start == NULL)
    {
        V = Cudd_ReadLogicZero(manager);
        Cudd_Ref(V);
    }
    else
    {
        V = Cudd_bddAnd(manager, head_start, Cudd_ReadOne(manager));
        Cudd_Ref(V);
    }
    DdNode* old_V = Cudd_ReadOne(manager);
    Cudd_Ref(old_V);
    DdNode* one = Cudd_ReadOne(manager);
    Cudd_Ref(one);
    DdNode* pre_set;
    DdNode* pg_pre_set;
    DdNode* Z;
    DdNode* tmp;
    DdNode* Cw;
    if (mode >= WIN_CANDIDATE_SET)
    {
        Cw = Cudd_ReadLogicZero(manager);
        Cudd_Ref(Cw);
    }

    if (Cudd_EquivDC(manager, A, Cudd_ReadLogicZero(manager), Cudd_ReadLogicZero(manager)))
         A = one;
    if (Cudd_EquivDC(manager, B, Cudd_ReadLogicZero(manager), Cudd_ReadLogicZero(manager)))
         B = one;
    int first_set = 1;
    while (!Cudd_EquivDC(manager, V, old_V, Cudd_ReadLogicZero(manager)))
    {
        deref(old_V);
        old_V = V;
        pre_set = pre(sys, V, one, quant1, 'a');
        tmp = Cudd_bddSwapVariables(manager, pre_set, array_list(sys->mgr->s_in_vars), array_list(sys->mgr->s_out_vars), sys->mgr->s_var_num);
        Cudd_Ref(tmp);
        deref(pre_set);
        pre_set = tmp;
        DdNode** PGpre_in = PGpre(sys, V, A, quant1, WIN_SET);
        pg_pre_set = PGpre_in[0];
        mxFree(PGpre_in);
        Z = Cudd_bddOr(manager, pre_set, pg_pre_set);
        Cudd_Ref(Z);
        deref(pre_set);
        deref(pg_pre_set);

        DdNode** in = win_intermediate(sys, A, B, Z, C, C_num, quant1, first_set ? mode : WIN_SET);
        if (mode >= WIN_SET)
        {
            V = in[0];
            if (first_set && mode >= WIN_CANDIDATE_SET)
            {
                deref(Cw);
                Cw = in[1];
                first_set = 0;
            }
        }
        deref(Z);
        mxFree(in);
    }
    deref(old_V);

    if (mode >= WIN_CANDIDATE_SET)
    {
        pre_set = pre(sys, V, one, quant1, 'e');
        tmp = Cudd_bddSwapVariables(manager, pre_set, array_list(sys->mgr->s_out_vars), array_list(sys->mgr->s_in_vars), sys->mgr->s_var_num);
        Cudd_Ref(tmp);
        deref(pre_set);
        pre_set = tmp;
        tmp = Cudd_bddAnd(manager, pre_set, Cudd_Not(V));
        Cudd_Ref(tmp);
        deref(pre_set);
        DdNode* tmp2 = Cudd_bddOr(manager, Cw, tmp);
        Cudd_Ref(tmp2);
        deref(Cw);
        deref(tmp);
        Cw = tmp2;
    }
    deref(one);

    DdNode** out = mxMalloc(sizeof(DdNode*)*mode);
    if (mode >= WIN_SET)
    {
        out[0] = V;
        if (mode >= WIN_CANDIDATE_SET)
            out[1] = Cw;
    }

    if (C_created)
    {
        deref(C[0]);
        mxFree(C);
    }
    return out;
}
