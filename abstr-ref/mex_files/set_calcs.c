#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include "mex.h"
#include "cudd.h"
#include "system_types.h"
#include "set_calcs.h"

#define deref(system) Cudd_RecursiveDeref(get_mgr(sys), system);

#define WIN_SET 1
#define WIN_CANDIDATE_SET 2
#define WIN_CANDIDATE_CONT 3


OutS pre(BDDSys* sys, DdNode* X_prime, DdNode* A, char quant1, char quant2, int mode)
{
    DdManager* manager = get_mgr(sys);
    DdNode* T = sys->trans_sys;
    DdNode* valid_transitions;
    DdNode* tmp;
    DdNode* cont_set;

    DdNode* outCube = sys->mgr->s_out_cube;
    DdNode* aCube = sys->mgr->a_cube;

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
        if (mode >= WIN_CANDIDATE_CONT)
        {
            cont_set = Cudd_bddAnd(manager, valid_transitions, tmp);
            Cudd_Ref(cont_set);
        }
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
        if (mode >= WIN_CANDIDATE_CONT)
        {
            cont_set = Cudd_bddAnd(manager, valid_transitions, A);
            Cudd_Ref(cont_set);
            DdNode* tmp2 = Cudd_bddAnd(manager, cont_set, tmp);
            Cudd_Ref(tmp2);
            deref(cont_set);
            cont_set = tmp2;
        }
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
        if (mode >= WIN_CANDIDATE_CONT)
        {
            cont_set = Cudd_bddAnd(manager, valid_transitions, tmp);
            Cudd_Ref(cont_set);
        }
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
        if (mode >= WIN_CANDIDATE_CONT)
        {
            cont_set = Cudd_bddAnd(manager, valid_transitions, tmp);
            Cudd_Ref(cont_set);
        }
        Cudd_RecursiveDeref(manager, valid_transitions);
        valid_transitions = tmp;
    }
    else
    {
        char message[200];
        sprintf(message, "Invalid quantifiers! %c %c\n", quant1, quant2);
        mexErrMsgIdAndTxt("mexBDD:PreSetCalc", message);
    }

    OutS out;
    out.mode = mode;
    if (mode >= WIN_SET)
    {
        out.win_set = valid_transitions;
        if (mode >= WIN_CANDIDATE_CONT)
        {
            //mexPrintf("Controller created\n");
            //mexEvalString("drawnow;");
            Manager* mgr = sys->mgr;
            DdNode* transition_out = Cudd_bddSwapVariables(manager, valid_transitions,
                array_list(mgr->s_in_vars), array_list(mgr->s_out_vars), mgr->s_var_num);
            Cudd_Ref(transition_out);
            out.cont = make_simple_cont(mgr, transition_out, cont_set, "pre");
        }
    }

    return out;
}

OutS inv(BDDSys* sys, DdNode* Z, DdNode* B, DdNode* U, DdNode* G,
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
    DdNode* cont_set = NULL;

    OutS in = pre(sys, Z, U, quant, 'e', WIN_SET);
    DdNode* quick_test = in.win_set;
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
            if (mode >= WIN_CANDIDATE_CONT)
            {
                cont_set = Cudd_ReadLogicZero(manager);
                Cudd_Ref(cont_set);
            }
        }
        OutS out;
        out.mode = mode;
        if (mode >= WIN_SET)
        {
            out.win_set = Y;
            if (mode >= WIN_CANDIDATE_SET)
            {
                out.cand_set = Cw;
                if (mode >= WIN_CANDIDATE_CONT)
                {
                    DdNode* extra_zero = Cudd_ReadLogicZero(manager);
                    Cudd_Ref(extra_zero);
                    out.cont = make_simple_cont(sys->mgr, cont_set, extra_zero, "pginv");
                }
            }
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
        OutS pre_set_in = pre(sys, tmp, U, quant, 'a', WIN_SET);
        DdNode* pre_set = pre_set_in.win_set;
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

    OutS out;
    out.mode = mode;
    if (mode >= WIN_SET)
    {
        out.win_set = Y;
        if (mode >= WIN_CANDIDATE_SET)
        {
            out.cand_set = Cw;
            if (mode >= WIN_CANDIDATE_CONT)
            {
                //mexPrintf("Controller created\n");
                //mexEvalString("drawnow;");
                tmp = Cudd_bddOr(manager, Y, Z);
                Cudd_Ref(tmp);
                in = pre(sys, tmp, U, quant, 'a', mode);
                out.cont = in.cont;
                cont_restrict(out.cont, Y);
                cont_set_from(out.cont, "pginv");
                deref(in.win_set);
                deref(tmp);
            }
        }
    }

    return out;
}

OutS PGpre(BDDSys* sys, DdNode* Z, DdNode* B, char quant, int mode)
{
    DdManager* manager = get_mgr(sys);
    DdNode* U;
    DdNode* G;
    DdNode* PG_set = Cudd_bddAnd(manager, Z, Cudd_ReadOne(manager));
    Cudd_Ref(PG_set);
    DdNode* inv_set;
    DdNode* tmp;
    DdNode* Cw;
    BDDContList* cont_list;
    BDDlist* set_list;

    if (mode >= WIN_CANDIDATE_SET)
    {
        Cw = Cudd_ReadLogicZero(manager);
        Cudd_Ref(Cw);
        if (mode >= WIN_CANDIDATE_CONT)
        {
            array_init(list, BDDContList, 0, 10);
            array_init(list2, BDDlist, 0, 10);
            cont_list = list;
            set_list = list2;
            DdNode* tmp_zero = Cudd_ReadLogicZero(manager);
            Cudd_Ref(tmp_zero);
            DdNode* Z_copy_1 = Cudd_bddAnd(manager, Z, Cudd_ReadOne(manager));
            Cudd_Ref(Z_copy_1);
            BDDCont* first_cont = make_simple_cont(sys->mgr, Z_copy_1, tmp_zero, "");
            array_pushback(cont_list, first_cont);
            DdNode* Z_copy_2 = Cudd_bddAnd(manager, Z, Cudd_ReadOne(manager));
            Cudd_Ref(Z_copy_2);
            array_pushback(set_list, Z_copy_2);
        }
    }

    for (int i = 0; i < array_len(sys->pg_U); i++)
    {

        if (array_get(sys->pg_U, i) == NULL)
            continue;
        if (quant == 'a' && !Cudd_EquivDC(manager, array_get(sys->pg_U, i), sys->all_actions, Cudd_Not(sys->all_actions)))
        {
            //mexPrintf("This happens! Progress group skipped\n");
            continue;
        }

        U = array_get(sys->pg_U, i);
        G = array_get(sys->pg_G, i);

        OutS in = inv(sys, PG_set, B, U, G, quant, mode);

        if (mode >= WIN_SET)
        {
            inv_set = in.win_set;
            if (mode >= WIN_CANDIDATE_SET)
            {
                deref(Cw);
                Cw = in.cand_set;
                if (mode >= WIN_CANDIDATE_CONT)
                {
                    DdNode* inv_set_copy = Cudd_bddAnd(manager, inv_set, Cudd_ReadOne(manager));
                    Cudd_Ref(inv_set_copy);
                    array_pushback(set_list, inv_set_copy);
                    array_pushback(cont_list, in.cont);
                }
            }
        }

        tmp = Cudd_bddOr(manager, inv_set, PG_set);
        Cudd_Ref(tmp);
        deref(inv_set);
        deref(PG_set);
        PG_set = tmp;
    }

    OutS out;
    out.mode = mode;
    if (mode >= WIN_SET)
    {
        out.win_set = PG_set;
        if (mode >= WIN_CANDIDATE_SET)
        {
            out.cand_set = Cw;
            if (mode >= WIN_CANDIDATE_CONT)
            {
                out.cont = make_reach_cont(sys->mgr, set_list, cont_list, "pre_pg");
            }
        }
    }

    return out;
}

OutS win_until(BDDSys* sys, DdNode* Z, DdNode* B, char quant, int mode)
{
    DdManager* manager = get_mgr(sys);
    // calculates Win_q,a(B U Z)
    DdNode* PGpre_set;
    DdNode* pre_set;
    DdNode* tmp;
    DdNode* Cw;
    BDDContList* cont_list;
    BDDlist* set_list;
    if (mode >= WIN_CANDIDATE_SET)
    {
        Cw = Cudd_ReadLogicZero(manager);
        Cudd_Ref(Cw);
        if (mode >= WIN_CANDIDATE_CONT)
        {
            array_init(list, BDDContList, 0, 10);
            array_init(list2, BDDlist, 0, 10);
            cont_list = list;
            set_list = list2;
        }
    }
    DdNode* one = Cudd_ReadOne(manager);
    Cudd_Ref(one);
    DdNode* X = Cudd_ReadLogicZero(manager);
    Cudd_Ref(X);
    DdNode* old_X = Cudd_ReadOne(manager);
    Cudd_Ref(old_X);
    while(!Cudd_EquivDC(manager, X, old_X, Cudd_ReadLogicZero(manager)))
    {
        int dispose_PG = 1;
        deref(old_X);
        old_X = X;

        OutS pre_in = pre(sys, X, one, quant, 'a', mode);
        pre_set = pre_in.win_set;
        tmp = Cudd_bddSwapVariables(manager, pre_set, array_list(sys->mgr->s_in_vars), array_list(sys->mgr->s_out_vars), sys->mgr->s_var_num);
        Cudd_Ref(tmp);
        deref(pre_set);
        pre_set = tmp;
        tmp = Cudd_bddAnd(manager, B, pre_set);
        Cudd_Ref(tmp);
        X = tmp;
        tmp = Cudd_bddOr(manager, Z, X);
        Cudd_Ref(tmp);
        deref(X);
        X = tmp;

        OutS in;
        if (mode == WIN_SET)
        {
            in = PGpre(sys, old_X, B, quant, mode);
            PGpre_set = in.win_set;
        }
        else if (mode >= WIN_CANDIDATE_SET)
        {
            in = PGpre(sys, X, B, quant, mode);
            PGpre_set = in.win_set;
            deref(Cw);
            Cw = in.cand_set;
            if (mode >= WIN_CANDIDATE_CONT)
            {
                array_pushback(set_list, pre_set)
                array_pushback(cont_list, pre_in.cont);

                tmp = Cudd_bddAnd(manager, PGpre_set, Cudd_Not(X));
                Cudd_Ref(tmp);
                if (!Cudd_EquivDC(manager, tmp, Cudd_ReadLogicZero(manager), Cudd_ReadLogicZero(manager)))
                {
                    array_pushback(set_list, PGpre_set);
                    array_pushback(cont_list, in.cont);
                    dispose_PG = 0;
                }
            }
        }
        tmp = Cudd_bddOr(manager, PGpre_set, X);
        Cudd_Ref(tmp);
        deref(X);
        X = tmp;
        if (mode < WIN_CANDIDATE_CONT)
            deref(pre_set);
        if (dispose_PG)
            deref(PGpre_set);
    }

    if (mode >= WIN_CANDIDATE_SET)
    {
        OutS pre_in = pre(sys, X, one, quant, 'a', WIN_SET);
        DdNode* XPre = pre_in.win_set;
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

    OutS out;
    out.mode = mode;
    if (mode >= WIN_SET)
    {
        out.win_set = X;
        if (mode >= WIN_CANDIDATE_SET)
        {
            out.cand_set = Cw;
            if (mode >= WIN_CANDIDATE_CONT)
                out.cont = make_reach_cont(sys->mgr, set_list, cont_list, "win_until");
        }
    }
    return out;
}

OutS win_until_and_always(BDDSys* sys, DdNode* A, DdNode* B, DdNode* Z, char quant, int mode)
{
    DdManager* manager = get_mgr(sys);
    DdNode* always_A = Cudd_ReadOne(manager);
    Cudd_Ref(always_A);
    DdNode* win_until_set;
    DdNode* B_always;
    DdNode* Z_always;
    DdNode* Cw;
    BDDCont* cont;
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

        int interm_mode = mode < WIN_CANDIDATE_SET ? mode : WIN_CANDIDATE_SET;
        OutS in = win_intermediate(sys, sys->all_states, A, zero, &one, 1, quant, interm_mode);
        Cudd_RecursiveDeref(manager, zero);
        Cudd_RecursiveDeref(manager, one);
        if (mode >= WIN_SET)
        {
            always_A = in.win_set;
            if (mode >= WIN_CANDIDATE_SET)
            {
                deref(Cw);
                Cw = in.cand_set;
            }
        }
    }

    B_always = Cudd_bddAnd(manager, always_A, B);
    Cudd_Ref(B_always);
    Z_always = Cudd_bddAnd(manager, always_A, Z);
    Cudd_Ref(Z_always);
    deref(always_A);

    OutS in = win_until(sys, Z_always, B_always, quant, mode);
    Cudd_RecursiveDeref(manager, B_always);
    Cudd_RecursiveDeref(manager, Z_always);
    if (mode >= WIN_SET)
    {
        win_until_set = in.win_set;
        if (mode >= WIN_CANDIDATE_SET)
        {
            tmp = Cudd_bddOr(manager, Cw, in.cand_set);
            Cudd_Ref(tmp);
            deref(Cw);
            deref(in.cand_set);
            Cw = tmp;
            if (mode >= WIN_CANDIDATE_CONT)
            {
                cont = in.cont;
                cont_set_from(cont, "win_until_and_always");
            }
        }
    }

    OutS out;
    out.mode = mode;
    if (mode >= WIN_SET)
    {
        out.win_set = win_until_set;
        if (mode >= WIN_CANDIDATE_SET)
        {
            out.cand_set = Cw;
            if (mode >= WIN_CANDIDATE_CONT)
                out.cont = cont;
        }
    }
    return out;
}

OutS win_intermediate(BDDSys* sys, DdNode* A, DdNode* B, DdNode* Z,
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
    BDDContList* cont_list;
    if (mode >= WIN_CANDIDATE_SET)
    {
        Cw = Cudd_ReadLogicZero(manager);
        Cudd_Ref(Cw);
        if (mode >= WIN_CANDIDATE_CONT)
        {
            array_init(list, BDDContList, 0, 10);
            cont_list = list;
        }
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
        OutS pre_set_in = pre(sys, W, one, quant, 'a', WIN_SET);
        pre_set = pre_set_in.win_set;
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

            OutS in = win_until_and_always(sys, A, B, Z_iter[i], quant, mode);
            if (mode >= WIN_SET)
            {
                win_until_set = in.win_set;
                if (first_set && mode >= WIN_CANDIDATE_SET)
                {
                    tmp = Cudd_bddOr(manager, Cw, in.cand_set);
                    Cudd_Ref(tmp);
                    deref(Cw);
                    deref(in.cand_set);
                    Cw = tmp;
                }
                else if (mode >= WIN_CANDIDATE_SET)
                {
                    deref(in.cand_set);
                }
                if (mode >= WIN_CANDIDATE_CONT)
                {
                    if (i < array_len(cont_list))
                    {
                        array_set(cont_list, i, in.cont);
                    }
                    else
                        array_pushback(cont_list, in.cont);
                }
            }
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

    OutS out;
    out.mode = mode;
    if (mode >= WIN_SET)
    {
        out.win_set = W;
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
            out.cand_set = Cw;
            if (mode >= WIN_CANDIDATE_CONT)
            {
                array_init(set_list, BDDlist, 0, C_num+1);
                DdNode* W_copy = Cudd_bddAnd(manager, W, Cudd_ReadOne(manager));
                Cudd_Ref(W_copy);
                array_pushback(set_list, W_copy);
                for (int i = 0; i < C_num; i++)
                {
                    tmp = Cudd_bddAnd(manager, C[i], B);
                    Cudd_Ref(tmp);
                    array_pushback(set_list, tmp);
                }
                //mexPrintf("set_list = %d, c_list = %d\n", array_len(set_list), array_len(cont_list));
                out.cont = make_recurrence_cont(sys->mgr, set_list, cont_list, "win_intermediate");
            }
        }
    }

    return out;
}

OutS win_primal(BDDSys* sys, DdNode* A, DdNode* B, DdNode** C, int C_num,
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
    BDDContList* cont_list;
    BDDlist* set_list;
    if (mode >= WIN_CANDIDATE_SET)
    {
        Cw = Cudd_ReadLogicZero(manager);
        Cudd_Ref(Cw);
        if (mode >= WIN_CANDIDATE_CONT)
        {
            array_init(list, BDDContList, 0, 10);
            array_init(list2, BDDlist, 0, 10);
            cont_list = list;
            set_list = list2;
        }
    }

    if (Cudd_EquivDC(manager, A, Cudd_ReadLogicZero(manager), Cudd_ReadLogicZero(manager)))
         A = one;
    if (Cudd_EquivDC(manager, B, Cudd_ReadLogicZero(manager), Cudd_ReadLogicZero(manager)))
         B = one;
    int first_set = 1;
    int remove_old = 1;
    while (!Cudd_EquivDC(manager, V, old_V, Cudd_ReadLogicZero(manager)))
    {
        if (mode < WIN_CANDIDATE_CONT)
        {
            deref(old_V);
        }
        else
            remove_old = 0;
        old_V = V;
        OutS pre_set_in = pre(sys, V, one, quant1, 'a', WIN_SET);
        pre_set = pre_set_in.win_set;
        tmp = Cudd_bddSwapVariables(manager, pre_set, array_list(sys->mgr->s_in_vars), array_list(sys->mgr->s_out_vars), sys->mgr->s_var_num);
        Cudd_Ref(tmp);
        deref(pre_set);
        pre_set = tmp;
        OutS PGpre_in = PGpre(sys, V, A, quant1, WIN_SET);
        pg_pre_set = PGpre_in.win_set;
        Z = Cudd_bddOr(manager, pre_set, pg_pre_set);
        Cudd_Ref(Z);
        deref(pre_set);
        deref(pg_pre_set);

        OutS in = win_intermediate(sys, A, B, Z, C, C_num, quant1, mode);
        if (mode >= WIN_SET)
        {
            V = in.win_set;
            if (first_set && mode >= WIN_CANDIDATE_SET)
            {
                deref(Cw);
                Cw = in.cand_set;
                first_set = 0;
            }
            else if (mode >= WIN_CANDIDATE_SET)
            {
                deref(in.cand_set);
            }
            if (mode >= WIN_CANDIDATE_CONT)
            {
                array_pushback(set_list, V);
                array_pushback(cont_list, in.cont);
            }
        }
        deref(Z);
    }
    if (remove_old)
        deref(old_V);

    if (mode >= WIN_CANDIDATE_SET)
    {
        OutS in = pre(sys, V, one, quant1, 'e', WIN_SET);
        pre_set = in.win_set;
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

    OutS out;
    out.mode = mode;
    if (mode >= WIN_SET)
    {
        out.win_set = V;
        if (mode >= WIN_CANDIDATE_SET)
        {
            out.cand_set = Cw;
            if (mode >= WIN_CANDIDATE_CONT)
                out.cont = make_reach_cont(sys->mgr, set_list, cont_list, "win_primal");
        }
    }

    if (C_created)
    {
        deref(C[0]);
        mxFree(C);
    }
    return out;
}
