#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdint.h>
#include<string.h>
#include "mex.h"
#include "cudd.h"
#include "DynArray.h"
#include "system_types.h"

#define MAX(a,b) a > b ? a : b

BDDCont* make_simple_cont(Manager* manager, DdNode* set, DdNode* input, char* from)
{
    BDDCont cont;
    cont.mgr = manager;
    array_init(cont_sets, BDDlist, 1, 1);
    array_set(cont_sets, 0, set);
    cont.sets = cont_sets;
    array_make_persist(cont.sets);
    cont.subconts = NULL;
    cont.input = input;
    cont.type = SIMPLE;
    cont.mem_var = 0;
    int len = strlen(from);
    cont.from = mxMalloc(sizeof(char)*(len+1));
    mexMakeMemoryPersistent(cont.from);
    for (int i = 0; i < len; i++)
        cont.from[i] = from[i];
    cont.from[len] = '\0';
    mgr_incr(manager);

    BDDCont* ptr = mxMalloc(sizeof(BDDCont));
    *ptr = cont;
    mexMakeMemoryPersistent(ptr);
    return ptr;
}

BDDCont* make_reach_cont(Manager* manager, BDDlist* sets, BDDContList* subconts, char* from)
{
    if (array_len(sets) != array_len(subconts))
    {
        mexErrMsgIdAndTxt("mexBDD:ControllerCreation",
                          "reach controller: set_list and c_list must be same size");
    }
    BDDCont cont;
    cont.mgr = manager;
    cont.sets = sets;
    array_make_persist(cont.sets);
    cont.subconts = subconts;
    array_make_persist(cont.subconts);
    cont.input = NULL;
    cont.type = REACH;
    cont.mem_var = 0;
    int len = strlen(from);
    cont.from = mxMalloc(sizeof(char)*(len+1));
    mexMakeMemoryPersistent(cont.from);
    for (int i = 0; i < len; i++)
        cont.from[i] = from[i];
    cont.from[len] = '\0';
    mgr_incr(manager);

    BDDCont* ptr = mxMalloc(sizeof(BDDCont));
    *ptr = cont;
    mexMakeMemoryPersistent(ptr);
    return ptr;
}

BDDCont* make_recurrence_cont(Manager* manager, BDDlist* sets, BDDContList* subconts, char* from)
{
    if (array_len(sets) != array_len(subconts) + 1)
    {
        mexErrMsgIdAndTxt("mexBDD:ControllerCreation",
                          "recurrence controller: must have length(set_list) == length(c_list) + 1");
    }
    BDDCont cont;
    cont.mgr = manager;
    cont.sets = sets;
    array_make_persist(cont.sets);
    cont.subconts = subconts;
    array_make_persist(cont.subconts);
    cont.input = NULL;
    cont.type = RECURRENCE;
    cont.mem_var = 0;
    int len = strlen(from);
    cont.from = mxMalloc(sizeof(char)*(len+1));
    mexMakeMemoryPersistent(cont.from);
    for (int i = 0; i < len; i++)
        cont.from[i] = from[i];
    cont.from[len] = '\0';
    mgr_incr(manager);

    BDDCont* ptr = mxMalloc(sizeof(BDDCont));
    *ptr = cont;
    mexMakeMemoryPersistent(ptr);
    return ptr;
}

void cont_restrict(BDDCont* cont, DdNode* out_set)
{
    if (cont->type == SIMPLE)
    {
        Manager* mgr = cont->mgr;
        DdManager* manager = mgr->manager;
        DdNode* in_set = Cudd_bddSwapVariables(manager, out_set,
                array_list(mgr->s_out_vars), array_list(mgr->s_in_vars), mgr->s_var_num);
        Cudd_Ref(in_set);
        DdNode* cont_set = array_get(cont->sets, 0);
        DdNode* tmp = Cudd_bddAnd(manager, cont_set, out_set);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, cont_set);
        array_set(cont->sets, 0, tmp);
        tmp = Cudd_bddAnd(manager, cont->input, in_set);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, cont->input);
        Cudd_RecursiveDeref(manager, in_set);
        cont->input = tmp;
    }
    else
        mexErrMsgIdAndTxt("mexBDD:ControllerRestriction", "Complex controller cant be restricted");
}

void cont_set_from(BDDCont* cont, char* from)
{
    int len = strlen(from);
    cont->from = mxRealloc(cont->from, sizeof(char)*len);
    mexMakeMemoryPersistent(cont->from);
    for (int i = 0; i < len; i++)
        cont->from[i] = from[i];
}

DdNode* cont_get_input(BDDCont* cont, DdNode* states)
{
    Manager* mgr = cont->mgr;
    DdManager* ddmgr = mgr->manager;

    // states should be out-state BDD
    if (cont->type == SIMPLE)
    {
        if (!Cudd_EquivDC(ddmgr, array_get(cont->sets, 0), states, Cudd_Not(states)))
        {
            mexWarnMsgIdAndTxt("mexBDD:InputExtraction", "Set outside simple controller domain");
            return NULL;
        }
        return cont_get_simple_input(cont, states);
    }

    BDDCont* bottom_cont = cont;
    while (bottom_cont->type != SIMPLE)
    {
        if (bottom_cont->type == REACH)
        {
            int ind = MAX(0, bottom_cont->mem_var-1);
            if (Cudd_EquivDC(ddmgr, array_get(bottom_cont->sets, ind), states, Cudd_Not(states)))
            {
                /*
                Test for skipping the error of progress group sets (Zexiang)
                When creating the cont of pg, set V and a empty
                controller will always be added at first, which causes error
                when execute line 114 'a_list = bottom_cnt.subcontrollers(state);'
                */
                if (bottom_cont->mem_var != 0)
                    continue;
            }
            else if (!Cudd_EquivDC(ddmgr, array_get(bottom_cont->sets, cont->mem_var), states, Cudd_Not(states)))
            {
                // do whole search
                for (int i = 0; i < array_len(bottom_cont->sets); i++)
                {
                    if (Cudd_EquivDC(ddmgr, array_get(bottom_cont->sets, i), states, Cudd_Not(states)))
                    {
                        bottom_cont->mem_var = i;
                        break;
                    }
                }
            }
            if (!Cudd_EquivDC(ddmgr, array_get(bottom_cont->sets, cont->mem_var), states, Cudd_Not(states)))
            {
                mexWarnMsgIdAndTxt("mexBDD:InputExtraction", "Set outside controller domain");
                return NULL;
            }
        }
        else if (bottom_cont->type == RECURRENCE)
        {
            if (!Cudd_EquivDC(ddmgr, array_get(bottom_cont->sets, 0), states, Cudd_Not(states)))
            {
                mexWarnMsgIdAndTxt("mexBDD:InputExtraction", "Set outside controller domain");
                return NULL;
            }
            if (Cudd_EquivDC(ddmgr, array_get(bottom_cont->sets, bottom_cont->mem_var+1), states, Cudd_Not(states)))
            {
                if (bottom_cont->mem_var == array_len(bottom_cont->subconts)-1)
                {
                    bottom_cont->mem_var = 0;
                }
                else
                {
                    bottom_cont->mem_var++;
                }
            }
        }
        bottom_cont = array_get(bottom_cont->subconts, bottom_cont->mem_var);
    }
    return cont_get_simple_input(bottom_cont, states);
}

DdNode* cont_get_simple_input(BDDCont* cont, DdNode* out_states)
{
    Manager* mgr = cont->mgr;
    DdManager* ddmgr = mgr->manager;
    DdNode* in_states = Cudd_bddSwapVariables(ddmgr, out_states,
                        array_list(mgr->s_out_vars), array_list(mgr->s_in_vars), mgr->s_var_num);
    Cudd_Ref(in_states);
    DdNode* s_in_cube = Cudd_bddComputeCube(ddmgr, array_list(mgr->s_in_vars), NULL, mgr->s_var_num);
    Cudd_Ref(s_in_cube);
    DdNode* a_list = Cudd_bddAndAbstract(ddmgr, cont->input, in_states, s_in_cube);
    Cudd_Ref(a_list);
    Cudd_RecursiveDeref(ddmgr, in_states);
    Cudd_RecursiveDeref(ddmgr, s_in_cube);
    return a_list;
}

void free_manager(Manager* mgr)
{
    //mexPrintf("Clearing s_in_vars\n");
    array_free(mgr->s_in_vars);
    array_free(mgr->s_in_inds);
    //mexPrintf("Clearing s_out_vars\n");
    array_free(mgr->s_out_vars);
    array_free(mgr->s_out_inds);
    //mexPrintf("Clearing a_vars\n");
    array_free(mgr->a_vars);
    array_free(mgr->a_inds);
    //mexPrintf("Clearing s_encodings\n");
    //mexEvalString("drawnow;");
    for (int i = 0; i < array_len(mgr->s_encs); i++)
    {
       //mexPrintf("Clearing s_encoding %d\n", i);
       //mexEvalString("drawnow;");
       if (array_get(mgr->s_encs, i) != NULL)
            array_free(array_get(mgr->s_encs, i));
    }
    //mexPrintf("Clearing s_encoding list\n");
    //mexEvalString("drawnow;");
    array_free(mgr->s_encs);

    //mexPrintf("Clearing a_encodings\n");
    //mexEvalString("drawnow;");
    for (int i = 0; i < array_len(mgr->a_encs); i++)
    {
       //mexPrintf("Clearing a_encoding %d\n", i);
       //mexEvalString("drawnow;");
       if (array_get(mgr->a_encs, i) != NULL)
            array_free(array_get(mgr->a_encs, i));
    }
    //mexPrintf("Clearing a_encoding list\n");
    //mexEvalString("drawnow;");
    array_free(mgr->a_encs);
    //mexPrintf("Clearing enc_state_map\n");
    //mexEvalString("drawnow;");
    array_free(mgr->enc_state_map);
    //mexPrintf("Checking non-zero references: %d\n", Cudd_CheckZeroRef(mgr->manager));
    Cudd_Quit(mgr->manager);
}

void mgr_decr(Manager* mgr)
{
    mgr->struct_count--;
    //mexPrintf("Decreased struct count: %d\n", mgr->struct_count);
    if (mgr->struct_count <= 0)
        free_manager(mgr);
}

void mgr_incr(Manager* mgr)
{
    mgr->struct_count++;
    //mexPrintf("Increased struct count\n");
}
