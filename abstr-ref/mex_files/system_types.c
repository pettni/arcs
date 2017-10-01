#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdint.h>
#include "mex.h"
#include "cudd.h"
#include "DynArray.h"
#include "system_types.h"

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
    cont.mem_var = 1;
    cont.from = from;
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
    cont.subconts = subconts;
    cont.input = NULL;
    cont.type = REACH;
    cont.mem_var = 1;
    cont.from = from;

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
    cont.subconts = subconts;
    cont.input = NULL;
    cont.type = RECURRENCE;
    cont.mem_var = 1;
    cont.from = from;

    BDDCont* ptr = mxMalloc(sizeof(BDDCont));
    *ptr = cont;
    mexMakeMemoryPersistent(ptr);
    return ptr;
}

void free_manager(Manager* mgr)
{
    mexPrintf("Clearing s_in_vars\n");
    array_free(mgr->s_in_vars);
    array_free(mgr->s_in_inds);
    mexPrintf("Clearing s_out_vars\n");
    array_free(mgr->s_out_vars);
    array_free(mgr->s_out_inds);
    mexPrintf("Clearing a_vars\n");
    array_free(mgr->a_vars);
    array_free(mgr->a_inds);
    mexPrintf("Clearing s_encodings\n");
    mexEvalString("drawnow;");
    for (int i = 0; i < array_len(mgr->s_encs); i++)
    {
       mexPrintf("Clearing s_encoding %d\n", i);
       mexEvalString("drawnow;");
       if (array_get(mgr->s_encs, i) != NULL)
            array_free(array_get(mgr->s_encs, i));
    }
    mexPrintf("Clearing s_encoding list\n");
    mexEvalString("drawnow;");
    array_free(mgr->s_encs);

    mexPrintf("Clearing a_encodings\n");
    mexEvalString("drawnow;");
    for (int i = 0; i < array_len(mgr->a_encs); i++)
    {
       mexPrintf("Clearing a_encoding %d\n", i);
       mexEvalString("drawnow;");
       if (array_get(mgr->a_encs, i) != NULL)
            array_free(array_get(mgr->a_encs, i));
    }
    mexPrintf("Clearing a_encoding list\n");
    mexEvalString("drawnow;");
    array_free(mgr->a_encs);
    mexPrintf("Clearing enc_state_map\n");
    mexEvalString("drawnow;");
    array_free(mgr->enc_state_map);
    Cudd_Quit(mgr->manager);
}

void mgr_decr(Manager* mgr)
{
    mgr->struct_count--;
    if (mgr->struct_count <= 0)
        free_manager(mgr);
}

void mgr_incr(Manager* mgr)
{
    mgr->struct_count++;
}
