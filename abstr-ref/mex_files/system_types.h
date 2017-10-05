#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdint.h>
#include "cudd.h"
#include "DynArray.h"

#ifndef SYSTEM_TYPES
#define SYSTEM_TYPES

typedef uint32_t uint;
typedef struct Manager Manager;
typedef struct BDDSystem BDDSys;
typedef struct BDDController BDDCont;
typedef struct BDDlist BDDlist;
typedef struct NumList NumList;
typedef struct EncList EncList;
typedef struct BDDSysList BDDSysList;
typedef struct BDDContList BDDContList;

array_declare(BDDlist, DdNode*);
array_declare(NumList, uint);
array_declare(EncList, NumList*);
array_declare(BDDSysList, BDDSys*);
array_declare(BDDContList, BDDCont*);

typedef enum {SIMPLE, REACH, RECURRENCE} cont_type;

#define get_mgr(structure) structure->mgr->manager

BDDCont* make_simple_cont(Manager* manager, DdNode* set, DdNode* input, char* from);
BDDCont* make_reach_cont(Manager* manager, BDDlist* sets, BDDContList* subconts, char* from);
BDDCont* make_recurrence_cont(Manager* manager, BDDlist* sets, BDDContList* subconts, char* from);
void cont_restrict(BDDCont* cont, DdNode* set);
void cont_set_from(BDDCont* cont, char* from);
DdNode* cont_get_input(BDDCont* cont, DdNode* states);
DdNode* cont_get_simple_input(BDDCont* cont, DdNode* states);
void free_manager(Manager* mgr);
void mgr_decr(Manager* mgr);
void mgr_incr(Manager* mgr);

struct Manager {
     DdManager* manager;
     uint struct_count;
     BDDlist*    s_in_vars; // free
     NumList*    s_in_inds;
     DdNode*     s_in_cube;
     BDDlist*    s_out_vars; // free
     NumList*    s_out_inds;
     DdNode*     s_out_cube;
     BDDlist*    a_vars; // free
     NumList*    a_inds;
     DdNode*     a_cube;
     uint        s_var_num;
     uint        a_var_num;
     EncList*    s_encs; // free one level
     EncList*    a_encs; // free one level
     NumList*    enc_state_map; // free
};

struct BDDSystem {
     Manager*    mgr;
     DdNode*     trans_sys;
     DdNode*     all_states;
     DdNode*     all_actions;
     BDDlist*    pg_U; // free one level
     BDDlist*    pg_G; // free one level
     uint        pg_num;
};

struct BDDController
{
     Manager*     mgr;
     uint          ID;

     BDDlist*     sets; // free
     BDDContList* subconts; // free
     DdNode*      input;
     cont_type    type;
     int          mem_var;
     char*        from;
};

#endif
