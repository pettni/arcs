#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdint.h>
#include "cudd.h"
#include "DynArray.h"

#ifndef SYSTEM_TYPES
#define SYSTEM_TYPES

typedef uint32_t uint;
typedef struct BDDSystem BDDSys;
typedef struct BDDlist BDDlist;
typedef struct NumList NumList;
typedef struct EncList EncList;
typedef struct BDDSysList BDDSysList;

array_declare(BDDlist, DdNode*);
array_declare(NumList, uint);
array_declare(EncList, NumList*);
array_declare(BDDSysList, BDDSys*);

struct BDDSystem {
     DdManager*  manager;
     DdNode*     trans_sys;
     BDDlist*    s_in_vars; // free
     NumList*    s_in_inds;
     BDDlist*    s_out_vars; // free
     NumList*    s_out_inds;
     BDDlist*    a_vars; // free
     NumList*    a_inds;
     uint        s_var_num;
     uint        a_var_num;
     DdNode*     all_states;
     DdNode*     all_actions;
     EncList*    s_encs; // free one level
     EncList*    a_encs; // free one level
     NumList*    enc_state_map; // free
     BDDlist*    pg_U; // free one level
     BDDlist*    pg_G; // free one level
     uint        pg_num;
};

#endif
