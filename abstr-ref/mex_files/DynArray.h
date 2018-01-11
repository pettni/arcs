#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "mex.h"

#ifndef DYNARRAY
#define DYNARRAY

#define array_declare(DYN_ARRAY_TYPE, TYPE_REF)\
    struct DYN_ARRAY_TYPE\
    {\
        TYPE_REF* list;\
        size_t allocated_size;\
        size_t size;\
        size_t allocation_buffer;\
        int is_persistent;\
    };\
    typedef struct DYN_ARRAY_TYPE DYN_ARRAY_TYPE\

#define array_init(ARRAY, ARRAY_TYPE, SIZE, buffer)\
    ARRAY_TYPE* ARRAY = NULL;\
    if (buffer > 0) \
    {\
        ARRAY = mxMalloc(sizeof(ARRAY_TYPE));\
        ARRAY->list = mxMalloc(sizeof(array_type(ARRAY))*(SIZE+buffer)); \
        ARRAY->allocated_size = SIZE+buffer;\
        ARRAY->size = SIZE;\
        ARRAY->allocation_buffer = buffer;\
        ARRAY->is_persistent = 0;\
    }
#define array_free(ARRAY) {mxFree(ARRAY->list); mxFree(ARRAY);}

#define array_make_persist(ARRAY)\
{\
    mexMakeMemoryPersistent(ARRAY->list);\
    mexMakeMemoryPersistent(ARRAY);\
    ARRAY->is_persistent = 1;\
}

#define array_cpy(ARRAY, IND, coll, size)\
{\
    memcpy(&(array_get(ARRAY, IND)), coll, size);\
}

#define array_pushback(ARRAY, ELEMENT) {\
    if (ARRAY->size+1 >= ARRAY->allocated_size)\
        array_grow(ARRAY);\
    ARRAY->list[ARRAY->size] = ELEMENT;\
    ARRAY->size++;\
}

#define array_grow(ARRAY)\
{\
    ARRAY->allocated_size += ARRAY->allocation_buffer;\
    ARRAY->list = mxRealloc(ARRAY->list, sizeof(*(ARRAY->list))*ARRAY->allocated_size);\
    if (ARRAY->is_persistent)\
        mexMakeMemoryPersistent(ARRAY->list);\
}

#define array_alloc_size(ARRAY) ARRAY->allocated_size

#define array_growby(ARRAY, LEN)\
{\
    ARRAY->allocated_size += LEN;\
    ARRAY->list = mxRealloc(ARRAY->list, sizeof(*(ARRAY->list))*ARRAY->allocated_size);\
    if (ARRAY->is_persistent)\
        mexMakeMemoryPersistent(ARRAY->list);\
}

#define array_set(ARRAY, IND, ELEMENT)\
{\
    if (IND < ARRAY->size)\
        ARRAY->list[IND] = ELEMENT;\
}
#define array_get(ARRAY, IND) ARRAY->list[IND]

#define array_insert(ARRAY, IND, ELEMENT)\
{\
    uint arr_len = array_len(ARRAY);\
    if (IND < arr_len)\
    {\
        if (arr_len + 1 >= ARRAY->allocated_size)\
            array_grow(ARRAY);\
        ARRAY->size++;\
        memcpy(&(ARRAY->list[IND+1]), &(ARRAY->list[IND]), array_len(ARRAY)-1-IND);\
        ARRAY->list[IND] = ELEMENT;\
    }\
    else if (IND == arr_len)\
        array_pushback(ARRAY, ELEMENT);\
}

#define array_remove(ARRAY, IND)\
{\
    for (int i = IND+1; i < array_len(ARRAY); i++)\
        ARRAY->list[i-1] = ARRAY->list[i];\
    ARRAY->size--;\
}

#define array_len(ARRAY) ARRAY->size
#define array_type(ARRAY) __typeof__(*ARRAY->list)
#define array_list(ARRAY) ARRAY->list
#endif
