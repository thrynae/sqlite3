/* -------------------------------------------------------------------------- */
/*                                                                            */
/* Version: 3.0.0                                                             */
/* Date:    2022-01-18                                                        */
/* Author:  H.J. Wisselink                                                    */
/* Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 )     */
/* Email = 'h_j_wisselink*alumnus_utwente_nl';                                */
/* Real_email = regexprep(Email,{'*','_'},{'@','.'})                          */
/*                                                                            */
/* -------------------------------------------------------------------------- */
/*                                                                            */
/* This file only required minor edits to make it conform to stricter         */
/* compiler standards.                                                        */
/*                                                                            */
/* -------------------------------------------------------------------------- */
/*                                                                            */
/* Original file (MIT license):                                               */
/* http://web.archive.org/web/202108id_/https://raw.githubusercontent.com/    */
/* rmartinjak/mex-sqlite3/master/structlist.c                                 */
/*                                                                            */
/* -------------------------------------------------------------------------- */

#include "mex.h"
#include "structlist_strict.h"

static struct node *nodecreate(mxArray *a) {
    struct node *n = mxMalloc(sizeof *n);
    n->array = mxDuplicateArray(a);
    n->next = NULL;
    return n;
}

static void nodefree(struct node *n) {
    if (!n) {
        return;
    }
    nodefree(n->next);
    mxFree(n);
}

void structlist_init(struct structlist *l)
{
    l->head = l->tail = NULL;
    l->num_fields = 0;
    l->size = 0;
}

void structlist_free(struct structlist *l)
{
    if (!l) {
        return;
    }
    nodefree(l->head);
    l->size = 0;
}


int structlist_add(struct structlist *l, mxArray *a)
{
    struct node *n = nodecreate(a);
    if (!n) {
        return -1;
    }
    if (!l->tail) {
        l->head = l->tail = n;
        l->num_fields = mxGetNumberOfFields(a);
    } else {
        if (mxGetNumberOfFields(a) != l->num_fields) {
            nodefree(n);
            return -2;
        }
        l->tail->next = n;
        l->tail = n;
    }
    l->size += mxGetN(a);
    return 0;
}

mxArray *structlist_collapse(struct structlist *l)
{
    int i, n;
    const char *name;
    size_t index;
    struct node *node;
    int field;
    mxArray *val;
    
    mxArray *a = mxCreateStructMatrix(1, l->size, 0, NULL);
    if (!l->size) {
        return a;
    }
    n = mxGetNumberOfFields(l->head->array);
    for (i = 0; i < n; i++) {
        name = mxGetFieldNameByNumber(l->head->array, i);
        mxAddField(a, name);
    }
    
    index = 0;
    for (node = l->head; node; node = node->next) {
        n = mxGetN(node->array);
        for (i = 0; i < n; i++) {
            for (field = 0; field < l->num_fields; field++) {
                val = mxGetFieldByNumber(node->array, i, field);
                mxSetFieldByNumber(a, index, field, mxDuplicateArray(val));
            }
            index++;
        }
        mxDestroyArray(node->array);
    }
    structlist_free(l);
    return a;
}

