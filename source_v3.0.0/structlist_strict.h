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
/* This file did not require any edits to make it conform to stricter         */
/* compiler standards. Only this header was added.                            */
/*                                                                            */
/* -------------------------------------------------------------------------- */
/*                                                                            */
/* Original file (MIT license):                                               */
/* http://web.archive.org/web/202108id_/https://raw.githubusercontent.com/    */
/* rmartinjak/mex-sqlite3/master/structlist.h                                 */
/*                                                                            */
/* -------------------------------------------------------------------------- */

#ifndef STRUCTLIST_H
#define STRUCTLIST_H

#include "mex.h"

struct structlist
{
    struct node {
        mxArray *array;
        struct node *next;
    } *head, *tail;
    int size;
    int num_fields;
};

void structlist_init(struct structlist *l);
void structlist_free(struct structlist *l);
int structlist_add(struct structlist *l, mxArray *a);
mxArray *structlist_collapse(struct structlist *l);
#endif

