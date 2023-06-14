/* -------------------------------------------------------------------------- */
/*                                                                            */
/* Version: 3.1.0                                                             */
/* Date:    2023-06-14                                                        */
/* Author:  H.J. Wisselink                                                    */
/* Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 )     */
/* Email = 'h_j_wisselink*alumnus_utwente_nl';                                */
/* Real_email = regexprep(Email,{'*','_'},{'@','.'})                          */
/*                                                                            */
/* -------------------------------------------------------------------------- */
/*                                                                            */
/* This file contains the main interface code between the SQLite3 C project   */
/* code and MATLAB/Octave.                                                    */
/* The original file by Robin Martinjak has been edited to let it conform to  */
/* the stricter standards of older compilers. It should now adhere to the     */
/* ANSI C standard and be compatible with all supported MATLAB compilers      */
/* (except the 32-bit LCC, which are missing required integer support).       */
/* While actual additional error handling has not been added, the granularity */
/* of the reported error has been improved by including the description from  */
/* the SQLite3 project file.                                                  */
/*                                                                            */
/* This mex file has also been adapted to return the version number if called */
/* with no inputs and one output.                                             */
/*                                                                            */
/* Naming rules for columns in SQLite are much more lenient than fielnames in */
/* MATLAB. To allow external handling of invalid fieldnames, the actual names */
/* are returned as a second output argument in a cell vector.                 */
/*                                                                            */
/* -------------------------------------------------------------------------- */
/*                                                                            */
/* Original file (MIT license):                                               */
/* http://web.archive.org/web/202108id_/https://raw.githubusercontent.com/    */
/* rmartinjak/mex-sqlite3/master/sqlite3.c                                    */
/*                                                                            */
/* -------------------------------------------------------------------------- */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "sqlite3.h"

#include "mex.h"

#include "structlist_strict.h"
#include "CharConversions.c"

mxClassID dummy;

/* Yes, this is terrible, but it works. */
int TYPECHECK(const mxArray *ARRAY, int nargin,
               const mxClassID arg0, const mxClassID arg1, const mxClassID arg2,
               const mxClassID arg3, const mxClassID arg4, const mxClassID arg5) {
    int arg,ok=0;
    mxClassID type_got = mxGetClassID(ARRAY);
    mxClassID type_expect[6];
    
    if (nargin>0) { type_expect[0]=arg0; }
    if (nargin>1) { type_expect[1]=arg1; }
    if (nargin>2) { type_expect[2]=arg2; }
    if (nargin>3) { type_expect[3]=arg3; }
    if (nargin>4) { type_expect[4]=arg4; }
    if (nargin>5) { type_expect[5]=arg5; }
    
    for (arg=0;arg<nargin;arg++)
    {
        if (type_got == type_expect[arg]) {
            ok = 1;
        }
    }
    if (!ok) {
        /* Unsupported data type detected. */
        return SQLITE_MISMATCH;
    }
    return SQLITE_OK;
}


int bind_string(sqlite3_stmt *stmt, int index, const mxArray *array)
{
    int res;
    char *s;
    res = TYPECHECK(array,1, mxCHAR_CLASS ,dummy,dummy,dummy,dummy,dummy);
    if (res != SQLITE_OK) {
        return res;
    }
    s = WRAPPERmxArrayToUTF8String(array);
    res = sqlite3_bind_text(stmt, index, s,
            -1,  /* ALL the bytes */
            SQLITE_TRANSIENT  /* make a copy */
            );
    free(s);
    return res;
}

int bind_double(sqlite3_stmt *stmt, int index, const mxArray *array)
{
    int res;
    double val;
    res = TYPECHECK(array,2, mxSINGLE_CLASS, mxDOUBLE_CLASS ,dummy,dummy,dummy,dummy);
    if (res != SQLITE_OK) {
        return res;
    }
    val = mxGetScalar(array);
    return sqlite3_bind_double(stmt, index, val);
}

int bind_int64(sqlite3_stmt *stmt, int index, const mxArray *array)
{
    int res;
    mxClassID cls;
    int64_t val = 0;
    res = TYPECHECK(array,6,
              mxINT16_CLASS, mxUINT16_CLASS,
              mxINT32_CLASS, mxUINT32_CLASS,
              mxINT64_CLASS, mxUINT64_CLASS);
    if (res != SQLITE_DONE) {
        return res;
    }

    cls = mxGetClassID(array);
    #define GET_VALUE(CLASS, TYPE) \
    if (cls == CLASS) { \
    TYPE *p = mxGetData(array); \
    val = *p; \
    }
    GET_VALUE(mxINT16_CLASS, int16_t);
    GET_VALUE(mxUINT16_CLASS, uint16_t);
    GET_VALUE(mxINT16_CLASS, int32_t);
    GET_VALUE(mxUINT16_CLASS, uint32_t);
    GET_VALUE(mxINT16_CLASS, int64_t);
    GET_VALUE(mxUINT16_CLASS, uint64_t);
    return sqlite3_bind_int64(stmt, index, val);
}

int bind_params(sqlite3_stmt* stmt, const mxArray *params, int column)
{
    int res, i, n = mxGetNumberOfFields(params);
    res = TYPECHECK(params,1, mxSTRUCT_CLASS ,dummy,dummy,dummy,dummy,dummy);
    if (res != SQLITE_OK) {
        return res;
    }
    for (i = 0; i < n; i++) {
        mxArray *array = mxGetFieldByNumber(params, column, i);
        mxClassID cls = mxGetClassID(array);
        int res;
        switch (cls) {
            case mxFUNCTION_CLASS:
                break;
                
            case mxCHAR_CLASS:
                res = bind_string(stmt, i + 1, array);
                break;
                
            case mxSINGLE_CLASS:
            case mxDOUBLE_CLASS:
                res = bind_double(stmt, i + 1, array);
                break;
                
            default:  /* anything else is an integer */
                res = bind_int64(stmt, i + 1, array);
        }
        if (res != SQLITE_OK) {
            return res;
        }
    }
    return SQLITE_OK;
}

int execute_many(sqlite3 *db, sqlite3_stmt *stmt, const mxArray *params)
{
    int res, i, n = mxGetN(params);
    for (i = 0; i < n; i++) {
        res = bind_params(stmt, params, i);
        if (res != SQLITE_OK) {
            return res;
        }

        res = sqlite3_step(stmt);
        if (res != SQLITE_DONE) {
            return res;
        }
        sqlite3_reset(stmt);
        sqlite3_clear_bindings(stmt);
    }
    return SQLITE_DONE;
}

mxArray *get_column(sqlite3_stmt *stmt, int index)
{
    int type = sqlite3_column_type(stmt, index);
    if (type == SQLITE_TEXT) {
        const char *text = sqlite3_column_text(stmt, index);
        return WRAPPERmxCreateCharMatrixFromStrings(1, &text);
    } else if (type == SQLITE_FLOAT) {
        return mxCreateDoubleScalar(sqlite3_column_double(stmt, index));
    } else if (type == SQLITE_INTEGER) {
        mxArray *a = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
        int64_t val = sqlite3_column_int64(stmt, index);
        memcpy(mxGetData(a), &val, sizeof val);
        return a;
    } else if (type == SQLITE_NULL) {
        /* This return a NULL as a NaN.*/
        return mxCreateDoubleScalar(mxGetNaN());
    } else if (type == SQLITE_BLOB) {
        /* This returns an error later. */
        return NULL;
    }
    /* This should never happen. */
    return NULL;
}

int fetch_row(sqlite3_stmt *stmt, int num_columns, mxArray **array)
{
    mxArray *data;
    int i, res = sqlite3_step(stmt);
    if (res == SQLITE_ROW) {
        for (i = 0; i < num_columns; i++) {
            data = get_column(stmt, i);
            if (data)
            {
                /* Write data to the mxArray. */
                mxSetFieldByNumber(*array, 0, i, data);
            }
            else
            {
                /* Throw "Unsupported column type" error. */
                return SQLITE_FORMAT;
            }
        }
    }
    return res;
}


/* This function is not used, but could be used in fetch_results() to detect */
/* invalid field names and generate a MATLAB-compatible field name. */
int validate_field_name(const char *name)
{
    /* Check whether the field name is actually allowed in Matlab. */
    int isValid,i,currentChar;
    currentChar = name[0];
    currentChar = (int)currentChar;
    if (!((currentChar>=48 && currentChar<=57) ||
          (currentChar>=65 && currentChar<=90) ||
          (currentChar>=97 && currentChar<=122)
         )) {
        return isValid=0;
    }
    for (i = 1; i < strlen(name); i++) {
        currentChar = name[i];
        currentChar = (int)currentChar;
        if (!((currentChar>=48 && currentChar<=57) ||
              (currentChar>=65 && currentChar<=90) ||
              (currentChar>=97 && currentChar<=122) ||
              currentChar==95)) {
            return isValid=0;
        }
    }
    return isValid=1;
}

int fetch_results(sqlite3_stmt *stmt, mxArray **results, mxArray **columnnames)
{
    /* Fetch the data from the data and store the results in the first output */
    /* argument. The second output argument is used to keep track of the      */
    /* SQLite column names in case they are invalid.                          */
    mxArray *row;mxArray *nameAsmxChar;int i, num_columns, res;
    struct structlist result_list;
    structlist_init(&result_list);
    num_columns = sqlite3_column_count(stmt);
    /* Create a struct element with the correct fields. */
    row = mxCreateStructMatrix(1, 1, 0, NULL);
    /* Create a cell vector to hold the SQLite column names. */
    *columnnames = mxCreateCellMatrix(num_columns, 1);
    for (i = 0; i < num_columns; i++) {
        const char *name = sqlite3_column_name(stmt, i);
        nameAsmxChar = WRAPPERmxCreateCharMatrixFromStrings(1, &name);
        mxSetCell(*columnnames, i, nameAsmxChar);
        sprintf((char *)name,"field_%d",i+1);
        mxAddField(row, (const char *)name);
    }
    
    /* Fetch the data and store it in a new row in the struct. */
    while (res = fetch_row(stmt, num_columns, &row),
            res == SQLITE_ROW) {
        structlist_add(&result_list, row);
    }
    *results = structlist_collapse(&result_list);
    
    return res;
}

void
        mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char *db_path, *query;
    sqlite3 *db;
    sqlite3_stmt *stmt;
    int res,num_columns;
    
    enum { RHS_FILENAME, RHS_QUERY, RHS_PARAMS };
    if (nlhs==1 && nrhs==0) {
        plhs[0] = mxCreateString("3.1.0");
        return;
    } else {
        (void) nlhs;
    }
    
    if (nrhs < 2) {
        mexErrMsgIdAndTxt("sqlite3:sqlite3",
                "usage: sqlite3(file, query[, params])");
    }
    
    db_path = WRAPPERmxArrayToUTF8String(prhs[RHS_FILENAME]);
    query = WRAPPERmxArrayToUTF8String(prhs[RHS_QUERY]);
    
    
    if (sqlite3_open(db_path, &db) != SQLITE_OK) {
        mexErrMsgIdAndTxt("sqlite3:open", "failed to open db");
    }
    
    res = sqlite3_prepare_v2(db, query, -1, &stmt, NULL);
    if (res != SQLITE_OK) {
        sqlite3_close(db);
        mexErrMsgIdAndTxt(
            "sqlite3:prepare",
            "Failed to prepare query.\n"
            "Reasons could be invalid syntax or creating a table that already exists.\n"
            "sqlite3_prepare_v2 returned status code %d",
            res);
    }
    
    num_columns = sqlite3_column_count(stmt);
    
    if (num_columns > 0) {
        /* Looks like a SELECT query */
        res = fetch_results(stmt, &plhs[0], &plhs[1]);
    } else {
        /* If the user passes a struct array, use this to execute the statement
         * for each column */
        if (nrhs > 2) {
            res = execute_many(db, stmt, prhs[RHS_PARAMS]);
        } else {
            res = sqlite3_step(stmt);
        }
    }
    
    if (res == SQLITE_DONE) {
        /* success! */
        sqlite3_finalize(stmt);
        sqlite3_close(db);
        return;
    }

    /* Close the database file before throwing an error. */
    sqlite3_close(db);

    if(res == SQLITE_ERROR)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_ERROR",
                "SQLite responded with error code %d:\n%s",
                 res,"Generic error");
    } else if(res == SQLITE_INTERNAL)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_INTERNAL",
                "SQLite responded with error code %d:\n%s",
                 res,"Internal logic error in SQLite");
    } else if(res == SQLITE_PERM)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_PERM",
                "SQLite responded with error code %d:\n%s",
                 res,"Access permission denied");
    } else if(res == SQLITE_ABORT)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_ABORT",
                "SQLite responded with error code %d:\n%s",
                 res,"Callback routine requested an abort");
    } else if(res == SQLITE_BUSY)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_BUSY",
                "SQLite responded with error code %d:\n%s",
                 res,"The database file is locked");
    } else if(res == SQLITE_LOCKED)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_LOCKED",
                "SQLite responded with error code %d:\n%s",
                 res,"A table in the database is locked");
    } else if(res == SQLITE_NOMEM)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_NOMEM",
                "SQLite responded with error code %d:\n%s",
                 res,"A malloc() failed");
    } else if(res == SQLITE_READONLY)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_READONLY",
                "SQLite responded with error code %d:\n%s",
                 res,"Attempt to write a readonly database");
    } else if(res == SQLITE_INTERRUPT)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_INTERRUPT",
                "SQLite responded with error code %d:\n%s",
                 res,"Operation terminated by sqlite3_interrupt()");
    } else if(res == SQLITE_IOERR)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_IOERR",
                "SQLite responded with error code %d:\n%s",
                 res,"Some kind of disk I/O error occurred");
    } else if(res == SQLITE_CORRUPT)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_CORRUPT",
                "SQLite responded with error code %d:\n%s",
                 res,"The database disk image is malformed");
    } else if(res == SQLITE_NOTFOUND)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_NOTFOUND",
                "SQLite responded with error code %d:\n%s",
                 res,"Unknown opcode in sqlite3_file_control()");
    } else if(res == SQLITE_FULL)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_FULL",
                "SQLite responded with error code %d:\n%s",
                 res,"Insertion failed because database is full");
    } else if(res == SQLITE_CANTOPEN)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_CANTOPEN",
                "SQLite responded with error code %d:\n%s",
                 res,"Unable to open the database file");
    } else if(res == SQLITE_PROTOCOL)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_PROTOCOL",
                "SQLite responded with error code %d:\n%s",
                 res,"Database lock protocol error");
    } else if(res == SQLITE_EMPTY)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_EMPTY",
                "SQLite responded with error code %d:\n%s",
                 res,"Internal use only");
    } else if(res == SQLITE_SCHEMA)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_SCHEMA",
                "SQLite responded with error code %d:\n%s",
                 res,"The database schema changed");
    } else if(res == SQLITE_TOOBIG)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_TOOBIG",
                "SQLite responded with error code %d:\n%s",
                 res,"String or BLOB exceeds size limit");
    } else if(res == SQLITE_CONSTRAINT)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_CONSTRAINT",
                "SQLite responded with error code %d:\n%s",
                 res,"Abort due to constraint violation");
    } else if(res == SQLITE_MISMATCH)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_MISMATCH",
                "SQLite responded with error code %d:\n%s",
                 res,"Data type mismatch");
    } else if(res == SQLITE_MISUSE)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_MISUSE",
                "SQLite responded with error code %d:\n%s",
                 res,"Library used incorrectly");
    } else if(res == SQLITE_NOLFS)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_NOLFS",
                "SQLite responded with error code %d:\n%s",
                 res,"Uses OS features not supported on host");
    } else if(res == SQLITE_AUTH)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_AUTH",
                "SQLite responded with error code %d:\n%s",
                 res,"Authorization denied");
    } else if(res == SQLITE_FORMAT)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_FORMAT",
                "SQLite responded with error code %d:\n%s",
                 res,"Unsupported column type");
    } else if(res == SQLITE_RANGE)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_RANGE",
                "SQLite responded with error code %d:\n%s",
                 res,"2nd parameter to sqlite3_bind out of range");
    } else if(res == SQLITE_NOTADB)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_NOTADB",
                "SQLite responded with error code %d:\n%s",
                 res,"File opened that is not a database file");
    } else if(res == SQLITE_NOTICE)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_NOTICE",
                "SQLite responded with error code %d:\n%s",
                 res,"Notifications from sqlite3_log()");
    } else if(res == SQLITE_WARNING)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_WARNING",
                "SQLite responded with error code %d:\n%s",
                 res,"Warnings from sqlite3_log()");
    } else if(res == SQLITE_ROW)  {
        mexErrMsgIdAndTxt("sqlite3:SQLITE_ROW",
                "SQLite responded with error code %d:\n%s",
                 res,"sqlite3_step() has another row ready");
    } else {
        mexErrMsgIdAndTxt("sqlite3:UnknownError",
                "SQLite responded with error code %d",res);
    }
}

