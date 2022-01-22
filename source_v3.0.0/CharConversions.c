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
/* This file contains conversion functions to convert the UTF-16 MATLAB char  */
/* arrays to UTF-8 C char arrays.                                             */
/* It contains a check to determine whether the runtime uses UTF-16 chars (as */
/* there is a possibility of Octave moving to UTF-16 as well at some point in */
/* the future). This check works by casting the numeric value of the euro     */
/* symbol to char and back to double. If this process is lossless, that means */
/* the runtime uses UTF-16. Any warnings are captured and the previous        */
/* warning state is restored.                                                 */
/*                                                                            */
/* Note that if sizeof(int)*CHAR_BIT is less than 32, this may cause issues   */
/* for 2-word UTF-16 entries.                                                 */
/* Also note that the conversion from UTF-8 C char arrays to UTF-16 MATLAB    */
/* char arrays currently does not work. Any char array will be returned as a  */
/* uint8 and will require conversion outside of this mex routine.             */
/* Both of these issues may be fixed in a future version.                     */
/*                                                                            */
/* -------------------------------------------------------------------------- */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <limits.h> /* This includes CHAR_BIT on Linux Octave. */

#include "mex.h"

/* The global describes what the runtime uses to encode chars:                */
/* UTF-8 (Octave) or UTF-16 (Matlab).                                         */
/*                                                                            */
/* The override variables can be used to force the use of either the          */
/* normal internal functions (set the global to 8) or the replacement         */
/* functions (set the global to 16).                                          */


/* All 3 will be set to 8 or 16 if 0. */
int GLOBAL_RuntimeCharIsUTF=0;
/*Allow granular control over the two conversions.*/
int DEBUG_OVERRIDE_mx_to_char=0;
int DEBUG_OVERRIDE_char_to_mx=0;

/* Wrapper function signatures:                     */
char * WRAPPERmxArrayToUTF8String(const mxArray*);
mxArray * WRAPPERmxCreateCharMatrixFromStrings(int,const char**);

/* Function signatures in this file:                */
char    * mxCharToUTF8(           const mxArray*    );
mxArray * UTF8TomxUTF16char(      const char**      );
int       CountCharsInUTF16(      int*          ,int);
int     * UTF16ToUnicode(         int*          ,int);
int       DetermineUTF8ByteCount( int*          ,int);
int     * AllocStrForUTF16(       mxArray*      ,int);
int     * AllocStrForUTF8(        int*          ,int);
void      DetermineUTFCharRuntime(                  );

/* Wrapper functions: */
char *WRAPPERmxArrayToUTF8String(const mxArray *pa){
    char *str;
    if (GLOBAL_RuntimeCharIsUTF==0){
        DetermineUTFCharRuntime(); }
    if (DEBUG_OVERRIDE_mx_to_char==0){
        DEBUG_OVERRIDE_mx_to_char=GLOBAL_RuntimeCharIsUTF; }
    
    if (DEBUG_OVERRIDE_mx_to_char==16) {
        str=mxCharToUTF8(pa);
    } else {
        str=mxArrayToString(pa);
    }
    return str;
}
mxArray *WRAPPERmxCreateCharMatrixFromStrings(int m, const char **str) {
    if (GLOBAL_RuntimeCharIsUTF==0){
        DetermineUTFCharRuntime(); }
    if (DEBUG_OVERRIDE_char_to_mx==0){
        DEBUG_OVERRIDE_char_to_mx=GLOBAL_RuntimeCharIsUTF; }
    
    if (DEBUG_OVERRIDE_char_to_mx==16) {
        /* Only m=1 is used, so that can be assumed. */
        return UTF8TomxUTF16char(str);
    } else {
        return mxCreateCharMatrixFromStrings(m,str);
    }
}


/* Now the actual functions can be defined: */

mxArray* UTF8TomxUTF16char(const char **pConstStr) {
    /* This function assumes a valid UTF-8 string. */
    int index_UTF8txt;
    int index_UTF16txt;
    int UTF8length;
    int count_UTF8bytes;
    int c_UTF16_w;
    int dummy1;int dummy2;
    int leading_UTF8;
    int UTF8buf[4];
    unsigned int W1;unsigned int W2;
    unsigned int *UTF16text;
    char *UTF8text;const char *ConstStr;
    mxArray *mxUTF16;
    mxArray *rhs[1];
    mxArray *lhs[1];
    
    
    /* Make a copy to remove the const qualifier. */
    ConstStr=*(pConstStr);
    UTF8length=(sizeof (char))*strlen(ConstStr);
    UTF8text=malloc( UTF8length+1 );
    strncpy(UTF8text,ConstStr,strlen(ConstStr));
    UTF8text[UTF8length]='\0';
    
    /* We could also copy the data to a uint8 vector and let Matlab sort out  */
    /* the conversion:                                                        */
    mxUTF16=mxCreateNumericMatrix(1,strlen(UTF8text),mxUINT8_CLASS,mxREAL);
    memcpy(mxGetData(mxUTF16), UTF8text, strlen(UTF8text)+1);
    return mxUTF16;
    /* Below is an attempt to do the conversion in C.                         */

    /* Allocate memory for the UTF16 string (possibly a factor 3 too large). */
    UTF16text=malloc(sizeof(int) * UTF8length);
    
    index_UTF16txt=0;count_UTF8bytes=0;
    for (index_UTF8txt=0;index_UTF8txt<UTF8length;dummy1++) {
        /* Use a while loop to process the bytes of the UTF-8 string. */
        /* If the string is not valid UTF-8, the trailing bytes are   */
        /* set to 0 to prevent a segfault.                            */
        
        /* Increment based on previous iteration. */
        index_UTF8txt=index_UTF8txt+count_UTF8bytes;
        UTF8buf[1]=0;UTF8buf[2]=0;UTF8buf[3]=0;/*clear the buffer*/
        
        leading_UTF8=0xFF&(UTF8text[index_UTF8txt]);
        UTF8buf[0]=leading_UTF8;
        count_UTF8bytes=1;                   /* at least 1 byte */
        if (leading_UTF8>=0xC0){             /* at least 2 bytes*/
            if ((index_UTF8txt+1)>UTF8length){
                UTF8buf[count_UTF8bytes]=0; /* prevent segfault */
            } else {
                UTF8buf[count_UTF8bytes]=0xFF&(UTF8text[index_UTF8txt+1]);
            }
            count_UTF8bytes++;dummy1=UTF8buf[0];dummy2=0x3F;
            UTF8buf[0]=dummy1&dummy2; /*set first 2 bits to 0*/
        }
        if (leading_UTF8>=0xE0){             /* at least 3 bytes*/
            if ((index_UTF8txt+2)>UTF8length){
                UTF8buf[count_UTF8bytes]=0; /* prevent segfault */
            } else {
                UTF8buf[count_UTF8bytes]=0xFF&(UTF8text[index_UTF8txt+2]);
            }
            count_UTF8bytes++;dummy1=UTF8buf[0];dummy2=0x1F;
            UTF8buf[0]=dummy1&dummy2; /*set first 3 bits to 0*/
        }
        if (leading_UTF8>=0xF0){             /* 4 bytes         */
            if ((index_UTF8txt+2)>UTF8length){
                UTF8buf[count_UTF8bytes]=0; /* prevent segfault */
            } else {
                UTF8buf[count_UTF8bytes]=0xFF&(UTF8text[index_UTF8txt+3]);
            }
            count_UTF8bytes++;dummy1=UTF8buf[0];dummy2=0x0F;
            UTF8buf[0]=dummy1&dummy2; /*set first 4 bits to 0*/
        }
        
        /* Now we can determine the Unicode from UTF8 buffer. */
        if (count_UTF8bytes==1) {
            UTF16text[index_UTF16txt]=UTF8buf[0];
            index_UTF16txt++;
        } else if (count_UTF8bytes==2) {
            /* Strip the UTF-8 metadata bits and shift up. */
            dummy1=UTF8buf[0];dummy2=0x3F;UTF8buf[0]=(dummy1&dummy2)<<6;
            dummy1=UTF8buf[1];dummy2=0x3F;UTF8buf[1]=(dummy1&dummy2)   ;

            UTF16text[index_UTF16txt]=(UTF8buf[0])+(UTF8buf[1]);
            index_UTF16txt++;
        } else if (count_UTF8bytes==3) {
            /* Strip the UTF-8 metadata bits and shift up. */
            dummy1=UTF8buf[0];dummy2=0x3F;UTF8buf[0]=(dummy1&dummy2)<<(6*2);
            dummy1=UTF8buf[1];dummy2=0x3F;UTF8buf[1]=(dummy1&dummy2)<<(6  );
            dummy1=UTF8buf[2];dummy2=0x3F;UTF8buf[2]=(dummy1&dummy2)       ;

            UTF16text[index_UTF16txt]=UTF8buf[0]+UTF8buf[1]+UTF8buf[2];
            index_UTF16txt++;
        } else  {
            /* Strip the UTF-8 metadata. */
            dummy1=UTF8buf[0];dummy2=0x3F;UTF8buf[0]=dummy1&dummy2;
            dummy1=UTF8buf[1];dummy2=0x3F;UTF8buf[1]=dummy1&dummy2;
            dummy1=UTF8buf[2];dummy2=0x3F;UTF8buf[2]=dummy1&dummy2;
            dummy1=UTF8buf[3];dummy2=0x3F;UTF8buf[3]=dummy1&dummy2;
            
            /* byte 1: 11110xaa */
            /* byte 2: 10bbbbbb */
            /* byte 3: 10ccdddd */
            /* byte 4: 10eeeeee */
            /* word 1: 110110aabbbbbbcc */
            /* word 2: 110111ddddeeeeee */
            W1=0xD800+(UTF8buf[0]<<8)+(UTF8buf[1]<<2)+(UTF8buf[2]>>4);
            W1=W1-0x40; /* perform U'=U-0x10000 */
            W2=0xDC00+ ((0x0F&(UTF8buf[2]))<<6) + UTF8buf[3];
            
            UTF16text[index_UTF16txt]=W1; index_UTF16txt++;
            UTF16text[index_UTF16txt]=W2; index_UTF16txt++;
        }
    }
    c_UTF16_w=index_UTF16txt-1;
    
    if ( (sizeof(int)*CHAR_BIT)==32 ) {
        /* The first c_UTF16_w elements contain the 32 bit UTF-16 data. */
        mxUTF16=mxCreateNumericMatrix(1,c_UTF16_w,mxUINT32_CLASS,mxREAL);
    } else {
        /* The first c_UTF16_w elements contain the 16 bit UTF-16 data. */
        mxUTF16=mxCreateNumericMatrix(1,c_UTF16_w,mxUINT16_CLASS,mxREAL);
    }
    memcpy(mxGetData(mxUTF16), UTF16text, sizeof(int)*c_UTF16_w);
    rhs[0]=mxUTF16;
    mexCallMATLAB(1,lhs,1,rhs,"char");
    mxUTF16=lhs[0];
    
    free(UTF16text); /* prevent memory leak */
    return mxUTF16;
}

char * mxCharToUTF8_(const mxArray* input) {
    /* This function contains a bug. It results in a segfault, as it          */
    /* sometimes improperly terminates the string, leading to extra           */
    /* characters. However, this function does properly deal with 16 bit int, */
    /* so if it can be fixed, this version would be the better choice.        */
    mxArray *rhs[1];mxArray *lhs[1];
    mxArray *mxChar;
    int mxCharCount,bytes_to_copy,dummy,i_IntStr,i_buffer,curr,next;
    long long_curr,long_next,U,long_buffer[4];
    int *IntStr;
    char *buffer;char *output;
    
    /* First convert the input inside Matlab to make the data copy easier.    */
    /* Technically, this should be uint, but since the sign bit will never be */
    /* used, using signed int is fine. */
    rhs[0]=mxDuplicateArray(input);
    if (        (sizeof(int)*CHAR_BIT)==64 ) {
        mexCallMATLAB(1,lhs,1,rhs,"int64");
    } else if ( (sizeof(int)*CHAR_BIT)==32 ) {
        mexCallMATLAB(1,lhs,1,rhs,"int32");
    } else {
        mexCallMATLAB(1,lhs,1,rhs,"int16");
    }
    mxDestroyArray(rhs[0]); /* Avoid memory leak. */
    mxChar=lhs[0];
    
    
    /* Copy the data from the Matlab variable to a C int array. */
    mxCharCount=mxGetNumberOfElements(mxChar);
    bytes_to_copy=sizeof(int)*mxCharCount;
    IntStr=malloc(bytes_to_copy);
    memcpy(IntStr,mxGetData(mxChar),bytes_to_copy);
    
    /* Allocate an array that will hold the output. It will be too large, but */
    /* we will copy the result to another array and de-allocate this one      */
    /* before returning.                                                      */
    /* A 2-word UTF-16 will convert to a 4 byte UTF-8 (i.e. a factor of 2),   */
    /* so the worst case scenario is that every UTF-16 character is a 3 byte  */
    /* UTF-8 character.                                                       */
    buffer = malloc(sizeof(char)*3*mxCharCount);
    memset(buffer,'\0',sizeof(buffer)); /* Null-terminate string. */
    
    dummy=0;
    i_buffer=(-1); /* Start off negative (this indexes the buffer). */
    /* Now we loop through the IntStr. */
    for (i_IntStr=0;i_IntStr<mxCharCount;i_IntStr++) {
        /* UTF-8 encoding scheme:                                         */
        /* 1 byte: 0xxxxxxx                                               */
        /* 2 byte: 110xxxxx 10xxxxxx                                      */
        /* 3 byte: 1110xxxx 10xxxxxx 10xxxxxx                             */
        /* 4 byte: 11110xxx 10xxxxxx 10xxxxxx 10xxxxxx                    */
        /* Note that UTF-16 doesn't encode a code point U, but U-0x10000. */
        curr=IntStr[i_IntStr];
        if (curr>0xD7FF && curr<0xE000) {                        /*4 bytes*/
            /* This is a surrogate pair. Load the second one as well. */
            if (i_IntStr==mxCharCount) { next=0;i_IntStr++;}
            else { next=IntStr[i_IntStr];i_IntStr++; }
            
            /* Convert everything to a long so we have room for UTF-32. */
            long_curr=(long)curr;long_next=(long)next;
            /* Keep only encoding bits. */
            /* UTF-16: 110110aabbbbbbcc 110111ddddeeeeee */
            long_curr=long_curr&0x3FF;long_next=long_next&0x3FF;
            /* Calculate the UTF-32 value. */
            U=(long_curr*0x400) + long_next + 0x10000;
            
            /* Bitshift to determine the UTF-8 encoding bits. */
            long_buffer[0]= (U>>18)&0x03;
            long_buffer[1]= (U>>12)&0x3F;
            long_buffer[2]= (U>>6 )&0x3F;
            long_buffer[3]= (U    )&0x3F;
            
            /* Add leading bits and store in output buffer. */
            buffer[i_buffer+1]=0xF0 + (int)long_buffer[0];
            buffer[i_buffer+2]=0x80 + (int)long_buffer[1];
            buffer[i_buffer+3]=0x80 + (int)long_buffer[2];
            buffer[i_buffer+4]=0x80 + (int)long_buffer[3];
            i_buffer=i_buffer+4;
        } else if (curr>0x07FF){                                 /*3 bytes*/
            buffer[i_buffer+1]=0xE0+((curr>>12)&0x0F);
            buffer[i_buffer+2]=0x80+((curr>>6 )&0x3F);
            buffer[i_buffer+3]=0x80+( curr     &0x3F);
            i_buffer=i_buffer+3;
        } else if (curr>0x007F){                                 /*2 bytes*/
            buffer[i_buffer+1]=0xC0+((curr>>6 )&0x1F);
            buffer[i_buffer+2]=0x80+( curr     &0x3F);
            i_buffer=i_buffer+2;
        } else {                                                 /*1 byte*/
            buffer[i_buffer+1]=curr;
            i_buffer++;
        }
    }

    /* Copy the contents of the buffer to the output and return. */
    output = malloc(sizeof(char)*(i_buffer+2));
    memset(output,'\0',sizeof(output)); /* Null-terminate string. */
    strcpy(output,buffer);
    
    free(buffer);free(IntStr); /* Avoid memory leak. */
    return output;
}

void DetermineUTFCharRuntime() {
    /* This checks if 8364==double(char(8364)). */
    /* The value of the global will be set to either 8 or 16. */
    mxArray *lhs_0[1];  mxArray *rhs_0[1];
    mxArray *lhs_1[1];  mxArray *rhs_1[1];
    mxArray *lhs_2[2];  mxArray *rhs_2[2];
    mxArray *w;mxArray *wmsg;mxArray *wID;mxArray *c;mxArray *d;mxArray *tf;
    mxArray *L1;mxArray *L2;mxArray *L3;
    
    /* literals:                              /* L1='off';L2='all';L3=8364; */
    L1=mxCreateString("off");
    L2=mxCreateString("all");
    L3=mxCreateDoubleScalar(8364);
    
    lhs_1[0]=w;rhs_2[0]=L1;rhs_2[1]=L2;
    mexCallMATLAB(1,lhs_1,2,rhs_2,"warning"); /* w=warning('off','all'); */
    w=lhs_1[0];
    
    lhs_2[0]=wmsg;lhs_2[1]=wID;
    mexCallMATLAB(2,lhs_2,0,rhs_0,"lastwarn");/* [wmsg,wID]=lastwarn; */
    wmsg=lhs_2[0];wID=lhs_2[1];
    
    lhs_1[0]=c;rhs_1[0]=L3;
    mexCallMATLAB(1,lhs_1,1,rhs_1,"char");    /* c=char(L3); */
    c=lhs_1[0];
    
    lhs_1[0]=d;rhs_1[0]=c;
    mexCallMATLAB(1,lhs_1,1,rhs_1,"double");  /* d=double(c); */
    d=lhs_1[0];
    
    lhs_1[0]=tf;rhs_2[0]=L3;rhs_2[1]=d;
    mexCallMATLAB(1,lhs_1,2,rhs_2,"isequal"); /* tf=isequal(L3,d); */
    tf=lhs_1[0];
    
    rhs_1[0]=w;
    mexCallMATLAB(0,lhs_0,1,rhs_1,"warning"); /* warning(w); */
    
    rhs_2[0]=wmsg;rhs_2[1]=wID;
    mexCallMATLAB(0,lhs_0,2,rhs_2,"isequal"); /* lastwarn(wmsg,wID); */
    
    
    if (0==(int)mxGetScalar(tf)) {
        GLOBAL_RuntimeCharIsUTF=8;
    } else {
        GLOBAL_RuntimeCharIsUTF=16;
    }
    
    /* Do the housekeeping to avoid memory leaks. */
    mxDestroyArray(w);
    mxDestroyArray(wmsg);
    mxDestroyArray(wID);
    mxDestroyArray(c);
    mxDestroyArray(d);
    mxDestroyArray(tf);
    mxDestroyArray(L1);
    mxDestroyArray(L2);
    mxDestroyArray(L3);
}

char * mxCharToUTF8(const mxArray* input) {
    mxArray *rhs[1];mxArray *lhs[1];
    mxArray *mxChar;
    int *CharStr;int *UnicodeStr;int *UTF8Str;
    char *buffer;
    char *output;
    char SingleCharBuffer[40];
    int UnicodeCount,mxCharCount,UTF8Count,i;
    
    /* First convert to int32 inside Matlab to make the data copy easier. */
    rhs[0]=mxDuplicateArray(input);
    mexCallMATLAB(1,lhs,1,rhs,"int32");
    mxDestroyArray(rhs[0]); /* Avoid memory leak. */
    mxChar=lhs[0];
    
    /* Create an int array with the Unicode code point. */
    mxCharCount=mxGetNumberOfElements(mxChar);
    CharStr=AllocStrForUTF16(mxChar,mxCharCount);
    UnicodeCount=CountCharsInUTF16(CharStr,mxCharCount);
    UnicodeStr=UTF16ToUnicode(CharStr,UnicodeCount);
    free(CharStr); /* No longer needed. */
    
    /* Create an int array with the UTF8 encoded characters. */
    UTF8Count=DetermineUTF8ByteCount(UnicodeStr,UnicodeCount);
    UTF8Str=AllocStrForUTF8(UnicodeStr,UnicodeCount);
    free(UnicodeStr); /* No longer needed. */
    
    buffer=malloc(sizeof(char) * (UTF8Count+2) ); /* +2 to NULL-terminate */
    output=malloc(sizeof(char) * (UTF8Count+2) ); /* +2 to NULL-terminate */
    memset(output,'\0',UTF8Count+2);
    for (i=0;i<UTF8Count;i++) {
        sprintf(SingleCharBuffer,"%c",UTF8Str[i]);
        buffer[i]=SingleCharBuffer[0];
    }
    strncpy(output,buffer,UTF8Count);
    free(buffer);
    return output;
}

int CountCharsInUTF16(int *str, int ElemCount) {
    int count,curr;
    size_t i,n;
    count=0;
    for (i = 0; i < ElemCount; i++) {
        curr=str[i];count++;
        if (curr < 0xD7FF || curr > 0xE000) {
            /* 0000...D7FF and E000...FFFF are single-element characters. */
        } else {
            i++;
        }
    }
    return count;
}

int * UTF16ToUnicode(int *str, int CharCount) {
    int W1, W2, curr, i, i2;
    int *u;
    
    u = malloc (sizeof (int) * CharCount);
    
    i2=0;/*i2 indexes the str array, i indexes the Unicode array*/
    for (i = 0; i < CharCount; i++) {
        /* 0000...D7FF and E000...FFFF are single-element characters. */
        if (str[i2] < 0xD7FF || str[i2] > 0xE000) {
            /* single-word character */
            curr=str[i2];
            i2++;
        } else {
            /* parse two-word character */
            W1=str[i2]-0xD800;  i2++;
            W2=str[i2]-0xDC00;  i2++;
            curr=(W1*0x400)+W2+0x10000;
        }
        u[i]=curr;
    }
    return u;
}

int DetermineUTF8ByteCount(int *UnicodeStr,int UnicodeCount) {
    int i,UTF8ByteCount=0;
    /* Count number of bytes required for the given Unicode string. */
    for (i = 0; i < UnicodeCount; i++) {
        UTF8ByteCount++;                           /* at least 1 byte */
        if(UnicodeStr[i]>0x007F){UTF8ByteCount++;} /* at least 2 bytes*/
        if(UnicodeStr[i]>0x07FF){UTF8ByteCount++;} /* at least 3 bytes*/
        if(UnicodeStr[i]>0xFFFF){UTF8ByteCount++;} /* 4 bytes         */
    }
    return UTF8ByteCount;
}

int * AllocStrForUTF16(mxArray *str, int ElemCount) {
    int i,bytes_to_copy;
    int *p,*u;
    bytes_to_copy = sizeof (int) * ElemCount;
    u = malloc (bytes_to_copy);
    memcpy(u,mxGetData(str),bytes_to_copy);
    return u;
}

int * AllocStrForUTF8(int *UnicodeStr,int UnicodeCount) {
    /* Allocate memory for the UTF-8 string and convert from Unicode. */
    int * str;
    int UTF8ByteCount,iOut,iIn,curr,Last6Bits;
    
    /* Allocate memory for the output array. */
    str=malloc(sizeof(int)*DetermineUTF8ByteCount(UnicodeStr,UnicodeCount)+1);
    
    /* Do the conversion to UTF8 here. */
    iOut=0;
    for (iIn = 0; iIn < UnicodeCount; iIn++) {
        curr=UnicodeStr[iIn];
        
        if        (curr>0xFFFF){ /*4 bytes*/
            /*trailing bytes*/
            Last6Bits=curr&0x3F;curr=curr-Last6Bits;curr=curr>>6;
            str[iOut+3]=Last6Bits+0x80;
            Last6Bits=curr&0x3F;curr=curr-Last6Bits;curr=curr>>6;
            str[iOut+2]=Last6Bits+0x80;
            Last6Bits=curr&0x3F;curr=curr-Last6Bits;curr=curr>>6;
            str[iOut+1]=Last6Bits+0x80;
            /*leading byte*/
            str[iOut]=curr+0xF0;
            iOut=iOut+4;
        } else if (curr>0x07FF){ /*3 bytes*/
            /*trailing bytes*/
            Last6Bits=curr&0x3F;curr=curr-Last6Bits;curr=curr>>6;
            str[iOut+2]=Last6Bits+0x80;
            Last6Bits=curr&0x3F;curr=curr-Last6Bits;curr=curr>>6;
            str[iOut+1]=Last6Bits+0x80;
            /*leading byte*/
            str[iOut]=curr+0xE0;
            iOut=iOut+3;
        } else if (curr>0x007F){ /*2 bytes*/
            /*trailing byte*/
            Last6Bits=curr&0x3F;curr=curr-Last6Bits;curr=curr>>6;
            str[iOut+1]=Last6Bits+0x80;
            /*leading byte*/
            str[iOut]=curr+0xC0;
            iOut=iOut+2;
        } else {                 /*1 byte*/
            str[iOut]=curr;
            iOut=iOut+1;
        }
    }
    return str;
}

