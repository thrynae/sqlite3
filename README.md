[![View sqlite3 on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/68298-sqlite3)

This function is a wrapper for a command line interface and a mex interface. The mex files were compiled for most operating systems for Matlab. For Octave an internal function will download the source files and compile it, since the resulting mex files are less portable than Matlab mex files.  

Since this is a wrapper for different implementations, it is possible that more complex commands have inconsistent effects depending on whether the mex implementation was used or the CLI.  

The CLI is removed from this version, but will likely be back in a future update. When it is back, the CLI may have inconsistent effects across different operating systems. The test suite will only test basic functionality for the CLI.  

Sources:  
The basis for the interface is the SQLite3 project itself. The sqlite3.c and sqlite3.h files can be downloaded from an archived zip file [here](http://web.archive.org/web/202108id_/https://www.sqlite.org/2021/sqlite-amalgamation-3360000.zip).
The originals for sqlite3_interface.c, structlist.c, and structlist.h can be found on [GitHub](https://github.com/rmartinjak/mex-sqlite3). These 3 files were edited to make them conform to the stricter standards of older compilers and to remove a message: `mexPrintf("binding params %d of %zu\n", i, mxGetM(params));'
Additionally a converter file was written to deal with UTF-16 char encoding.

Licence: CC by-nc-sa 4.0