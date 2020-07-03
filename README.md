This function is a wrapper for mex interfaces that were compiled for most operating systems, and for both Matlab and Octave. All Matlab releases from R14SP3 (v7.1) and later should work. Octave mex files can be compiled from the source, which is downloaded by this function itself.  
Matlab releases older than R14SP3 on Windows (and R2011a on Ubuntu) will use a command line interface (CLI), which imposes several restrictions on syntax and may yield inconsistent results. If the SQL statement returns output, the raw output from the system call is sent as the second output argument to allow custom parsing.  
Note that only the non-CLI Matlab implementations support char values outside of the 0-255 range. If you plan on using Octave or Matlab 6.5 you should make sure the input is valid. Input and output are not sanitized to reflect this, in case it does work as expected.  

A use demo is included.

Sources:  
The Matlab interface is actually mksqlite version 2.5, see SourceForge for the compiled binaries included ([direct link](https://sourceforge.net/projects/mksqlite/files/mksqlite-2.5.zip/download)).  
The Octave mex files included were compiled from the files listed in the help text. The original files can also be downloaded from [here](https://github.com/rmartinjak/mex-sqlite3), [here](http://sqlite.org/2018/sqlite-amalgamation-3230100.zip) and [here](https://github.com/LuaDist/lsqlite3).  
The CLI (command line interface) is in the sqlite-tools-win32-x86-3230100.zip file on sqlite.org ([capture to the Wayback Machine](http://web.archive.org/web/20180515193517if_/http://sqlite.org/2018/sqlite-tools-win32-x86-3230100.zip))  

Licence: CC by-nc-sa 4.0