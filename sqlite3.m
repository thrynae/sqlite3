function [response_struct,cli_output]=sqlite3(db_or_opts,SQL_statement,input_struct)
%Matlab and Octave interface to the SQLite engine
%
% With SQLite you can access an SQL database, without the need of a server. This function is a
% wrapper for a command line interface and a mex interface. The mex files were compiled for most
% operating systems for Matlab. For Octave an internal function will download the source files and
% compile it, since the resulting mex files are less portable than Matlab mex files.
% Since this is a wrapper for different implementations, it is possible that more complex commands
% have inconsistent effects depending on whether the mex implementation was used or the CLI.
%
% The CLI is removed from this version, but will likely be back in a future update. When it is
% back, the CLI may have inconsistent effects across different operating systems. The test suite
% will only test basic functionality for the CLI.
%
% Syntax:
%   sqlite3(filename,SQL_statement)
%   sqlite3(optionstruct,SQL_statement)
%   sqlite3(___,input_struct)
%   response_struct=sqlite3(___)
%   [response_struct,cli_output]=sqlite3(___)
%
% Input/output arguments:
% response_struct:
%   This is an Nx1 struct array, with each field name corresponding to the column name in the
%   selected table. Other commands might result in different types of output. A very limited number
%   of commands and syntax options has been tested.
%   For CLI interactions, an attempt will be made to parse the returned output to a struct and to
%   convert to an appropriate data type. For fields where this fails, the fields will be char
%   arrays. The "PRAGMA table_info('table name')" command is used to parse the data types of
%   output_struct for CLI calls.
% cli_output_string:
%   This is a char vector containing the raw output returned by the system function.
%   This variable defaults to an empty char array if not available.
%   NB: the current version of this function does not support a CLI interaction, so this variable
%   will always be an empty char. It is only here to avoid syntax errors for users updating this
%   function.
% filename:
%   The file name of the database file to be edited. Using a full or relative path is optional.
% optionstruct:
%   A struct with the parameters. Any missing field will be filled from the defaults listed below.
%   Using incomplete parameter names or incorrect capitalization is allowed, as long as there is a
%   unique match.
%   NB: with this syntax the filename must be included as a parameter.
% SQL_statement:
%   A char vector containing a valid SQL statement. Please note that the validity of the SQL
%   statement is mostly determined by the mex implementation used. There may be differences between
%   implementations.
% input_struct:
%   If the SQL statement contains placeholders, the values are filled from this parameter. This
%   input argument must be a struct and the number of fields should match the number of
%   placeholders.
%   The number of unnamed parameters ('?') must be equal to the number of fields left after
%   matching named (:NAME, @NAME, $NAME) or indexed (?NNN) parameters.
%
% Name,Value parameters:
%   filename:
%      Since the optionstruct replaces the filename, it has to be entered as an optional parameter.
%      Not specifying this will result in an error.
%      [default=error;]
%   use_CLI:
%      Use a system-installed command line interface. This will result in an error if no CLI is
%      available. For this version of the function the CLI has been removed, so this will error if
%      you set this flag to true. The CLI will probably be back in a future update.
%      On prior versions on Windows an exe file was downloaded, while Linux assumed a system
%      install. The complete statement was written to a file in the tempdir. There may be
%      additional limitations compared to mex implementations. [default=false;]
%   print_to_con:
%      An attempt is made to also use this parameter for warnings or errors during input parsing.
%      A logical that controls whether warnings and other output will be printed to the command
%      window. Errors can't be turned off. [default=true;] if either print_to_fid, print_to_obj, or
%      print_to_fcn is specified then [default=false]
%   print_to_fid:
%      An attempt is made to also use this parameter for warnings or errors during input parsing.
%      The file identifier where console output will be printed. Errors and warnings will be
%      printed including the call stack. You can provide the fid for the command window (fid=1) to
%      print warnings as text. Errors will be printed to the specified file before being actually
%      thrown. [default=[];]
%      If print_to_fid, print_to_obj, and print_to_fcn are all empty, this will have the effect of
%      suppressing every output except errors. 
%      Array inputs are allowed.
%   print_to_obj:
%      An attempt is made to also use this parameter for warnings or errors during input parsing.
%      The handle to an object with a String property, e.g. an edit field in a GUI where console
%      output will be printed. Messages with newline characters (ignoring trailing newlines) will
%      be returned as a cell array. This includes warnings and errors, which will be printed
%      without the call stack. Errors will be written to the object before the error is actually
%      thrown. [default=[];]
%      If print_to_fid, print_to_obj, and print_to_fcn are all empty, this will have the effect of
%      suppressing every output except errors.
%      Array inputs are allowed.
%   print_to_fcn:
%      An attempt is made to also use this parameter for warnings or errors during input parsing.
%      A struct with a function handle, anonymous function or inline function in the 'h' field and
%      optionally additional data in the 'data' field. The function should accept three inputs: a
%      char array (either 'warning' or 'error'), a struct with the message, id, and stack, and the
%      optional additional data. The function(s) will be run before the error is actually thrown.
%      [default=[];]
%      If print_to_fid, print_to_obj, and print_to_fcn are all empty, this will have the effect of
%      suppressing every output except errors.
%      Array inputs are allowed.
%
% The basis for the interface is the SQLite3 project itself. The sqlite3.c and sqlite3.h files
% can be downloaded from an archived zip file here:
% http://web.archive.org/web/202108id_/https://www.sqlite.org/2021/sqlite-amalgamation-3360000.zip
% The originals for sqlite3_interface.c, structlist.c, and structlist.h can be found on GitHub
% (https://github.com/rmartinjak/mex-sqlite3). These 3 files were edited to make them conform to
% the stricter standards of older compilers and to remove a message ('mexPrintf("binding params %d
% of %zu\n", i, mxGetM(params));').
% Additionally a converter file was written to deal with UTF-16 char encoding.
%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 3.0.0                                                         |%
%|  Date:    2022-01-22                                                    |%
%|  Author:  H.J. Wisselink                                                |%
%|  Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 ) |%
%|  Email = 'h_j_wisselink*alumnus_utwente_nl';                            |%
%|  Real_email = regexprep(Email,{'*','_'},{'@','.'})                      |%
%|                                                                         |%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%
% Tested on several versions of Matlab (ML 6.5 and onward) and Octave (4.4.1 and onward), and on
% multiple operating systems (Windows/Ubuntu/MacOS). For the full test matrix, see the HTML doc.
% Compatibility considerations:
% - This is expected to work on most releases. Forward compatibility of Matlab mex files is
%   excellent, but backward compatibility is often an issue. If you use an OS/realease-combination
%   for which the current binaries don't work, you're welcome to send me a compiled mex file by
%   email, so I can add it to a future update.

%Verify number of input arguments
if nargout>=2
    cli_output='';%set default
end
if nargin<2 || nargin>3
    error('HJW:sqlite3:nargin','Incorrect number of input argument.')
end
if nargout>2
    error('HJW:sqlite3:nargout','Incorrect number of output argument.')
end

% Linearize a struct array input so later rows are processed as well.
args={db_or_opts,SQL_statement};if nargin>=3,args{end+1}=reshape(input_struct,1,[]);end
[success,opts,ME]=sqlite3_parse_inputs(args{:});
if ~success
    error_(opts.print_to,ME)
else
    print_to=opts.print_to;
    if nargin==3
        args={opts.filename,opts.SQL_statement,opts.input_struct};
    else
        args={opts.filename,opts.SQL_statement};
    end
end

%Determine the function call type.
persistent fun_handle
if isempty(fun_handle)
    % Store the handle to the correct mex file in a persistent variable.
    % Octave requires a compile specific to at least the release and architecture, so it will be
    % compiled from the source files. If any files need to be downloaded, this will require
    % internet access.
    fun_handle=get_sqlite3_mex(print_to);
end

% Actually call the mex implementation.
try ME=[]; %#ok<NASGU>
    c=cell(1,logical(nargout));
    [c{:}]=feval(fun_handle,args{:});
    if nargout,response_struct=c{1};end
catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
    error_(print_to,ME)
end

% A normal clear would wipe the persistent.
if nargout==0                    ,return,end
if ~isa(response_struct,'struct'),return,end

% On runtimes where chars are encoded as UTF-16, text is returned as uint8. That means any text
% should be re-encoded. While it is technically possible to do this inside the mex implementation,
% the code to do so seems to have a bug. Any help in fixing this/these bug/bugs is welcome.
if ~CharIsUTF8
    fields=fieldnames(response_struct);
    for n=1:numel(fields)
        for m=1:numel(response_struct)
            x=response_struct(m).(fields{n});
            if ~isa(x,'uint8'),continue,end
            % Attempt to convert to UTF-16. If the content is not valid UTF-8, this will not
            % trigger an error, but return the content as uint8. This should not happen.
            [U,tf]=UTF8_to_unicode(x);
            if tf,response_struct(m).(fields{n})=unicode_to_char(U);end
        end
    end
end
end
function out=bsxfun_plus(in1,in2)
%Implicit expansion for plus(), but without any input validation.
try
    out=in1+in2;
catch
    try
        out=bsxfun(@plus,in1,in2);
    catch
        sz1=size(in1);                    sz2=size(in2);
        in1=repmat(in1,max(1,sz2./sz1));  in2=repmat(in2,max(1,sz1./sz2));
        out=in1+in2;
    end
end
end
function c=char2cellstr(str,LineEnding)
%Split char or uint32 vector to cell (1 cell element per line). Default splits are for CRLF/CR/LF.
%The input data type is preserved.
%
%Since the largest valid Unicode codepoint is 0x10FFFF (i.e. 21 bits), all values will fit in an
%int32 as well. This is used internally to deal with different newline conventions.
%
%The second input is a cellstr containing patterns that will be considered as newline encodings.
%This will not be checked for any overlap and will be processed sequentially.

returnChar=isa(str,'char');
str=int32(str);%convert to signed, this should not crop any valid Unicode codepoints.

if nargin<2
    %Replace CRLF, CR, and LF with -10 (in that order). That makes sure that all valid encodings of
    %newlines are replaced with the same value. This should even handle most cases of files that
    %mix the different styles, even though such mixing should never occur in a properly encoded
    %file. This considers LFCR as two line endings.
    if any(str==13)
        str=PatternReplace(str,int32([13 10]),int32(-10));
        str(str==13)=-10;
    end
    str(str==10)=-10;
else
    for n=1:numel(LineEnding)
        str=PatternReplace(str,int32(LineEnding{n}),int32(-10));
    end
end

%Split over newlines.
newlineidx=[0 find(str==-10) numel(str)+1];
c=cell(numel(newlineidx)-1,1);
for n=1:numel(c)
    s1=(newlineidx(n  )+1);
    s2=(newlineidx(n+1)-1);
    c{n}=str(s1:s2);
end

%Return to the original data type.
if returnChar
    for n=1:numel(c),c{n}=  char(c{n});end
else
    for n=1:numel(c),c{n}=uint32(c{n});end
end
end
function tf=CharIsUTF8
%This provides a single place to determine if the runtime uses UTF-8 or UTF-16 to encode chars.
%The advantage is that there is only 1 function that needs to change if and when Octave switches to
%UTF-16. This is unlikely, but not impossible.
persistent persistent_tf
if isempty(persistent_tf)
    if exist('OCTAVE_VERSION','builtin')~=0 %Octave
        %Test if Octave has switched to UTF-16 by looking if the Euro symbol is losslessly encoded
        %with char.
        w=struct('w',warning('off','all'));[w.msg,w.ID]=lastwarn;
        persistent_tf=~isequal(8364,double(char(8364)));
        warning(w.w);lastwarn(w.msg,w.ID); % Reset warning state.
    else
        persistent_tf=false;
    end
end
tf=persistent_tf;
end
function date_correct=check_date(outfilename,opts,SaveAttempt)
%Check if the date of the downloaded file matches the requested date.
%
%There are two strategies. Strategy 1 is guaranteed to be correct, but isn't always possible.
%Strategy 2 could give an incorrect answer, but is possible in more situations. In the case of
%non-web page files (like e.g. an image), both will fail. This will trigger a missing date error,
%for which you need to input a missing date response (m_date_r).
%
%Strategy 1:
%Rely on the html for the header to provide the date of the currently viewed capture.
%Strategy 2:
%Try a much less clean version: don't rely on the top bar, but look for links that indicate a link
%to the same date in the Wayback Machine. The most common occurring date will be compared with
%date_part.

if ~exist(outfilename,'file')
    date_correct=false;return
    %If the file doesn't exist (not even as a 0 byte file), evidently something went wrong, so
    %retrying or alerting the user is warranted.
end
[m_date_r,date_bounds,print_to]=deal(opts.m_date_r,opts.date_bounds,opts.print_to);
%Loading an unsaved page may result in a capture of the live page (but no save in the WBM). If this
%happens the time in the file will be very close to the current time if this is the case. If the
%save was actually triggered this is valid, but if this is the result of a load attempt, it is
%unlikely this is correct, in which case it is best to trigger the response to an incorrect date:
%attempt an explicit save.
%Save the time here so any time taken up by file reading and processing doesn't bias the estimation
%of whether or not this is too recent.
if ~SaveAttempt
    currentTime=WBM_getUTC_local;
end

%Strategy 1:
%Rely on the html for the header to provide the date of the currently viewed capture.
StringToMatch='<input type="hidden" name="date" value="';
data=readfile(outfilename);
%A call to ismember would be faster, but it can result in a memory error in ML6.5. The
%undocumented ismembc function only allows numeric, logical, or char inputs (and Octave lacks it),
%so we can't use that on our cellstr either. That is why we need the while loop here.
pos=0;
while pos<=numel(data) && (pos==0 || ~strcmp(stringtrim(data{pos}),'<td class="u" colspan="2">'))
    %equivalent to pos=find(ismember(data,'<td class="u" colspan="2">'));
    pos=pos+1;
end
if numel(data)>=(pos+1)
    line=data{pos+1};
    idx=strfind(line,StringToMatch);
    idx=idx+length(StringToMatch)-1;
    date_as_double=str2double(line(idx+(1:14)));
    date_correct= date_bounds.double(1)<=date_as_double && date_as_double<=date_bounds.double(2);
    return
end
%Strategy 2:
%Try a much less clean version: don't rely on the top bar, but look for links that indicate a link
%to the same date in the Wayback Machine. The most common occurring date will be compared with
%date_part.
%The file was already loaded with data=readfile(outfilename);
data=data(:)';data=cell2mat(data);
%The data variable is now a single long string.
idx=strfind(data,'/web/');
if numel(idx)==0
    if m_date_r==0     %ignore
        date_correct=true;
        return
    elseif m_date_r==1 %warning
        warning_(print_to,'HJW:WBM:MissingDateWarning',...
            'No date found in file, unable to check date, assuming it is correct.')
        date_correct=true;
        return
    elseif m_date_r==2 %error
        error_(print_to,'HJW:WBM:MissingDateError',...
            ['Could not find date. This can mean there is an ',...
            'error in the save. Try saving manually.'])
    end
end
datelist=zeros(size(idx));
data=[data 'abcdefghijklmnopqrstuvwxyz'];%Avoid error in the loop below.
if exist('isstrprop','builtin')
    for n=1:length(idx)
        for m=1:14
            if ~isstrprop(data(idx(n)+4+m),'digit')
                break
            end
        end
        datelist(n)=str2double(data(idx(n)+4+(1:m)));
    end
else
    for n=1:length(idx)
        for m=1:14
            if ~any(double(data(idx(n)+4+m))==(48:57))
                break
            end
        end
        datelist(n)=str2double(data(idx(n)+4+(1:m)));
    end
end
[a,ignore_output,c]=unique(datelist);%#ok<ASGLU> ~
%In some future release, histc might not be supported anymore.
try
    [ignore_output,c2]=max(histc(c,1:max(c)));%#ok<HISTC,ASGLU>
catch
    [ignore_output,c2]=max(accumarray(c,1)); %#ok<ASGLU>
end
date_as_double=a(c2);
date_correct= date_bounds.double(1)<=date_as_double && date_as_double<=date_bounds.double(2);

if ~SaveAttempt
    %Check if the time in the file is too close to the current time to be an actual loaded capture.
    %Setting this too high will result in too many save triggers, but setting it too low will lead
    %to captures being missed on slower systems/networks. 15 seconds seems a reasonable middle
    %ground.
    %One extreme situation to be aware of: it is possible for a save to be triggered, the request
    %arrives successfully and the page is saved, but the response from the server is wrong or
    %missing, triggering an HTTP error. This may then lead to a load attempt. Now we have the
    %situation where there is a save of only several seconds old, but the the SaveAttempt flag is
    %false. The time chosen here must be short enough to account for this situation.
    %Barring such extreme circumstances, page ages below a minute are suspect.
    
    if date_as_double<1e4%Something is wrong
        %Trigger missing date response. This shouldn't happen, so offer a graceful exit.
        if m_date_r==0     %ignore
            date_correct=true;
            return
        elseif m_date_r==1 %warning
            warning_(print_to,'HJW:WBM:MissingDateWarning',...
                'No date found in file, unable to check date, assuming it is correct.')
            date_correct=true;
            return
        elseif m_date_r==2 %error
            error_(print_to,'HJW:WBM:MissingDateError',...
                ['Could not find date. This can mean there is an error in the save.',...
                char(10),'Try saving manually.']) %#ok<CHARTEN>
        end
    end
    
    %convert the date found to a format that the ML6.5 datenum supports
    line=sprintf('%014d',date_as_double);
    line={line(1:4),line(5:6),line(7:8),...  %date
        line(9:10),line(11:12),line(13:14)}; %time
    line=str2double(line);
    timediff=(currentTime-datenum(line))*24*60*60;
    if timediff<10%seconds
        date_correct=false;
    elseif timediff<60%seconds
        warning_(print_to,'HJW:WBM:LivePageStored',...
            ['The live page might have been saved instead of a capture.',char(10),...
            'Check on the WBM if a capture exists.']) %#ok<CHARTEN>
    end
end
end
function outfilename=check_filename(filename,outfilename)
%It can sometimes happen that the outfilename provided by urlwrite is incorrect. Therefore, we need
%to check if either the outfilename file exists, or the same file, but inside the current
%directory. It is unclear when this would happen, but it might be that this only happens when the
%filename provided only contains a name, and not a full or relative path.
outfilename2=[pwd filesep filename];
if ~strcmp(outfilename,outfilename2) && ~exist(outfilename,'file') && exist(outfilename2,'file')
	outfilename=outfilename2;
end
end
function [tf,ME]=CheckMexCompilerExistence
% Returns true if a mex compiler is expected to be installed.
% The method used for R2008a and later is fairly slow, so the flag is stored in a file. Run
% ClearMexCompilerExistenceFlag() to reset this test.
%
% This function may result in false positives (e.g. by detecting an installed compiler that doesn't
% work, or if a compiler is required for a specific language).
% False negatives should be rare.
%
% Based on: http://web.archive.org/web/2/http://www.mathworks.com/matlabcentral/answers/99389
% (this link will redirect to the URL with the full title)
%
% The actual test will be performed in a separate function. That way the same persistent can be
% used for different functions containing this check as a subfunction.

persistent tf_ ME_
if isempty(tf_)
    fun_handle='CheckMexCompilerExistence_persistent';

    %In some release-runtime combinations addpath has a permanent effect, in others it doesn't. By
    %putting this code in this block, we are trying to keep these queries to a minimum.
    [p,ME]=CreatePathFolder__CheckMexCompilerExistence_persistent;
    if ~isempty(ME),tf=false;return,end

    fn=fullfile(p,'ClearMexCompilerExistenceFlag.m');
    txt={...
        'function ClearMexCompilerExistenceFlag';
        'fn=create_fn';
        'if exist(fn,''file''),delete(fn),end';
        'end';
        'function fn=create_fn';
        'v=version;v=v(regexp(v,''[a-zA-Z0-9()\.]''));';
        'if ~exist(''OCTAVE_VERSION'', ''builtin'')';
        '    runtime=''MATLAB'';';
        '    type=computer;';
        'else';
        '    runtime=''OCTAVE'';';
        '    arch=computer;arch=arch(1:(min(strfind(arch,''-''))-1));';
        '    if ispc';
        '        if strcmp(arch,''x86_64'')  ,type= ''win_64'';';
        '        elseif strcmp(arch,''i686''),type= ''win_i686'';';
        '        elseif strcmp(arch,''x86'') ,type= ''win_x86'';';
        '        else                      ,type=[''win_'' arch];';
        '        end';
        '    elseif isunix && ~ismac %Essentially this is islinux';
        '        if strcmp(arch,''i686'')      ,type= ''lnx_i686'';';
        '        elseif strcmp(arch,''x86_64''),type= ''lnx_64'';';
        '        else                        ,type=[''lnx_'' arch];';
        '        end';
        '    elseif ismac';
        '        if strcmp(arch,''x86_64''),type= ''mac_64'';';
        '        else                    ,type=[''mac_'' arch];';
        '        end';
        '    end';
        'end';
        'type=strrep(strrep(type,''.'',''''),''-'','''');';
        'flag=[''flag_'' runtime ''_'' v ''_'' type ''.txt''];';
        'fn=fullfile(fileparts(mfilename(''fullpath'')),flag);';
        'end';
        ''};
    fid=fopen(fn,'wt');fprintf(fid,'%s\n',txt{:});fclose(fid);

    fn=fullfile(p,[fun_handle '.m']);
    txt={...
        'function [tf,ME]=CheckMexCompilerExistence_persistent';
        '% Returns true if a mex compiler is expected to be installed.';
        '% The method used for R2008a and later is fairly slow, so the flag is stored in';
        '% a file. Run ClearMexCompilerExistenceFlag() to reset this test.';
        '%';
        '% This function may result in false positives (e.g. by detecting an installed';
        '% compiler that doesn''t work, or if a compiler is required for a specific';
        '% language). False negatives should be rare.';
        '%';
        '% Based on: http://web.archive.org/web/2/';
        '% http://www.mathworks.com/matlabcentral/answers/99389';
        '% (this link will redirect to the URL with the full title)';
        '';
        'persistent tf_ ME_';
        'if isempty(tf_)';
        '    ME_=create_ME;';
        '    fn =create_fn;';
        '    if exist(fn,''file'')';
        '        str=fileread(fn);';
        '        tf_=strcmp(str,''compiler found'');';
        '    else';
        '        % Use evalc to suppress anything printed to the command window.';
        '        [txt,tf_]=evalc(func2str(@get_tf)); %#ok<ASGLU>';
        '        fid=fopen(fn,''w'');';
        '        if tf_,fprintf(fid,''compiler found'');';
        '        else , fprintf(fid,''compiler not found'');end';
        '        fclose(fid);';
        '    end';
        '    ';
        'end';
        'tf=tf_;ME=ME_;';
        'end';
        'function fn=create_fn';
        'v=version;v=v(regexp(v,''[a-zA-Z0-9()\.]''));';
        'if ~exist(''OCTAVE_VERSION'', ''builtin'')';
        '    runtime=''MATLAB'';';
        '    type=computer;';
        'else';
        '    runtime=''OCTAVE'';';
        '    arch=computer;arch=arch(1:(min(strfind(arch,''-''))-1));';
        '    if ispc';
        '        if strcmp(arch,''x86_64'')  ,type= ''win_64'';';
        '        elseif strcmp(arch,''i686''),type= ''win_i686'';';
        '        elseif strcmp(arch,''x86'') ,type= ''win_x86'';';
        '        else                      ,type=[''win_'' arch];';
        '        end';
        '    elseif isunix && ~ismac %Essentially this is islinux';
        '        if strcmp(arch,''i686'')      ,type= ''lnx_i686'';';
        '        elseif strcmp(arch,''x86_64''),type= ''lnx_64'';';
        '        else                        ,type=[''lnx_'' arch];';
        '        end';
        '    elseif ismac';
        '        if strcmp(arch,''x86_64''),type= ''mac_64'';';
        '        else                    ,type=[''mac_'' arch];';
        '        end';
        '    end';
        'end';
        'type=strrep(strrep(type,''.'',''''),''-'','''');';
        'flag=[''flag_'' runtime ''_'' v ''_'' type ''.txt''];';
        'fn=fullfile(fileparts(mfilename(''fullpath'')),flag);';
        'end';
        'function ME_=create_ME';
        'msg={...';
        '    ''No selected compiler was found.'',...';
        '    ''Please make sure a supported compiler is installed and set up.'',...';
        '    ''Run mex(''''-setup'''') for version-specific documentation.'',...';
        '    '''',...';
        '    ''Run ClearMexCompilerExistenceFlag() to reset this test.''};';
        'msg=sprintf(''\n%s'',msg{:});msg=msg(2:end);';
        'ME_=struct(...';
        '    ''identifier'',''HJW:CheckMexCompilerExistence:NoCompiler'',...';
        '    ''message'',msg);';
        'end';
        'function tf=get_tf';
        '[isOctave,v_num]=ver_info;';
        'if isOctave';
        '        % Octave normally comes with a compiler out of the box, but for some';
        '        % methods of installation an additional package may be required.';
        '        tf=~isempty(try_file_compile);';
        '    elseif v_num>=706 % ifversion(''>='',''R2008a'')';
        '        % Just try to compile a MWE. Getting the configuration is very slow.';
        '        [cc,TryNormalCheck]=try_file_compile;';
        '        if TryNormalCheck';
        '            % Something strange happened, so try the normal check anyway.';
        '            try cc=mex.getCompilerConfigurations;catch,cc=[];end';
        '        end';
        '        tf=~isempty(cc);';
        '    else';
        '        if ispc,ext=''.bat'';else,ext=''.sh'';end';
        '        tf=exist(fullfile(prefdir,[''mexopts'' ext]),''file'');';
        'end';
        'end';
        'function [isOctave,v_num]=ver_info';
        '% This is a compact and feature-poor equivalent of ifversion.';
        '% To save space this can be used as an alternative.';
        '% Example: R2018a is 9.4, so v_num will be 904.';
        'isOctave=exist(''OCTAVE_VERSION'', ''builtin'');';
        'v_num=version;';
        'ii=strfind(v_num,''.'');if numel(ii)~=1,v_num(ii(2):end)='''';ii=ii(1);end';
        'v_num=[str2double(v_num(1:(ii-1))) str2double(v_num((ii+1):end))];';
        'v_num=v_num(1)+v_num(2)/100;v_num=round(100*v_num);';
        'end';
        'function [cc,TryNormalCheck]=try_file_compile';
        'TryNormalCheck=false;';
        'try';
        '    [p,n]=fileparts(tempname);e=''.c'';';
        '    n=n(regexp(n,''[a-zA-Z0-9_]'')); % Keep only valid characters.';
        '    n=[''test_fun__'' n(1:min(15,end))];';
        '    fid=fopen(fullfile(p,[n e]),''w'');';
        '    fprintf(fid,''%s\n'',...';
        '        ''#include "mex.h"'',...';
        '        ''void mexFunction(int nlhs, mxArray *plhs[],'',...';
        '        ''  int nrhs, const mxArray *prhs[]) {'',...';
        '        ''    plhs[0]=mxCreateString("compiler works");'',...';
        '        ''    return;'',...';
        '        ''}'');';
        '    fclose(fid);';
        'catch';
        '    % If there is a write error in the temp dir, something is wrong.';
        '    % Just try the normal check.';
        '    cc=[];TryNormalCheck=true;return';
        'end';
        'try';
        '    current=cd(p);';
        'catch';
        '    % If the cd fails, something is wrong here. Just try the normal check.';
        '    cc=[];TryNormalCheck=true;return';
        'end';
        'try';
        '    mex([n e]);';
        '    cc=feval(str2func(n));';
        '    clear(n); % Clear to remove file lock.';
        '    cd(current);';
        'catch';
        '    % Either the mex or the feval failed. That means we can safely assume no';
        '    % working compiler is present. The normal check should not be required.';
        '    cd(current);';
        '    cc=[];TryNormalCheck=false;return';
        'end';
        'end';
        ''};
    fid=fopen(fn,'wt');fprintf(fid,'%s\n',txt{:});fclose(fid);
    fun_handle=str2func(fun_handle);

    [tf_,ME_]=feval(fun_handle);
end
tf=tf_;ME=ME_;
end
function str=convert_from_codepage(str,inverted)
%Convert from the Windows-1252 codepage.
persistent or ta
if isempty(or)
    %This list is complete for all characters (up to 0xFFFF) that can be encoded with ANSI.
    CPwin2UTF8=[338 140;339 156;352 138;353 154;376 159;381 142;382 158;402 131;710 136;732 152;
        8211 150;8212 151;8216 145;8217 146;8218 130;8220 147;8221 148;8222 132;8224 134;8225 135;
        8226 149;8230 133;8240 137;8249 139;8250 155;8364 128;8482 153];
    or=CPwin2UTF8(:,2);ta=CPwin2UTF8(:,1);
end
if nargin>1 && inverted
    origin=ta;target=or;
else
    origin=or;target=ta;
end
str=uint32(str);
for m=1:numel(origin)
    str=PatternReplace(str,origin(m),target(m));
end
end
function [p,ME]=CreatePathFolder__CheckMexCompilerExistence_persistent
%Try creating a folder in either the tempdir or a persistent folder and try adding it to the path
%(if it is not already in there). If the folder is not writable, the current folder will be used.
try
    ME=[];
    p=fullfile(GetWritableFolder,'FileExchange','CheckMexCompilerExistence');
    if isempty(strfind([path ';'],[p ';'])) %#ok<STREMP>
        %This means f is not on the path.
        if ~exist(p,'dir'),mkdir(p);end
        addpath(p,'-end');
    end
catch
    ME=struct('identifier','HJW:CheckMexCompilerExistence:PathFolderFail',...
        'message','Creating a folder on the path to store compiled mex files failed.');
end
end
function p=CreatePathFolder__sqlite3(print_to)
%Try creating a folder in either the tempdir or a persistent folder and try adding it to the path
%(if it is not already in there). If the folder is not writable, the current folder will be used.
try
    p=fullfile(GetWritableFolder,'FileExchange','sqlite3');
    if isempty(strfind([path ';'],[p ';'])) %#ok<STREMP>
        %This means f is not on the path.
        if ~exist(p,'dir'),mkdir(p);end
        addpath(p,'-end');
    end
catch
    error_(print_to,'HJW:sqlite3:PathFolderFail',...
        'Creating a folder on the path to store compiled mex files failed.')
end
end
function error_(options,varargin)
%Print an error to the command window, a file and/or the String property of an object.
%The error will first be written to the file and object before being actually thrown.
%
%Apart from controlling the way an error is written, you can also run a specific function. The
%'fcn' field of the options must be a struct (scalar or array) with two fields: 'h' with a function
%handle, and 'data' with arbitrary data passed as third input. These functions will be run with
%'error' as first input. The second input is a struct with identifier, message, and stack as
%fields. This function will be run with feval (meaning the function handles can be replaced with
%inline functions or anonymous functions).
%
%The intention is to allow replacement of every error(___) call with error_(options,___).
%
% NB: the error trace that is written to a file or object may differ from the trace displayed by
% calling the builtin error function. This was only observed when evaluating code sections.
%
%options.boolean.con: if true throw error with rethrow()
%options.fid:         file identifier for fprintf (array input will be indexed)
%options.boolean.fid: if true print error to file
%options.obj:         handle to object with String property (array input will be indexed)
%options.boolean.obj: if true print error to object (options.obj)
%options.fcn          struct (array input will be indexed)
%options.fcn.h:       handle of function to be run
%options.fcn.data:    data passed as third input to function to be run (optional)
%options.boolean.fnc: if true the function(s) will be run
%
%syntax:
%  error_(options,msg)
%  error_(options,msg,A1,...,An)
%  error_(options,id,msg)
%  error_(options,id,msg,A1,...,An)
%  error_(options,ME)               %equivalent to rethrow(ME)
%
%examples options struct:
%  % Write to a log file:
%  opts=struct;opts.fid=fopen('log.txt','wt');
%  % Display to a status window and bypass the command window:
%  opts=struct;opts.boolean.con=false;opts.obj=uicontrol_object_handle;
%  % Write to 2 log files:
%  opts=struct;opts.fid=[fopen('log2.txt','wt') fopen('log.txt','wt')];

persistent this_fun
if isempty(this_fun),this_fun=func2str(@error_);end

%Parse options struct, allowing an empty input to revert to default.
if isempty(options),options=validate_print_to__options(struct);end
options                   =parse_warning_error_redirect_options(  options  );
[id,msg,stack,trace,no_op]=parse_warning_error_redirect_inputs( varargin{:});
if no_op,return,end
ME=struct('identifier',id,'message',msg,'stack',stack);

%Print to object.
if options.boolean.obj
    msg_=msg;while msg_(end)==10,msg_(end)='';end%Crop trailing newline.
    if any(msg_==10)  % Parse to cellstr and prepend 'Error: '.
        msg_=char2cellstr(['Error: ' msg_]);
    else              % Only prepend 'Error: '.
        msg_=['Error: ' msg_];
    end
    for OBJ=options.obj(:).'
        try set(OBJ,'String',msg_);catch,end
    end
end

%Print to file.
if options.boolean.fid
    for FID=options.fid(:).'
        try fprintf(FID,'Error: %s\n%s',msg,trace);catch,end
    end
end

%Run function.
if options.boolean.fcn
    if ismember(this_fun,{stack.name})
        %To prevent an infinite loop, trigger an error.
        error('prevent recursion')
    end
    for FCN=options.fcn(:).'
        if isfield(FCN,'data')
            try feval(FCN.h,'error',ME,FCN.data);catch,end
        else
            try feval(FCN.h,'error',ME);catch,end
        end
    end
end

%Actually throw the error.
rethrow(ME)
end
function [valid,filename]=filename_is_valid(filename)
%Check if the file name and path are valid (non-empty char or scalar string).
valid=true;
persistent forbidden_names
if isempty(forbidden_names)
    forbidden_names={'CON','PRN','AUX','NUL','COM1','COM2','COM3','COM4','COM5','COM6','COM7',...
        'COM8','COM9','LPT1','LPT2','LPT3','LPT4','LPT5','LPT6','LPT7','LPT8','LPT9'};
    %ismember is case sensitive, so replicate to include lower case as well
    forbidden_names=[forbidden_names lower(forbidden_names)];
end
if isa(filename,'string') && numel(filename)==1
    %Convert a scalar string to a char array.
    filename=char(filename);
end
if ~isa(filename,'char') || numel(filename)==0
    valid=false;return
else
    %File name is indeed a char. Check if the folder exists (if specified), and do a check if there
    %are characters that can't exist in a normal file name. The method used here is not fool-proof,
    %but should cover most use cases and operating systems.
    [fullpath,fn,ext]=fileparts(filename);
    if ~isempty(fullpath) && ~exist(fullpath,'dir')
        valid=false;return
    end
    fn=[fn,ext];
    
    if any(ismember([char(0:31) '<>:"/\|?*'],fn)) || ...
            any(ismember(forbidden_names,fn)) || ...
            any(fn(end)=='. ')
        valid=false;return
    end
end
end
function varargout=findND(X,varargin)
%Find non-zero elements in ND-arrays. Replicates all behavior from find.
%
% The syntax is equivalent to the built-in find, but extended to multi-dimensional input.
%
% [...] = findND(X,K) returns at most the first K indices. K must be a positive scalar of any type.
%
% [...] = findND(X,K,side) returns either the first K or the last K indices. The input side  must
% be a char, either 'first' or 'last'. The default behavior is 'first'.
%
% [I1,I2,I3,...,In] = findND(X,...) returns indices along all the dimensions of X.
%
% [I1,I2,I3,...,In,V] = findND(X,...) returns indices along all the dimensions of X, and
% additionally returns a vector containing the values.
%
% Note for Matlab 6.5:
% The syntax with more than one input is present in the online doc for R14 (Matlab 7.0), so this
% might be the latest release without support for this syntax.
%
%  _______________________________________________________________________
% | Compatibility | Windows 10  | Ubuntu 20.04 LTS | MacOS 10.15 Catalina |
% |---------------|-------------|------------------|----------------------|
% | ML R2020a     |  works      |  not tested      |  not tested          |
% | ML R2018a     |  works      |  works           |  not tested          |
% | ML R2015a     |  works      |  works           |  not tested          |
% | ML R2011a     |  works      |  works           |  not tested          |
% | ML 6.5 (R13)  |  works      |  not tested      |  not tested          |
% | Octave 5.2.0  |  works      |  works           |  not tested          |
% | Octave 4.4.1  |  works      |  not tested      |  works               |
% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
%
% Version: 1.2.1
% Date:    2020-07-06
% Author:  H.J. Wisselink
% Licence: CC by-nc-sa 4.0 ( https://creativecommons.org/licenses/by-nc-sa/4.0 )
% Email = 'h_j_wisselink*alumnus_utwente_nl';
% Real_email = regexprep(Email,{'*','_'},{'@','.'})

%Parse inputs
if ~(isnumeric(X) || islogical(X)) || numel(X)==0
    error('HJW:findND:FirstInput',...
        'Expected first input (X) to be a non-empty numeric or logical array.')
end
switch nargin
    case 1 %[...] = findND(X);
        side='first';
        K=inf;
    case 2 %[...] = findND(X,K);
        side='first';
        K=varargin{1};
        if ~(isnumeric(K) || islogical(K)) || numel(K)~=1 || any(K<0)
            error('HJW:findND:SecondInput',...
                'Expected second input (K) to be a positive numeric or logical scalar.')
        end
    case 3 %[...] = FIND(X,K,'first');
        K=varargin{1};
        if ~(isnumeric(K) || islogical(K)) || numel(K)~=1 || any(K<0)
            error('HJW:findND:SecondInput',...
                'Expected second input (K) to be a positive numeric or logical scalar.')
        end
        side=varargin{2};
        if ~isa(side,'char') || ~( strcmpi(side,'first') || strcmpi(side,'last'))
            error('HJW:findND:ThirdInput','Third input must be either ''first'' or ''last''.')
        end
        side=lower(side);
    otherwise
        error('HJW:findND:InputNumber','Incorrect number of inputs.')
end

%parse outputs
nDims=length(size(X));
%allowed outputs: 0, 1, nDims, nDims+1
if nargout>1 && nargout<nDims
    error('HJW:findND:Output','Incorrect number of output arguments.')
end

persistent OldSyntax
if isempty(OldSyntax)
    OldSyntax=ifversion('<',7,'Octave','<',3);
end

varargout=cell(nargout,1);
if OldSyntax
    %The find(X,k,side) syntax was introduced between 6.5 and 7.
    if nargout>nDims
        [ind,col_index_equal_to_one,val]=find(X(:));%#ok no tilde pre-R2009b
        %X(:) converts X to a column vector. Treating X(:) as a matrix forces val to be the actual
        %value, instead of the column index.
        if length(ind)>K
            if strcmp(side,'first') %Select first K outputs.
                ind=ind(1:K);
                val=val(1:K);
            else                    %Select last K outputs.
                ind=ind((end-K+1):end);
                val=val((end-K+1):end);
            end
        end
        [varargout{1:(end-1)}] = ind2sub(size(X),ind);
        varargout{end}=val;
    else
        ind=find(X);
        if length(ind)>K
            if strcmp(side,'first')
                %Select first K outputs.
                ind=ind(1:K);
            else
                %Select last K outputs.
                ind=ind((end-K):end);
            end
        end
        [varargout{:}] = ind2sub(size(X),ind);
    end
else
    if nargout>nDims
        [ind,col_index_equal_to_one,val]=find(X(:),K,side);%#ok<ASGLU>
        %X(:) converts X to a column vector. Treating X(:) as a matrix forces val to be the actual
        %value, instead of the column index.
        [varargout{1:(end-1)}] = ind2sub(size(X),ind);
        varargout{end}=val;
    else
        ind=find(X,K,side);
        [varargout{:}] = ind2sub(size(X),ind);
    end
end
end
function fun_handle=get_sqlite3_mex(print_to)
% This will attempt to return a handle to a compiled mex interface.
% It will try by doing any of the following:
% - Re-add a persistent folder with the mex to the path.
% - Copy the pre-compiled binary to the persistent folder.
% - Download the pre-compiled binary.
% - Copy the source files and compile a mex binary.
% - Download the source files and compile a mex binary.

persistent persistent_fun_handle
if isempty(persistent_fun_handle)
    %In some release-runtime combinations addpath has a permanent effect, in others it doesn't. By
    %putting this code in this block, we are trying to keep these queries to a minimum.
    p=CreatePathFolder__sqlite3(print_to);

    % Octave mex files must be compiled for each machine, but Matlab mex files are a bit more
    % portable. Therefore, the pre-compiled Matlab mex files can be downloaded directly, but for
    % Octave we should download the source files and compile it.
    DoCompile=false;
    switch mexext
        case 'mex' % Octave
            % The mexname function will generate a file name unique to the platform and version of
            % Octave. This has the benefit of being able to use the same folder for multiple
            % compiled mex files.
            DoCompile=true;
            funname='sqlite3';
            [mexfilename,funname]=mexname(funname);
        case 'dll'
            % Compiled with Borland C++ Builder version 6.0
            funname='sqlite3_MATLAB_v06_05_PCWIN';   % R13
            mexfilename=[funname '.' mexext];
        case 'mexw32'
            % Compiled with Open Watcom 2
            funname='sqlite3_MATLAB_v07_01_PCWIN';   % R14SP3
            mexfilename=[funname '.' mexext];
        case 'mexa64'
            % Compiled with gcc version 9.3.0-17ubuntu1~20.04
            funname='sqlite3_MATLAB_v07_14_GLNXA64'; % R2012a
            mexfilename=[funname '.' mexext];
        case 'mexw64'
            % Compiled with Microsoft SDK 7.1
            funname='sqlite3_MATLAB_v08_02_PCWIN64'; % R2013b
            mexfilename=[funname '.' mexext];
        case 'mexmaci64'
            % Compiled with Xcode 11
            funname='sqlite3_MATLAB_v09_03_MACI64';  % R2017b
            mexfilename=[funname '.' mexext];
        otherwise
            ID='HJW:sqlite3:PlatformNotFound';
            warning_(print_to,ID,...
                ['The platform was not recognized.',char(10),...
                'A compilation of the sqlite3 mex file will be attempted.',char(10),...
                'Use warning(''off'',''%s'') to disable this warning.'],ID) %#ok<CHARTEN>
            DoCompile=true;
            funname='sqlite3';
            [mexfilename,funname]=mexname(funname);
    end
    mexfilename_full=fullfile(p,mexfilename);
    
    % Store the settings for the WBM calls in a persistent.
    get_sqlite3_mex___WBM_settings(struct(...
        'flag','id','print_to_fid',print_to.fid,'print_to_con',false));
    
    % An updater can be implemented here. The function can be called once and cleared (to prevent
    % the file being in use during a deletion attempt).
    %   code_version=feval(str2func(funname));clear(funname);
    % Now the file can be deleted and replaced by a newer version.
    if 3~=exist(mexfilename_full,'file')
        % Some files are included when the sqlite3 function is downloaded. These files should be
        % copied from the subdirectory to the folder that was just added to the path. If the files
        % aren't there, that means they should be downloaded. Source files only need to be
        % downloaded if a compile is required.
        base=fileparts(mfilename('fullpath'));
        if ~DoCompile
            subdir=fullfile(base,'compiled');
            if exist(fullfile(subdir,mexfilename),'file')
                copyfile(fullfile(subdir,mexfilename),mexfilename_full);
            end
            subdir=fullfile(base,'bin','sqlite3','compiled');
            if exist(fullfile(subdir,mexfilename),'file')
                copyfile(fullfile(subdir,mexfilename),mexfilename_full);
            end
            if 3~=exist(mexfilename_full,'file')
                baseurl='https://hjwisselink.nl/FEXsubmissiondata/68298-sqlite3/compiled/';
                fn=WBM(mexfilename_full,[baseurl mexfilename],get_sqlite3_mex___WBM_settings);
                if isempty(fn),delete(fn);end
            end
            if 3~=exist(mexfilename_full,'file')
                DoCompile=true;
            end
        end

        % Attempt the compilation.
        isOctave=exist('OCTAVE_VERSION', 'builtin');
        if DoCompile
            [CompilerExists,ME]=CheckMexCompilerExistence;
            if ~CompilerExists,error_(print_to,ME),end
            %- download required files to a subfolder inside the temp folder
            %- cd to that folder
            %- compile mex file
            %- cd back
            %- copy mex file to the persistent folder
            %- delete temp subfolder
            try
                tmpfolder=[tmpname('compile_sqlite3_mex_tmp_folder') filesep];
                mkdir(tmpfolder);
                get_sqlite3_mex___download_sources(tmpfolder,print_to)
                current=cd(tmpfolder);
                get_sqlite3_mex___compile(mexfilename)
                copyfile(mexfilename,mexfilename_full);
                cd(current);
            catch
            end
            try
                if isOctave,oldval=confirm_recursive_rmdir(false);end
                rmdir(tmpfolder,'s');
                if isOctave,confirm_recursive_rmdir(oldval);end
            catch
            end
        end
    end
    if 3~=exist(mexfilename_full,'file')
        d=fullfile(fileparts(mfilename('fullpath')),get_sqlite3_mex___src_dir);
        error_(print_to,'HJW:sqlite3:CompilationFailed',...
            ['The download/compilation of the sqlite3 mex file failed.',char(10),...
            'If you want to manually download the compiled mex file, please put it here:',...
            char(10),p,char(10),...
            'If you have a working compiler, put ther source files here for auto-compile:',...
            char(10),d]) %#ok<CHARTEN>
    end
    persistent_fun_handle=str2func(funname);
end
fun_handle=persistent_fun_handle;
end
function p=get_sqlite3_mex___src_dir
p='source_v3.0.0';
end
function out=get_sqlite3_mex___WBM_settings(in)
persistent s
if nargin==1,s=in;end
out=s;
end
function get_sqlite3_mex___download_sources(target_folder)
% This function will download the C source files required for the mex interface.
%
% The C source files come from 3 different sources:
%
% sqlite3.c and sqlite3.h:
% http://web.archive.org/web/202108id_/
% https://www.sqlite.org/2021/sqlite-amalgamation-3360000.zip
%
% Originals for sqlite3_interface.c, structlist.c, and structlist.h (edited to make them conform to
% the stricter standards of older compilers and remove a mexPrintf call about param binding):
% http://web.archive.org/web/202108id_/
% https://raw.githubusercontent.com/rmartinjak/mex-sqlite3/master/sqlite3.c
% https://raw.githubusercontent.com/rmartinjak/mex-sqlite3/master/structlist.c
% https://raw.githubusercontent.com/rmartinjak/mex-sqlite3/master/structlist.h
%
% CharConversions.c:
% This converter file was written to deal with UTF-16 char encoding. It isn't perfect, but it will
% get most of the way there.

files={...
    'sqlite3.c',...
    'sqlite3.h',...
    'sqlite3_interface_strict.c',...
    'structlist_strict.c',...
    'structlist_strict.h',...
    'CharConversions.c'};
baseurl=['https://hjwisselink.nl/FEXsubmissiondata/68298-sqlite3/' get_sqlite3_mex___src_dir '/'];
base_local=fileparts(mfilename('fullpath'));
for n=1:numel(files)
    fn=fullfile(target_folder,files{n});
    if exist(fn,'file'),continue,end

    local=fullfile(base_local,get_sqlite3_mex___src_dir,files{n});
    if exist(local,'file'),copyfile(local,fn);continue,end
    local=fullfile(base_local,'bin','sqlite3',get_sqlite3_mex___src_dir,files{n});
    if exist(local,'file'),copyfile(local,fn);continue,end

    % The local file exists in neither location, so download it.
    WBM(fn,[baseurl files{n}],get_sqlite3_mex___WBM_settings);
end
end
function get_sqlite3_mex___compile(mexfilename)
v=version;ind=strfind(v,'.');v(ind(2):end)='';v=str2double(v);
isOctave=exist('OCTAVE_VERSION', 'builtin');

mexfilename=strrep(mexfilename,'.dll','');
c_files={'sqlite3_interface_strict.c','sqlite3.c','structlist_strict.c'};

if ispc
    UTF16flag={'-DSQLITE_OMIT_UTF16 #1 '};
    if isOctave || v<9
        %suppress a warning on new compilers and avoid an error on old ones
        UTF16flag{1}=strrep(UTF16flag{1},' #1 ','#1');
    end
else
    UTF16flag={'-DSQLITE_OMIT_UTF16=1'};
end

if ispc
    mex('-output',mexfilename,c_files{:},UTF16flag{:})
else
    %both Mac and Linux
    mex('-output',mexfilename,c_files{:},'-ldl')
end

if isOctave
    files=strrep(c_files,'.c','.o');
    for n=1:numel(files)
        if exist(files{n},'file'),delete(files{n});end
    end
end
end
function [str,stack]=get_trace(skip_layers,stack)
if nargin==0,skip_layers=1;end
if nargin<2, stack=dbstack;end
stack(1:skip_layers)=[];

%Parse the ML6.5 style of dbstack (the name field includes full file location).
if ~isfield(stack,'file')
    for n=1:numel(stack)
        tmp=stack(n).name;
        if strcmp(tmp(end),')')
            %Internal function.
            ind=strfind(tmp,'(');
            name=tmp( (ind(end)+1):(end-1) );
            file=tmp(1:(ind(end)-2));
        else
            file=tmp;
            [ignore,name]=fileparts(tmp); %#ok<ASGLU>
        end
        [ignore,stack(n).file]=fileparts(file); %#ok<ASGLU>
        stack(n).name=name;
    end
end

%Parse Octave style of dbstack (the file field includes full file location).
persistent IsOctave,if isempty(IsOctave),IsOctave=exist('OCTAVE_VERSION','builtin');end
if IsOctave
    for n=1:numel(stack)
        [ignore,stack(n).file]=fileparts(stack(n).file); %#ok<ASGLU>
    end
end

%Create the char array with a (potentially) modified stack.
s=stack;
c1='>';
str=cell(1,numel(s)-1);
for n=1:numel(s)
    [ignore_path,s(n).file,ignore_ext]=fileparts(s(n).file); %#ok<ASGLU>
    if n==numel(s),s(n).file='';end
    if strcmp(s(n).file,s(n).name),s(n).file='';end
    if ~isempty(s(n).file),s(n).file=[s(n).file '>'];end
    str{n}=sprintf('%c In %s%s (line %d)\n',c1,s(n).file,s(n).name,s(n).line);
    c1=' ';
end
str=horzcat(str{:});
end
function atomTime=getUTC(override)
%Returns the UTC time in the Matlab datenum format
%
%example syntax:
% disp(datestr(getUTC))
%
% There are several methods implemented in this function:
% - An implementation that requires a C mex function.
%   This method requires write access to a folder and a working C compiler. The compilation result
%   will be stored to a subdirectory of a folder similar to the AddOn path, or the tempdir, or the
%   current folder. Write permission is tested in that order.
% - An implementation using https://www.utctime.net/utc-timestamp.
%   The NIST has a server that returns the time, but it currently blocks API access.
%   This method requires internet access.
% - The local time and timezone offset can be determined with the wmic command or the get-date
%   Powershell function (Windows) or the date command (Linux and Mac).
%   To speed up the usage of this method, you can cache the difference with now() in a persistent
%   variable, that way you avoid the need for a slow system call.
% 
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 2.1.0                                                         |%
%|  Date:    2022-01-22                                                    |%
%|  Author:  H.J. Wisselink                                                |%
%|  Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 ) |%
%|  Email = 'h_j_wisselink*alumnus_utwente_nl';                            |%
%|  Real_email = regexprep(Email,{'*','_'},{'@','.'})                      |%
%|                                                                         |%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%
% Tested on several versions of Matlab (ML 6.5 and onward) and Octave (4.4.1 and onward), and on
% multiple operating systems (Windows/Ubuntu/MacOS). For the full test matrix, see the HTML doc.
% Compatibility considerations:
% - Some older releases don't support the web implementation.
% - The normal system call hangs on ML7.1 on XP. Since ML6.5 works fine on Windows 10, it seems
%   reasonable to assume that the OS is the cause of the hang. For XP (and older) there is an
%   alternative strategy in place, but this has a higher likelyhood to fail.
% - Similarly, as wmic has been deprecated, a Powershell alternative should be used on newer
%   versions of Windows. The need for this is automatically detected.

if nargin==0
    %Normal flow: first try the cmd method, then the C method, then the web method.
    UTC_epoch_seconds=getUTC_cmd;
    if isempty(UTC_epoch_seconds)
        UTC_epoch_seconds=getUTC_c;
    end
    if isempty(UTC_epoch_seconds)
        UTC_epoch_seconds=getUTC_web;
    end
    if isempty(UTC_epoch_seconds)
        error('HJW:getUTC:TimeReadFailed',...
            ['All methods of retrieving the UTC timestamp failed.\nEnsure you ',...
            'have write access to the current folder and check your internet connection.'])
    end
else
    %Override for debug/test, this will not throw an error on fail.
    if override==1
        UTC_epoch_seconds=getUTC_c(false);
    elseif override==2
        UTC_epoch_seconds=getUTC_web;
    elseif override==3
        UTC_epoch_seconds=getUTC_cmd;
    else
        error('non-implemented override')
    end
end
UTC_offset=UTC_epoch_seconds/(24*60*60);
atomTime=UTC_offset+datenum(1970,1,1);
end
function UTC_epoch_seconds=getUTC_c(allow_rethrow)
%Use a C implementation, which requires write permission in a folder.
if nargin==0,allow_rethrow=true;end
persistent utc_time_c tempdir_f funname utc_time_fun_handle Compile_attempts_remaining mexfilename
if isempty(utc_time_c)
    %Try creating a folder in either the tempdir or a persistent folder and adding it to the path
    %(if it is not already in there). If the folder is not writable, the current folder will be
    %used.
    %In some release-runtime combinations addpath has a permanent effect, in others it doesn't. By
    %putting this code in this block, we are trying to keep these queries to a minimum.
    tempdir_f=fullfile(GetWritableFolder,'FileExchange','getUTC');
    try
        if isempty(strfind([path ';'],[tempdir_f ';'])) %#ok<STREMP>
            %This means f is not on the path.
            if ~exist(tempdir_f,'dir'),mkdir(tempdir_f);end
            addpath(tempdir_f,'-end');
        end
    catch
    end
    
    funname='utc_time';
    [mexfilename,funname]=mexname(funname);
    try utc_time_fun_handle=str2func(funname);catch,end %The try-catch is required for Octave.
    
    %Only allow a few compilation attempts, so this function doesn't cause a lot of disk I/O if
    %there is no working compiler.
    Compile_attempts_remaining=5;
    
    %prepare to write this to a file and compile
    utc_time_c={'#include "mex.h"';
        '#include "time.h"';
        '';
        '/* Abraham Cohn,  3/17/2005 */';
        '/* Philips Medical Systems */';
        '';
        'void mexFunction(int nlhs, mxArray *plhs[], int nrhs,';
        '                 const mxArray *prhs[])';
        '{';
        '  time_t utc;';
        '  ';
        '  if (nlhs > 1) {';
        '    mexErrMsgTxt("Too many output arguments");';
        '  }';
        '  ';
        '  /* Here is a nice ref: www.cplusplus.com/ref/ctime/time.html */';
        '  time(&utc);';
        '  /* mexPrintf("UTC time in local zone: %s",ctime(&utc)); */';
        '  /* mexPrintf("UTC time in GMT: %s",asctime(gmtime(&utc))); */';
        '  ';
        '  /* Create matrix for the return argument. */';
        '  plhs[0] = mxCreateDoubleScalar((double)utc);';
        '   ';
        '}'};
    %(the original had mxCreateScalarDouble)
end

try
    UTC_epoch_seconds=feval(utc_time_fun_handle);
catch
    if exist(mexfilename,'file')
        if allow_rethrow
            ME=lasterror; %#ok<LERR>
            rethrow(ME);
        else
            UTC_epoch_seconds=[];return
        end
    end
    
    %Check if there is compiler in the first place.
    if ~CheckMexCompilerExistence
        Compile_attempts_remaining=0;
        UTC_epoch_seconds=[];return
    end
    
    %Try building missing C file.
    Compile_attempts_remaining=Compile_attempts_remaining-1;
    if Compile_attempts_remaining<0 %Don't endlessly try to compile.
        UTC_epoch_seconds=[];return
    end
    
    if TestFolderWritePermission(tempdir_f)
        f=tempdir_f; %Use the folder in the tempdir to store the mex.
    else
        f=pwd; %Revert to current folder.
    end
    current_folder=cd(f);
    try
        if ~exist(fullfile(f,[funname '.c']),'file')
            fid=fopen(fullfile(f,[funname '.c']),'w');
            for line=1:numel(utc_time_c)
                fprintf(fid,'%s\n',utc_time_c{line});
            end
            fclose(fid);
        end
        
        try
            mex([funname '.c']);
        catch
        end
        %Perform cleanup.
        for ext={'c','o'}
            file=fullfile(f,[funname '.' ext{1}]);if exist(file,'file'),delete(file),end
        end
    catch
    end
    cd(current_folder);
    if exist(mexfilename,'file')
        utc_time_fun_handle=str2func(funname); %Refresh the handle.
        UTC_epoch_seconds=getUTC_c(allow_rethrow); %Use recursion to catch errors.
    else
        %The compiling of the mex function failed
        UTC_epoch_seconds=[];
    end
end
end
function UTC_epoch_seconds=getUTC_cmd
%Use a command line implementation.
%Should return an empty array instead of an error if it fails.
try
    call_type=getUTC_cmd_call_type;
catch
    warning('determination of call type failed')
    UTC_epoch_seconds=[];return
end

try
   switch call_type
       case 'Unix'
           UTC_epoch_seconds=getUTC_cmd_Unix;
       case 'WMIC_sys'
           UTC_epoch_seconds=getUTC_cmd_wmic_sys;
       case 'WMIC_bat'
           UTC_epoch_seconds=getUTC_cmd_wmic_bat;
       case 'PS_get_date'
           UTC_epoch_seconds=getUTC_cmd_PS_get_date;
   end
catch
    UTC_epoch_seconds=[];
end
end
function call_type=getUTC_cmd_call_type
persistent call_type_
if isempty(call_type_)
    if ~ispc
        call_type_='Unix';
    else
        if WinVer<=5
            %The normal system call hangs on ML7.1 on XP. Since ML6.5 works fine on Windows 10, I'm
            %making the assumption that the OS is the cause of the hang.
            call_type_='WMIC_bat';
        elseif WinVer>=10
            %The cmd WMIC interface has been deprecated and seems to have been removed from Windows
            %10 21H1. On older W10 versions we can still use WMIC. Since every version of Windows
            %will either have the WMIC or PowerShell, a PowerShell call seems like a good fallback.
            [wmic_not_available,ignore]=system('wmic /?'); %#ok<ASGLU>
            if wmic_not_available
                call_type_='PS_get_date';
            else
                call_type_='WMIC_sys';
            end
        else
            call_type_='WMIC_sys';
        end
    end
end
call_type=call_type_;
end
function UTC_epoch_seconds=getUTC_cmd_PS_get_date
%Use Powershell to get the UTC time.
[s,str]=system(['powershell $a=get-date;',...
    '$a.ToUniversalTime().ToString(''yyyyMMddHHmmss'')']); %#ok<ASGLU>
str(str<48 | str>57)='';%Remove trailing newline by keeping only digits.
date=mat2cell(str,1,[4 2 2,2 2 2]);date=num2cell(str2double(date));
UTC_epoch_seconds=(datenum(date{:})-datenum(1970,1,1))*24*60*60;
end
function UTC_epoch_seconds=getUTC_cmd_Unix
[status,str]=system('date +%s'); %#ok<ASGLU>
UTC_epoch_seconds=str2double(str);
end
function UTC_epoch_seconds=getUTC_cmd_wmic_bat
%If a normal system call would hang, use a temporary bat file.
pausetime=1;
fn1=[tempname,'.bat'];fn2=[tempname,'.txt'];
fid=fopen(fn1,'w');
%This command returns YYYYMMDDHHMMSS.milliseconds+UTC_Offset_in_minutes
fprintf(fid,'%%systemroot%%\\system32\\wbem\\wmic os get LocalDateTime /value > "%s"\r\nexit',fn2);
fclose(fid);
system(['start /min "" cmd /c "' fn1 '"']);
then=now;%Store the current time to adjust for the pausetime and the function time.
pause(pausetime) % Wait for the system call to finish.
str=fileread(fn2);
try delete(fn1);catch,end,try delete(fn2);catch,end % Delete temp files.
UTC_epoch_seconds=getUTC_cmd_wmic_parse_str(str)+(now-then)*24*60*60;
end
function UTC_epoch_seconds=getUTC_cmd_wmic_sys
%Use a direct system call (instead of a temporary bat file).
[status,str]=system('wmic os get LocalDateTime /value'); %#ok<ASGLU>
%This command returns YYYYMMDDHHMMSS.milliseconds+UTC_Offset_in_minutes
UTC_epoch_seconds=getUTC_cmd_wmic_parse_str(str);
end
function UTC_epoch_seconds=getUTC_cmd_wmic_parse_str(str)
str=str(str>=43 & str<=57); % Strip irrelevant parts.
date=mat2cell(str(1:21),1,[4 2 2,2 2 2+1+6]);date=str2double(date);
date(5)=date(5)-str2double(str(22:end)); % Add offset.
date=num2cell(date);
UTC_epoch_seconds=(datenum(date{:})-datenum(1970,1,1))*24*60*60;
end
function UTC_epoch_seconds=getUTC_web
%Read the timestamp from a web server.
%For some reason this fails for old Matlab releases (at least ML6.5-R2013b).

persistent UseWebread%Only search the path once per session.
if isempty(UseWebread)
    try UseWebread=~isempty(which(func2str(@webread)));catch,UseWebread=false;end
end

%Skip this function if there is no internet connection (3 timeouts will take a lot of time).
if ~isnetavl,UTC_epoch_seconds=[];return,end

for tries=1:3
    try
        if UseWebread
            data=webread('http://www.utctime.net/utc-timestamp');
        else
            %This probably only works for Octave.
            data=urlread('http://www.utctime.net/utc-timestamp'); %#ok<URLRD>
        end
        break
    catch
    end
end
try
    data(data==' ')='';
    pat='vartimestamp=';
    ind1=strfind(data,pat)+numel(pat);
    ind2=strfind(data,';')-1;
    ind2(ind2<ind1)=[];
    UTC_epoch_seconds=str2double(data(ind1:ind2(1)));
catch
    UTC_epoch_seconds=[];
end
end
function [f,status]=GetWritableFolder(varargin)
%Return a folder with write permission
%
% If the output folder doesn't already exist, this function will attempt to create it. This
% function should provide a reliable and repeatable location to write files.
%
% Syntax:
%   f=GetWritableFolder
%   [f,status]=GetWritableFolder
%   [__]=GetWritableFolder(Name,Value)
%   [__]=GetWritableFolder(options)
%
% Input/output arguments:
% f:
%   Char array with the full path to the writable folder. This does not contain a trailing filesep.
% status:
%   A scalar double ranging from 0 to 3. 0 denotes a failure to find a folder, 1 means the folder
%   is in a folder close to the AddOn folder, 2 that it is a folder in the tempdir, 3 mean that the
%   returned path is a folder in the current directory.
% options:
%   A struct with Name,Value parameters. Missing parameters are filled with the defaults listed
%   below. Using incomplete parameter names or incorrect capitalization is allowed, as long as
%   there is a unique match.
%
% Name,Value parameters:
%   ForceStatus:
%      Retrieve the path corresponding to the status value. Using 0 allows an automatic
%      determination of the location (default=0;).
%    ErrorOnNotFound:
%      Throw an error when failing to find a writeable folder (default=true;).
%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 1.0.0                                                         |%
%|  Date:    2021-02-19                                                    |%
%|  Author:  H.J. Wisselink                                                |%
%|  Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 ) |%
%|  Email = 'h_j_wisselink*alumnus_utwente_nl';                            |%
%|  Real_email = regexprep(Email,{'*','_'},{'@','.'})                      |%
%|                                                                         |%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%
% Tested on several versions of Matlab (ML 6.5 and onward) and Octave (4.4.1 and onward), and on
% multiple operating systems (Windows/Ubuntu/MacOS). For the full test matrix, see the HTML doc.
% Compatibility considerations:
% - The path returned with status=1 is mostly the same as the addonpath for most releases. Although
%   it is not correct for all release/OS combinations, it should still work. If you have a managed
%   account, this might result in strange behavior.

[success,options,ME]=GetWritableFolder_parse_inputs(varargin{:});
if ~success
    rethrow(ME)
else
    [ForceStatus,ErrorOnNotFound,root_folder_list]=deal(options.ForceStatus,...
        options.ErrorOnNotFound,options.root_folder_list);
end
root_folder_list{end}=pwd;% Set this default here to avoid storing it in a persistent.
if ForceStatus
    status=ForceStatus;f=fullfile(root_folder_list{status},'PersistentFolder');
    try if ~exist(f,'dir'),mkdir(f);end,catch,end
    return
end

% Option 1: use a folder similar to the AddOn Manager.
status=1;f=root_folder_list{status};
try if ~exist(f,'dir'),mkdir(f);end,catch,end
if ~TestFolderWritePermission(f)
    % If the Add-On path is not writable, return the tempdir. It will not be persistent, but it
    % will be writable.
    status=2;f=root_folder_list{status};
    try if ~exist(f,'dir'),mkdir(f);end,catch,end
    if ~TestFolderWritePermission(f)
        % The tempdir should always be writable, but if for some reason it isn't: return the pwd.
        status=3;f=root_folder_list{status};
    end
end

% Add 'PersistentFolder' to whichever path was determined above.
f=fullfile(f,'PersistentFolder');
try if ~exist(f,'dir'),mkdir(f);end,catch,end

if ~TestFolderWritePermission(f)
    % Apparently even the pwd isn't writable, so we will either return an error, or a fail state.
    if ErrorOnNotFound
        error('HJW:GetWritableFolder:NoWritableFolder',...
            'This function was unable to find a folder with write permissions.')
    else
        status=0;f='';
    end
end
end
function [success,options,ME]=GetWritableFolder_parse_inputs(varargin)
%Parse the inputs of the GetWritableFolder function
% This function returns a success flag, the parsed options, and an ME struct.
% As input, the options should either be entered as a struct or as Name,Value pairs. Missing fields
% are filled from the default.

% Pre-assign outputs.
success=false;
ME=struct('identifier','','message','');

persistent default
if isempty(default)
    %Set defaults for options.
    default.ForceStatus=false;
    default.ErrorOnNotFound=false;
    default.root_folder_list={...
        GetPseudoAddonpath;
        fullfile(tempdir,'MATLAB');
        ''};% Overwrite this last element with pwd when called.
end

if nargin==2
    options=default;
    success=true;
    return
end

% Actually parse the Name,Value pairs (or the struct).
[options,replaced]=parse_NameValue(default,varargin{:});

% Test the optional inputs.
for k=1:numel(replaced)
    curr_option=replaced{k};
    item=options.(curr_option);
    ME.identifier=['HJW:GetWritableFolder:incorrect_input_opt_' lower(curr_option)];
    switch curr_option
        case 'ForceStatus'
            try
                if ~isa(default.root_folder_list{item},'char')
                    % This ensures an error for item=[true false true]; as well.
                    error('the indexing must have failed, trigger error')
                end
            catch
                ME.message=sprintf('Invalid input: expected a scalar integer between 1 and %d.',...
                    numel(default.root_folder_list));
                return
            end
        case 'ErrorOnNotFound'
            [passed,options.ErrorOnNotFound]=test_if_scalar_logical(item);
            if ~passed
                ME.message='ErrorOnNotFound should be either true or false.';
                return
            end
        otherwise
            ME.message=sprintf('Name,Value pair not recognized: %s.',curr_option);
            ME.identifier='HJW:GetWritableFolder:incorrect_input_NameValue';
            return
    end
end
success=true;ME=[];
end
function f=GetPseudoAddonpath
% This is mostly the same as the addonpath. Technically this is not correct for all release/OS
% combinations and the code below should be used:
%     addonpath='';
%     try s = Settings;addonpath=get(s.matlab.addons,'InstallationFolder');end %#ok<TRYNC>
%     try s = Settings;addonpath=get(s.matlab.apps,'AppsInstallFolder');end %#ok<TRYNC>
%     try s = settings;addonpath=s.matlab.addons.InstallationFolder.ActiveValue;end %#ok<TRYNC>
%
% However, this returns an inconsistent output:
%     R2011a          <pref doesn't exist>
%     R2015a Ubuntu  $HOME/Documents/MATLAB/Apps
%            Windows %HOMEPATH%\MATLAB\Apps
%     R2018a Ubuntu  $HOME/Documents/MATLAB/Add-Ons
%            Windows %HOMEPATH%\MATLAB\Add-Ons
%     R2020a Windows %APPDATA%\MathWorks\MATLAB Add-Ons
%
% To make the target folder consistent, only one of these options is chosen.
if ispc
    [ignore,appdata]=system('echo %APPDATA%');appdata(appdata<14)=''; %#ok<ASGLU> (remove LF/CRLF)
    f=fullfile(appdata,'MathWorks','MATLAB Add-Ons');
else
    [ignore,home_dir]=system('echo $HOME');home_dir(home_dir<14)=''; %#ok<ASGLU> (remove LF/CRLF)
    f=fullfile(home_dir,'Documents','MATLAB','Add-Ons');
end
end
function tf=ifversion(test,Rxxxxab,Oct_flag,Oct_test,Oct_ver)
%Determine if the current version satisfies a version restriction
%
% To keep the function fast, no input checking is done. This function returns a NaN if a release
% name is used that is not in the dictionary.
%
% Syntax:
%   tf=ifversion(test,Rxxxxab)
%   tf=ifversion(test,Rxxxxab,'Octave',test_for_Octave,v_Octave)
%
% Input/output arguments:
% tf: 
%   If the current version satisfies the test this returns true. This works similar to verLessThan.
% Rxxxxab:
%   A char array containing a release description (e.g. 'R13', 'R14SP2' or 'R2019a') or the numeric
%   version (e.g. 6.5, 7, or 9.6).
% test:
%   A char array containing a logical test. The interpretation of this is equivalent to
%   eval([current test Rxxxxab]). For examples, see below.
%
% Examples:
% ifversion('>=','R2009a') returns true when run on R2009a or later
% ifversion('<','R2016a') returns true when run on R2015b or older
% ifversion('==','R2018a') returns true only when run on R2018a
% ifversion('==',9.9) returns true only when run on R2020b
% ifversion('<',0,'Octave','>',0) returns true only on Octave
% ifversion('<',0,'Octave','>=',6) returns true only on Octave 6 and higher
%
% The conversion is based on a manual list and therefore needs to be updated manually, so it might
% not be complete. Although it should be possible to load the list from Wikipedia, this is not
% implemented.
%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 1.1.0                                                         |%
%|  Date:    2021-09-26                                                    |%
%|  Author:  H.J. Wisselink                                                |%
%|  Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 ) |%
%|  Email = 'h_j_wisselink*alumnus_utwente_nl';                            |%
%|  Real_email = regexprep(Email,{'*','_'},{'@','.'})                      |%
%|                                                                         |%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%
% Tested on several versions of Matlab (ML 6.5 and onward) and Octave (4.4.1 and onward), and on
% multiple operating systems (Windows/Ubuntu/MacOS). For the full test matrix, see the HTML doc.
% Compatibility considerations:
% - This is expected to work on all releases.

%The decimal of the version numbers are padded with a 0 to make sure v7.10 is larger than v7.9.
%This does mean that any numeric version input needs to be adapted. multiply by 100 and round to
%remove the potential for float rounding errors.
%Store in persistent for fast recall (don't use getpref, as that is slower than generating the
%variables and makes updating this function harder).
persistent  v_num v_dict octave
if isempty(v_num)
    %test if Octave is used instead of Matlab
    octave=exist('OCTAVE_VERSION', 'builtin');
    
    %get current version number
    v_num=version;
    ii=strfind(v_num,'.');if numel(ii)~=1,v_num(ii(2):end)='';ii=ii(1);end
    v_num=[str2double(v_num(1:(ii-1))) str2double(v_num((ii+1):end))];
    v_num=v_num(1)+v_num(2)/100;v_num=round(100*v_num);
    
    %get dictionary to use for ismember
    v_dict={...
        'R13' 605;'R13SP1' 605;'R13SP2' 605;'R14' 700;'R14SP1' 700;'R14SP2' 700;
        'R14SP3' 701;'R2006a' 702;'R2006b' 703;'R2007a' 704;'R2007b' 705;
        'R2008a' 706;'R2008b' 707;'R2009a' 708;'R2009b' 709;'R2010a' 710;
        'R2010b' 711;'R2011a' 712;'R2011b' 713;'R2012a' 714;'R2012b' 800;
        'R2013a' 801;'R2013b' 802;'R2014a' 803;'R2014b' 804;'R2015a' 805;
        'R2015b' 806;'R2016a' 900;'R2016b' 901;'R2017a' 902;'R2017b' 903;
        'R2018a' 904;'R2018b' 905;'R2019a' 906;'R2019b' 907;'R2020a' 908;
        'R2020b' 909;'R2021a' 910;'R2021b' 911};
end

if octave
    if nargin==2
        warning('HJW:ifversion:NoOctaveTest',...
            ['No version test for Octave was provided.',char(10),...
            'This function might return an unexpected outcome.']) %#ok<CHARTEN>
        if isnumeric(Rxxxxab)
            v=0.1*Rxxxxab+0.9*fix(Rxxxxab);v=round(100*v);
        else
            L=ismember(v_dict(:,1),Rxxxxab);
            if sum(L)~=1
                warning('HJW:ifversion:NotInDict',...
                    'The requested version is not in the hard-coded list.')
                tf=NaN;return
            else
                v=v_dict{L,2};
            end
        end
    elseif nargin==4
        % Undocumented shorthand syntax: skip the 'Octave' argument.
        [test,v]=deal(Oct_flag,Oct_test);
        % Convert 4.1 to 401.
        v=0.1*v+0.9*fix(v);v=round(100*v);
    else
        [test,v]=deal(Oct_test,Oct_ver);
        % Convert 4.1 to 401.
        v=0.1*v+0.9*fix(v);v=round(100*v);
    end
else
    % Convert R notation to numeric and convert 9.1 to 901.
    if isnumeric(Rxxxxab)
        v=0.1*Rxxxxab+0.9*fix(Rxxxxab);v=round(100*v);
    else
        L=ismember(v_dict(:,1),Rxxxxab);
        if sum(L)~=1
            warning('HJW:ifversion:NotInDict',...
                'The requested version is not in the hard-coded list.')
            tf=NaN;return
        else
            v=v_dict{L,2};
        end
    end
end
switch test
    case '==', tf= v_num == v;
    case '<' , tf= v_num <  v;
    case '<=', tf= v_num <= v;
    case '>' , tf= v_num >  v;
    case '>=', tf= v_num >= v;
end
end
function [connected,timing]=isnetavl(use_HTML_test_only)
%Check for an internet connection by pinging Google
%
% Syntax:
%   [connected,timing]=isnetavl(use_HTML_test_only)
%
% Input/output arguments:
% connected:
%   A logical denoting the connectivity. Note that this may return unexpected results if Matlab has
%   separate proxy settings and/or if your DNS is having issues.
% timing:
%   This contains the ping time as a double and defaults to 0 if there is no connection.
%   Note that this value will be larger than the true ping time if the HTML fallback is used.
% use_HTML_test_only:
%   If the input is convertible to true the HTML method will be used, which has the benefit of
%   testing if Matlab (or Octave) is able to connect to the internet. This might be relevant if
%   there are proxy settings or firewall rules specific to Matlab/Octave.
%   Note that the ping value will be larger than the true value if the HTML fallback is used.
%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 2.0.0                                                         |%
%|  Date:    2021-05-19                                                    |%
%|  Author:  H.J. Wisselink                                                |%
%|  Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 ) |%
%|  Email = 'h_j_wisselink*alumnus_utwente_nl';                            |%
%|  Real_email = regexprep(Email,{'*','_'},{'@','.'})                      |%
%|                                                                         |%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%
% Tested on several versions of Matlab (ML 6.5 and onward) and Octave (4.4.1 and onward), and on
% multiple operating systems (Windows/Ubuntu/MacOS). For the full test matrix, see the HTML doc.
% Compatibility considerations:
% - This is expected to work on all releases.
% - If you have proxy/firewall settings specific to Matlab/Octave, make sure to use the HTML method
%   by providing an input.

if nargin==0,                                             tf=false;
else,[tf1,tf2]=test_if_scalar_logical(use_HTML_test_only);tf=tf1&&tf2;
end
if tf
    %(the timing is not reliable)
    [connected,timing]=isnetavl___ping_via_html;
    return
end

tf=isnetavl__ICMP_is_blocked;
if isempty(tf)
    %Unable to determine if ping is allowed, the connection must be down.
    connected=0;
    timing=0;
else
    if tf
        %Ping is not allowed.
        %(the timing is not reliable)
        [connected,timing]=isnetavl___ping_via_html;
    else
        %Ping is allowed.
        [connected,timing]=isnetavl___ping_via_system;
    end
end
end
function [connected,timing]=isnetavl___ping_via_html
%Ping is blocked by some organizations. As an alternative, the google.com page can be loaded as a
%normal HTML, which should work as well, although it is slower. This also means the ping timing is
%no longer reliable.
persistent UseWebread
if isempty(UseWebread)
    try no_webread=isempty(which(func2str(@webread)));catch,no_webread=true;end
	UseWebread=~no_webread;
end
try
    then=now;
    if UseWebread
        str=webread('http://google.com'); %#ok<NASGU>
    else
        str=urlread('http://google.com'); %#ok<NASGU,URLRD>
    end
    connected=1;
    timing=(now-then)*24*3600*1000;
catch
    connected=0;
    timing=0;
end
end
function [connected,timing]=isnetavl___ping_via_system
if ispc
    try
        %                                   8.8.4.4 will also work
        [ignore_output,b]=system('ping -n 1 8.8.8.8');%#ok<ASGLU> ~
        stats=b(strfind(b,' = ')+3);
        stats=stats(1:3);%[sent received lost]
        if ~strcmp(stats,'110')
            error('trigger error')
        else
            %This branch will error for 'destination host unreachable'.
            connected=1;
            %This assumes there is only one place with '=[digits]ms' in the response, but this code
            %is not language-specific.
            [ind1,ind2]=regexp(b,' [0-9]+ms');
            timing=b((ind1(1)+1):(ind2(1)-2));
            timing=str2double(timing);
        end
    catch
        connected=0;
        timing=0;
    end
elseif isunix
    try
        %                                   8.8.4.4 will also work
        [ignore_output,b]=system('ping -c 1 8.8.8.8');%#ok<ASGLU> ~
        ind=regexp(b,', [01] ');
        if b(ind+2)~='1'
            %This branch includes 'destination host unreachable' errors.
            error('trigger error')
        else
            connected=1;
            %This assumes the first place with '=[digits] ms' in the response contains the ping
            %timing. This code is not language-specific.
            [ind1,ind2]=regexp(b,'=[0-9.]+ ms');
            timing=b((ind1(1)+1):(ind2(1)-2));
            timing=str2double(timing);
        end
    catch
        connected=0;
        timing=0;
    end
else
    error('How did you even get Matlab to work?')
end
end
function [tf,connected,timing]=isnetavl__ICMP_is_blocked
%Check if ICMP 0/8/11 is blocked.
%
%tf is empty if both methods fail.

persistent output
if ~isempty(output)
    tf=output;return
end

%First check if ping works.
[connected,timing]=isnetavl___ping_via_system;
if connected
    %Ping worked and there is an internet connection.
    output=false;
    tf=false;
    return
end

%There are two options: no internet connection, or ping is blocked.
[connected,timing]=isnetavl___ping_via_html;
if connected
    %There is an internet connection, therefore ping must be blocked.
    output=true;
    tf=true;
    return
end

%Both methods failed, internet is down. Leave the value of tf (and the persistent variable) set to
%empty so it is tried next time.
tf=[];
end
function [object,ME]=JSON(str,varargin)
%This interprets char array as JSON and returns an object the same size and shape as the builtin
%
% For very small files this function may be faster than the builtin, but for large files this
% function may be much slower. See the performance section in the HTML doc for detailed timing
% information.
%
% Syntax:
%   object=JSON(str)
%   object=JSON(___,options)
%   object=JSON(___,Name,Value)
%   [object,ME]=JSON(___)
%
% object:
%   This contains the parsed object. This should closely match the output of the Matlab builtin
%   jsondecode (see below for details).
% ME:
%   Errors during parsing will be caught. If a second output argument is specified, the error
%   will not be rethrown, but the corresponding MException object is returned instead.
% str:
%   The JSON string to be parsed. This should be a char vector or a string/cellstr.
% options:
%   A struct with Name,Value parameters. Missing parameters are filled with the default.
%   Note that the parameters are not validated.
%
% Name,Value parameters:
%   EnforceValidNumber:
%      With this boolean you can turn off the check if a number conforms to the JSON
%      specifications. This will cause str2double to determine the validity. No error will be
%      thrown in case of NaN output. [default=true;]
%   ThrowErrorForInvalid:
%      If this is false, no error will be throw if parsing fails. Instead, an empty array is
%      returned. [default=nargout<2;]
%   MaxRecursionDepth:
%      This function is a recursive function. Under some rare conditions, Matlab/Octave might crash
%      when the maximum recursion depth is reached, instead of throwing an error. This parameter
%      allows you to stay on the safe side.
%      The value can be set to inf to effectively remove the limit and only rely on the builtin
%      safeguards. [default=101-numel(dbstack);]
%
% This was implemented using
% https://www.ecma-international.org/wp-content/uploads/ECMA-404_2nd_edition_december_2017.pdf
% The failing and passing cases were validated using the test cases from a JSON test suite on
% GitHub (http://github.com/nst/JSONTestSuite), containing over 300 cases of possibly ambiguous
% syntax. Because the standard is not explicit for every situation, there are also test cases left
% to the implementation.
%
% Implementation details for compatibility with the jsondecode Matlab function:
% - If possible an array of arrays is treated as a row-major matrix.
%   * all data types should match
%   * all elements must be vectors of the same size
% - The null literal is an empty double, unless it is an element of an array, in which case it is
%   parsed as a NaN.
% - The name of an object is converted to a valid field name with a function similar to genvarname
%   or matlab.lang.makeValidName. This means characters after some whitespace characters are
%   converted to uppercase, all whitespace is removed, invalid characters are replaced by an
%   underscore, and an x is used as a prefix if the resulting name is empty or starts with a number
%   or underscore. In case of duplicate object names, an underscore and a counter is added. A
%   char(0) will cause the name to be cropped.
% - An empty array ('[]') is encoded as an empty double (val=[]).
% - An empty array array of arrays ('[[]]') is encoded as an empty cell (val={[]}).
%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 1.0.1                                                         |%
%|  Date:    2022-01-22                                                    |%
%|  Author:  H.J. Wisselink                                                |%
%|  Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 ) |%
%|  Email = 'h_j_wisselink*alumnus_utwente_nl';                            |%
%|  Real_email = regexprep(Email,{'*','_'},{'@','.'})                      |%
%|                                                                         |%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%
% Tested on several versions of Matlab (ML 6.5 and onward) and Octave (4.4.1 and onward), and on
% multiple operating systems (Windows/Ubuntu/MacOS). For the full test matrix, see the HTML doc.
% Compatibility considerations:
% - The recursion depth is limited to 100. This will affect all combinations of nesting of arrays
%   and objects. Without this, the recursion limit may be reached, causing Matlab/Octave to crash.
%   Matlab/Octave should prevent this on their own, but this is in place in case that protection
%   fails. The value can be set to inf to effectively remove the limit and only rely on the builtin
%   safeguards.

opts=struct;
opts.EnforceValidNumber=true;
opts.ThrowErrorForInvalid=nargout<2;
opts.MaxRecursionDepth=101-numel(dbstack);
opts=parse_NameValue(opts,varargin{:});
w=struct;w.w=warning('off','REGEXP:multibyteCharacters');[w.msg,w.ID]=lastwarn;
try ME=[];
    notation=StrToNotation(str,opts);
    object=ParseValue(notation);
catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
end
warning(w.w);lastwarn(w.msg,w.ID);%Reset warning state.
if ~isempty(ME)%An error has occurred.
    if opts.ThrowErrorForInvalid,rethrow(ME)
    else,                        object=[]; end
end
end
function notation=StrToNotation(str,opts)
%Treat the notation as a class, with the following properties (i.e. struct fields):
%  str : the actual notation input
%  s_tokens : a copy of str with only the structural tokens
%  braces : an array encoding the pairs of braces (0 for non-brace characters)
%  depth : approximate current recursion depth
%  opts : a struct with options
%
% Any whitespace allowed by the specification (i.e. [9 10 13 32]) will be removed outside of the
% strings to facilitate easier parsing.

if isa(str,'string'),str=cellstr(str);end
if iscellstr(str),str=sprintf('%s\n',str{:});end
if ~isa(str,'char') || numel(str)~=length(str)
    throw_error('The input should be a char vector or a string/cellstr.','Input')
end
str=reshape(str,1,[]);%Ensure str is a row vector.

persistent args ws legacy
if isempty(args)
    ws=[9 10 13 32];%this is not equal to \s
    args={['[' char(ws) ']*([\[{\]}:,])[' char(ws) ']*'],'$1','tokenize'};
    
    
    %The 'tokenize' option became the default in R14 (v7).
    v=version;v(min(find(v=='.')):end)='';pre_v7=str2double(v)<7; %#ok<MXFND>
    isOctave=exist('OCTAVE_VERSION','builtin')~=0;
    if ~pre_v7 && ~isOctave
        args(end)=[];
    end
    
    %Check if the regex crops special characters like char(0).
    legacy= 3==numel(regexprep(['123' char([0 10])],args{:})) ;
end

%This is true for the text including the double quotes.
txt=FindLiteralTextPositions(str);

%Remove whitespace that will not be parsed.
s_tokens=str;s_tokens(txt)='_';
if ~legacy
    s_tokens=regexprep(s_tokens,args{:});
else
    %Since the regex crops characters like char(0), we need to apply the regex to a surrogate.
    tmp=s_tokens;tmp(:)='n';                             %Mark all non-relevant.
    tmp(ismember(double(s_tokens),ws))='w';              %Mark all whitespace
    tmp(ismember(double(s_tokens),double('[{}]:,')))='t';%Mark all tokens
    L=zeros(1,1+numel(tmp));%Extend size by 1 to make cumsum work.
    [s,e]=regexp(tmp,'w+t');
    if ~isempty(s),L(s  )=1;L(e  )=-1;end
    [s,e]=regexp(tmp,'tw+');
    if ~isempty(s),L(s+1)=1;L(e+1)=-1;end
    L=logical(cumsum(L));L(end)=[];
    s_tokens(L)='';
end
while numel(s_tokens)>0 && any(s_tokens(end)==ws),s_tokens(end)='';end%Deal with lone values.
while numel(s_tokens)>0 && any(s_tokens(1)==ws),s_tokens(1)='';end%Deal with lone values.

%Apply whitespace removal to str as well.
str2=s_tokens;str2(s_tokens=='_')=str(txt);str=str2;

%Leave only structural tokens.
s_tokens(~ismember(double(s_tokens),double('[{}]:,')))='_';

[braces,ArrayOfArrays]=PairBraces(s_tokens);

notation.str=str;
notation.s_tokens=s_tokens;
notation.braces=braces;
notation.ArrayOfArrays=ArrayOfArrays;
notation.opts=opts;
notation.depth=0;
end
function txt=FindLiteralTextPositions(str)
%Characters after an odd number and before an even number double quote are part of a literal text.
%Since this is only true if there are no escaped double quotes, we need to mask those.
str=strrep(str,'\"','__');%Mask escaped double quotes.
L=str=='"';%Find all double quotes.
x=cumsum(L);%Count prior double quotes.
txt=L | mod(x,2)~=0;%Find the location of literal text.
end
function [braces,ArrayOfArrays]=PairBraces(s_tokens)
%Pair the braces and brackets and determine which positions are part of an array of arrays.

braces=zeros(size(s_tokens));
L=s_tokens=='{';
braces(L)=1:sum(L);
try
    L2=find(L);
    for ind=find(s_tokens=='}')
        %Equivalent to match=find(L(1:ind),1,'last').
        ind2=find(L2<ind);ind2=ind2(end);match=L2(ind2);L2(ind2)=[];
        braces(ind)=-braces(match);
        L(match)=false;
    end
    if any(L),error('trigger'),end
catch
    msg='Unmatched braces found.';
    id='PairBraces';
    throw_error(msg,id)
end
L=s_tokens=='[';
braces(L)=max(braces)+(1:sum(L));
try
    L2=find(L);
    for ind=find(s_tokens==']')
        %Equivalent to match=find(L(1:ind),1,'last').
        ind2=find(L2<ind);ind2=ind2(end);match=L2(ind2);L2(ind2)=[];
        braces(ind)=-braces(match);
        L(match)=false;
    end
    if any(L),error('trigger'),end
catch
    msg='Unmatched braces found.';
    id='PairBraces';
    throw_error(msg,id)
end

ArrayOfArrays=false(1,numel(s_tokens));
%Keep only the brackets and create a reduced form of the indices of the braces.
x=s_tokens(ismember(double(s_tokens),double('[{}]')));
br=braces(braces~=0);
%Find the starting positions of the array of arrays (i.e. the inner-most '[[').
[s,e]=regexp(x,'(\[)+\[');e=e-1;
for n=1:numel(s)
    %Mark each array of arrays with a counter.
    ArrayOfArrays(braces==br(e(n)))=true;
end
end
function notation=index(notation,indices)
%If the second input is a scalar, this is interpreted as indices:end.
if numel(indices)==1,indices(2)=numel(notation.str);end
indices=indices(1):indices(2);
notation.str=notation.str(indices);
notation.s_tokens=notation.s_tokens(indices);
notation.braces=notation.braces(indices);
notation.ArrayOfArrays=notation.ArrayOfArrays(indices);
end
function throw_error(msg,id)
error(['HJW:JSON:' id],msg)
end
function val=ParseValue(notation)
persistent number
if isempty(number)
    number=num2cell('-0123456789');
end
if numel(notation.str)==0,notation.str=' ';end%This will trigger an error later.
%Limit recursion depth to avoid crashes.
notation.depth=notation.depth+1;
if notation.depth>notation.opts.MaxRecursionDepth
    throw_error('Recursion limit reached, exiting to avoid crashes.','Recursion')
end

switch notation.str(1)
    case '{'
        val=ParseObject(notation);
    case '['
        val=ParseArray(notation);
    case number
        val=ParseNumber(notation.str,notation.opts.EnforceValidNumber);
    case '"'
        val=ParseString(notation.str);
        %Avoid 1x0 chars.
        if numel(val)==0,val='';end
    case 't'
        if ~strcmp(notation.str,'true')
            msg='Unexpected literal, expected ''true''.';
            id='Literal';
            throw_error(msg,id)
        end
        val=true;
    case 'f'
        if ~strcmp(notation.str,'false')
            msg='Unexpected literal, expected ''false''.';
            id='Literal';
            throw_error(msg,id)
        end
        val=false;
    case 'n'
        if ~strcmp(notation.str,'null')
            msg='Unexpected literal, expected ''null''.';
            id='Literal';
            throw_error(msg,id)
        end
        val=[];
    otherwise
        msg='Unexpected character, expected a brace, bracket, number, string, or literal.';
        id='Literal';
        throw_error(msg,id)
end
end
function val=ParseObject(notation)
ind=find(notation.braces==-notation.braces(1));%Find the matching closing brace.
if numel(ind)~=1 || ind~=numel(notation.str) || ~strcmp(notation.str([1 end]),'{}')
    msg='Unexpected end of object.';
    id='Object';
    throw_error(msg,id)
end

val=struct;
if ind==2,return,end%Empty object: '{}'.

%Select the part of the notation between the braces.
to_parse=index(notation,[2 numel(notation.str)-1]);

%Split over the commas that are not inside braces.
c_ind=find(cumsum(to_parse.braces)==0 & to_parse.s_tokens==',');
c_ind=[0 c_ind numel(to_parse.str)+1];
c_ind=[c_ind(1:(end-1))+1;c_ind(2:end)-1];

%Split each pair in the string part and a value.
brace_content=cell(size(c_ind));
for n=1:size(c_ind,2)
    pair=index(to_parse,c_ind(:,n));
    L=pair.s_tokens==':';
    if ~any(L),throw_error('No colon found in object definition.','Object'),end
    ind=find(L);ind=ind(1);
    try
        brace_content{1,n}=ParseString(pair.str(1:(ind-1)));
    catch
        throw_error('Invalid key in object definition.','Object')
    end
    try
        brace_content{2,n}=ParseValue(index(pair,ind+1));
    catch
        throw_error('Invalid value in object definition.','Object')
    end
    if isa(brace_content{2,n},'cell')
        %Wrap in a scalar cell to avoid creating a struct array later.
        brace_content{2,n}=brace_content(2,n);
    end
end

%Determine the fieldnames.
persistent RE dyn_expr
if isempty(RE)
    ws=char([9 10 11 12 13 32]);
    RE=['[' ws ']+([^' ws '])'];
    
    %Dynamically check if the dynamic expression replacement is available. This is possible since
    %R2006a (v7.2), and has not been implemented yet in Octave 6.2.0.
    dyn_expr=strcmp('fooBar',regexprep('foo bar',RE,'${upper($1)}'));
end
for n=1:size(brace_content,2)
    fn=brace_content{1,n};
    if dyn_expr
        fn=regexprep(fn,RE,'${upper($1)}');%Convert characters after whitespace to upper case.
    else
        [s,e]=regexp(fn,RE);
        if ~isempty(s)
            fn(e)=upper(fn(e));
            L=zeros(size(fn));L(s)=1;L(e)=-1;L=logical(cumsum(L));
            fn(L)='';
        end
    end
    fn=regexprep(fn,'\s','');%Remove all remaining whitespace.
    x=find(double(fn)==0);if ~isempty(x),fn(x(1):end)='';end%Gobble null characters.
    fn=regexprep(fn,'[^0-9a-zA-Z_]','_');%Replace remaining invalid characters.
    if isempty(fn)||any(fn(1)=='_0123456789')
        fn=['x' fn]; %#ok<AGROW>
    end
    fn_=fn;counter=0;
    while ismember(fn,brace_content(1,1:(n-1)))
        counter=counter+1;fn=sprintf('%s_%d',fn_,counter);
    end
    brace_content{1,n}=fn;
end

%Store to struct
val=struct(brace_content{:});
end
function str=ParseString(str)
persistent dict_ dict_2 symbol_length
if isempty(dict_)
    symbol_length=zeros(1,double('u'));symbol_length(double('"\/bfnrtu'))=[ones(1,8) 5]+1;
    dict_={...
        '"','"';...
        '\','\';...
        '/','/';...
        'b',char(8);...
        'f',char(12);...
        'n',char(10);...
        'r',char(13);...
        't',char(9)}; %#ok<CHARTEN>
    dict_2=double([cell2mat(dict_(:,1));'u']);
    for n=1:size(dict_,1),dict_{n,1}=['\' dict_{n,1}];end
end
if ~strcmp(str([1 end]),'""')
    msg='Unexpected end of string.';
    id='StringDelim';
    throw_error(msg,id)
end
if any(double(str)<32)
    msg='Unescaped control character in string.';
    id='StringControlChar';
    throw_error(msg,id)
end
str=str(2:(end-1));%Remove outer double quotes.
ind=regexp(str,'\\.');
if any(ind)
    %Create a unique list of all replacements. To prevent double replacement, split the str in a
    %cell and replace each element separately.
    
    %Check if there are only valid escaped characters.
    if ~all(ismember(double(str(ind+1)),dict_2))
        msg='Unexpected escaped character.';
        id='StringEscape';
        throw_error(msg,id)
    end
    
    %Find the true indices (this will properly deal with "\\\u000D").
    ind=regexp(str,'\\["\\/bfnrtu]');
    
    %Find all unicode replacements.
    ind2=strfind(str,'\u');
    if ~isempty(ind2)
        HasUnicode=true;
        try
            ind2=bsxfun_plus(ind2.',0:5);
            U=cellstr(unique(str(ind2),'rows'));
            for n=1:numel(U)
                U{n,2}=unicode_to_char(hex2dec(U{n,1}(3:end)));
            end
        catch
            msg='Unexpected escaped character.';
            id='StringEscape';
            throw_error(msg,id)
        end
    else
        HasUnicode=false;
        U=cell(0,2);
    end
    dict=[U;dict_];
    
    %Encode the length of each segment.
    len=symbol_length(str(ind+1));
    len=[0 len;diff([0,ind,1+numel(str)])-[1 len]];
    str_c=mat2cell(str,1,len(:));
    str_c=reshape(str_c,2,[]);
    
    %Check for \ in normal text, this is probably already caught.
    if any([str_c{2,:}]=='\')
        msg='Unexpected escaped character.';
        id='StringEscape';
        throw_error(msg,id)
    end
    
    %Keep track of processed elements so we don't do unnecessary work.
    done=true(size(str_c));
    done(1,:)=false;
    done(cellfun('prodofsize',str_c)==0)=true;
    for n=1:size(dict,1)
        L=~done & ismember(str_c,dict{n,1});
        str_c(L)=dict(n,2);
        done(L)=true;
        if all(done),break,end
    end
    if any(~done)
        msg='Unexpected escaped character.';
        id='StringEscape';
        throw_error(msg,id)
    end
    str=horzcat(str_c{:});
    if HasUnicode && CharIsUTF8
        %Fix the mis-encoding for surrogate pairs in UTF-16. This is only relevant on runtimes
        %where char is encoded as UTF-8. Currently, that means this only applies to Octave.
        [str,ignore_flag]=UTF8_to_unicode(str); %#ok<ASGLU>
        str=unicode_to_char(UTF16_to_unicode(str));
    end
end
end
function num=ParseNumber(str,EnforceValidNumber)
if ~EnforceValidNumber
    isValidJSON=true;
else
    %The complete regexp below doesn't work for some reason, so do it in two passes.
    %
    %       ['-?((0)|([1-9]+\d*))(\.\d+)?([eE]' '[\+-]' '?[0-9]+)?']
    expr{1}=['-?((0)|([1-9]+\d*))(\.\d+)?([eE]'  '\+'   '?[0-9]+)?'];
    expr{2}=['-?((0)|([1-9]+\d*))(\.\d+)?([eE]'    '-'  '?[0-9]+)?'];
    if true
        [s,e]=regexp(str,expr{1},'once');
        isValidJSON=~isempty(s) && s==1 && e==numel(str);
    end
    if ~isValidJSON
        [s,e]=regexp(str,expr{2},'once');
        isValidJSON=~isempty(s) && s==1 && e==numel(str);
    end
end
if ~isValidJSON
    msg='Invalid number format.';
    id='Number';
    throw_error(msg,id)
end
num=str2double(str);
end
function val=ParseArray(notation)
ind=find(notation.braces==-notation.braces(1));%Find the matching closing bracket.
if numel(ind)~=1 || ind~=numel(notation.str) || ~strcmp(notation.str([1 end]),'[]')
    msg='Unexpected end of array.';
    id='Array';
    throw_error(msg,id)
end

if strcmp(notation.str,'[]')
    %Empty array: '[]'.
    val=[];return
end
if strcmp(notation.str,'[[]]')
    %Empty array of arrays: '{[]}'.
    val={[]};return
end

%Avoid indexing (to select the part of the notation between the brackets).
br=notation.braces;br(1)=0;

%Split over the commas that are not inside brackets.
c_ind=find(cumsum(br)==0 & notation.s_tokens==',');
c_ind=[1 c_ind numel(notation.str)];
c_ind=[c_ind(1:(end-1))+1;c_ind(2:end)-1];
if any(diff(c_ind,1,1)<0)
    msg='Empty array element.';
    id='Array';
    throw_error(msg,id)
end
val=cell(size(c_ind,2),1);
for n=1:size(c_ind,2)
    val{n}=ParseValue(index(notation,c_ind(:,n)));
end

%An array of arrays should be treated as a row-major matrix.
%These are the requirements to parse to a matrix:
% - all data types should match
% - all elements must be vectors of the same size
%   (null being NaN instead of [] in the case of a numeric matrix)
if notation.ArrayOfArrays(1)
    try
        tmp=val;
        tmp(cellfun('isempty',tmp))={NaN};
        val=horzcat(tmp{:}).';
        return
    catch
        %Revert to cell vector.
    end
end
if ismember(class(val{n}),{'double','logical','struct'})
    tmp=val;
    if all(cellfun('isclass',val,'double'))
        tmp(cellfun('isempty',tmp))={NaN};
        if numel(unique(cellfun('prodofsize',tmp)))==1
            val=horzcat(tmp{:}).';
            return
        end
    elseif all(cellfun('isclass',val,'logical'))
        if numel(unique(cellfun('prodofsize',val)))==1
            val=horzcat(val{:}).';
            return
        end
    elseif all(cellfun('isclass',val,'struct'))
        %This will fail for dissimilar structs.
        try val=vertcat(val{:});catch,end
    else
        %Leave as cell vector.
    end
end
end
function [mex_filename,fun_name]=mexname(fun_name)
%Encode runtime version information in the function name.
% This can be useful if multiple versions of Matlab or Octave need to use the
% same folder to store compiled functions, while not being compatible.
%
% This function replaces a syntax like mex_filename=[fun_name '.' mexext].
%
% Syntax:
%   mex_filename=mexname(fun_name);
%   [mex_filename,updated_fun_name]=mexname(fun_name);
persistent append
if isempty(append)
    v=version;ind=[strfind(v,'.') numel(v)];
    v=sprintf('%02d.%02d',str2double({...
        v(1:(ind(1)-1)           ) ,...
        v(  (ind(1)+1):(ind(2)-1)) }));
    v=['v' strrep(v,'.','_')];
    if ~exist('OCTAVE_VERSION', 'builtin')
        runtime='MATLAB';
        type=computer;
    else
        runtime='OCTAVE';
        arch=computer;arch=arch(1:(min(strfind(arch,'-'))-1));
        if ispc
            if strcmp(arch,'x86_64')  ,type= 'win_64';
            elseif strcmp(arch,'i686'),type= 'win_i686';
            elseif strcmp(arch,'x86') ,type= 'win_x86';
            else                      ,type=['win_' arch];
            end
        elseif isunix && ~ismac %Essentially this is islinux
            if strcmp(arch,'i686')      ,type= 'lnx_i686';
            elseif strcmp(arch,'x86_64'),type= 'lnx_64';
            else                        ,type=['lnx_' arch];
            end
        elseif ismac
            if strcmp(arch,'x86_64'),type= 'mac_64';
            else                    ,type=['mac_' arch];
            end
        end
    end
    type=strrep(strrep(type,'.',''),'-','');
    append=cell(2,1);
    append{1}=['_' runtime '_' v '_' type];
    append{2}=[append{1} '.' mexext];
end

try %Test if fun_name is a valid name.
    if ~isvarname(fun_name),error('trigger catch block'),end
catch
    error('HJW:mexname:InvalidName',...
        'The provided input can''t be a function name')
end

mex_filename=[fun_name append{2}];
fun_name=[fun_name append{1}];
end
function [opts,replaced]=parse_NameValue(default,varargin)
%Match the Name,Value pairs to the default option, ignoring incomplete names and case.
%
%If this fails to find a match, this will throw an error with the offending name as the message.
%
%The second output is a cell containing the field names that have been set.

opts=default;replaced={};
if nargin==1,return,end

%Unwind an input struct to Name,Value pairs.
try
    struct_input=numel(varargin)==1 && isa(varargin{1},'struct');
    NameValue_input=mod(numel(varargin),2)==0 && all(...
        cellfun('isclass',varargin(1:2:end),'char'  ) | ...
        cellfun('isclass',varargin(1:2:end),'string')   );
    if ~( struct_input || NameValue_input )
        error('trigger')
    end
    if nargin==2
        Names=fieldnames(varargin{1});
        Values=struct2cell(varargin{1});
    else
        %Wrap in cellstr to account for strings (this also deals with the fun(Name=Value) syntax).
        Names=cellstr(varargin(1:2:end));
        Values=varargin(2:2:end);
    end
    if ~iscellstr(Names),error('trigger');end %#ok<ISCLSTR>
catch
    %If this block errors, that is either because a missing Value with the Name,Value syntax, or
    %because the struct input is not a struct, or because an attempt was made to mix the two
    %styles. In future versions of this functions an effort might be made to handle such cases.
    error('parse_NameValue:MixedOrBadSyntax',...
        'Optional inputs must be entered as Name,Value pairs or as struct.')
end

%Convert the real fieldnames to a char matrix.
d_Names1=fieldnames(default);d_Names2=lower(d_Names1);
len=cellfun('prodofsize',d_Names2);maxlen=max(len);
for n=find(len<maxlen).' %Pad with spaces where needed
    d_Names2{n}((end+1):maxlen)=' ';
end
d_Names2=vertcat(d_Names2{:});

%Attempt to match the names.
replaced=false(size(d_Names1));
for n=1:numel(Names)
    name=lower(Names{n});
    tmp=d_Names2(:,1:min(end,numel(name)));
    non_matching=numel(name)-sum(cumprod(double(tmp==repmat(name,size(tmp,1),1)),2),2);
    match_idx=find(non_matching==0);
    if numel(match_idx)~=1
        error('parse_NameValue:NonUniqueMatch',Names{n})
    end
    
    %Store the Value in the output struct.
    opts.(d_Names1{match_idx})=Values{n};
    
    replaced(match_idx)=true;
end
replaced=d_Names1(replaced);
end
function [id,msg,stack,trace,no_op]=parse_warning_error_redirect_inputs(varargin)
no_op=false;
if nargin==1
    %  error_(options,msg)
    %  error_(options,ME)
    if isa(varargin{1},'struct') || isa(varargin{1},'MException')
        ME=varargin{1};
        if numel(ME)==0
            no_op=true;
            [id,msg,stack,trace]=deal('');
            return
        end
        try
            stack=ME.stack;%Use the original call stack if possible.
            trace=get_trace(0,stack);
        catch
            [trace,stack]=get_trace(3);
        end
        id=ME.identifier;
        msg=ME.message;
        pat='Error using <a href="matlab:matlab.internal.language.introspective.errorDocCallback(';
        %This pattern may occur when using try error(id,msg),catch,ME=lasterror;end instead of
        %catching the MException with try error(id,msg),catch ME,end.
        %This behavior is not stable enough to robustly check for it, but it only occurs with
        %lasterror, so we can use that.
        if isa(ME,'struct') && numel(msg)>numel(pat) && strcmp(pat,msg(1:numel(pat)))
            %Strip the first line (which states 'error in function (line)', instead of only msg).
            msg(1:find(msg==10,1))='';
        end
    else
        [trace,stack]=get_trace(3);
        [id,msg]=deal('',varargin{1});
    end
else
    [trace,stack]=get_trace(3);
    if ~isempty(strfind(varargin{1},'%')) %The id can't contain a percent symbol.
        %  error_(options,msg,A1,...,An)
        id='';
        A1_An=varargin(2:end);
        msg=sprintf(varargin{1},A1_An{:});
    else
        %  error_(options,id,msg)
        %  error_(options,id,msg,A1,...,An)
        id=varargin{1};
        msg=varargin{2};
        if nargin>2
            A1_An=varargin(3:end);
            msg=sprintf(msg,A1_An{:});
        end
    end
end
end
function options=parse_warning_error_redirect_options(options)
%Fill the struct:
%options.boolean.con (this field is ignored in error_)
%options.boolean.fid
%options.boolean.obj
%options.boolean.fcn
if ~isfield(options,'boolean'),options.boolean=struct;end
if ~isfield(options.boolean,'con') || isempty(options.boolean.con)
    options.boolean.con=false;
end
if ~isfield(options.boolean,'fid') || isempty(options.boolean.fid)
    options.boolean.fid=isfield(options,'fid');
end
if ~isfield(options.boolean,'obj') || isempty(options.boolean.obj)
    options.boolean.obj=isfield(options,'obj');
end
if ~isfield(options.boolean,'fcn') || isempty(options.boolean.fcn)
    options.boolean.fcn=isfield(options,'fcn');
end
end
function out=PatternReplace(in,pattern,rep)
%Functionally equivalent to strrep, but extended to more data types.
out=in(:)';
if numel(pattern)==0
    L=false(size(in));
elseif numel(rep)>numel(pattern)
    error('not implemented (padding required)')
else
    L=true(size(in));
    for n=1:numel(pattern)
        k=find(in==pattern(n));
        k=k-n+1;k(k<1)=[];
        %Now k contains the indices of the beginning of each match.
        L2=false(size(L));L2(k)=true;
        L= L & L2;
        if ~any(L),break,end
    end
end
k=find(L);
if ~isempty(k)
    for n=1:numel(rep)
        out(k+n-1)=rep(n);
    end
    if numel(rep)==0,n=0;end
    if numel(pattern)>n
        k=k(:);%Enforce direction.
        remove=(n+1):numel(pattern);
        idx=bsxfun_plus(k,remove-1);
        out(idx(:))=[];
    end
end
end
function data=readfile(filename,varargin)
%Read a UTF-8 or ANSI (US-ASCII) file
%
% Syntax:
%   data=readfile(filename)
%   data=readfile(___,options)
%   data=readfile(___,Name,Value)
%
% Input/output arguments:
% data:
%   An n-by-1 cell (1 cell per line in the file, even empty lines).
% filename:
%   A char array with either relative or absolute path, or a URL.
% options:
%   A struct with Name,Value parameters. Missing parameters are filled with the defaults listed
%   below. Using incomplete parameter names or incorrect capitalization is allowed, as long as
%   there is a unique match.
%   Parameters related to warning/error redirection will be parsed first.
%
% Name,Value parameters:
%   err_on_ANSI:
%      If set to true, an error will be thrown when the input file is not recognized as UTF-8
%      encoded. This should normally not be an issue, as ANSI files can be read as well with this
%      function. [default=false;]
%   EmptyLineRule:
%      This contains a description of how empty lines should be handled. Lines that only contain
%      whitespace are considered empty as well, to conform to the behavior of readlines (this
%      therefore also depends on the Whitespace parameter). Valid values are 'read', 'skip',
%      'error', 'skipleading', and 'skiptrailing'.
%      The latter two are not available for readlines. Values can be entered as a scalar string or
%      as a char array. [default='read';]
%   WhitespaceRule:
%      This contains a description of how should leading and trailing whitespace be handled on each
%      line. Depending on the value of the Whitespace parameter this is equivalent to readlines.
%      Valid values are 'preserve', 'trim', 'trimleading', and 'trimtrailing'.
%      [default='preserve';]
%   LineEnding:
%      This parameter determines which characters are considered line ending characters. String
%      arrays and cell arrays of char vectors are parsed by sprintf, with each element being
%      considered a line break. String scalars and character vectors are treated as literal.
%      The default is {'\n','\r','\r\n'} meaning that \n\r is considered 2 line ends. This will not
%      be checked for any overlap and will be processed sequentially. The only is the default,
%      which will be sorted to {'\r\n','\n','\r'}. [default={'\n','\r','\r\n'};]
%   Whitespace:
%      This parameter determines which characters are treated as whitespace for the purposes of
%      EmptyLineRule and WhitespaceRule. This should be a char vector or a scalar string. Cell
%      arrays of char vectors are parsed by sprintf and concatenated. Note that the default for
%      readlines is sprintf(' \b\t'), but in this function this is expanded.
%      [default=[8 9 28:32 160 5760 8192:8202 8239 8287 12288];]
%   UseReadlinesDefaults:
%      Reproduce the default behavior of readlines as closely as possible. This includes
%      reproducing a bug which causes all characters that require 2 uint16 values to encode in
%      UTF-16 (everything outside the base multilingual plane, i.e. most emoji) to be converted to
%      char(26).
%      This will not convert the output to a string array. [default=false;]
%   print_to_con:
%      An attempt is made to also use this parameter for warnings or errors during input parsing.
%      A logical that controls whether warnings and other output will be printed to the command
%      window. Errors can't be turned off. [default=true;] if either print_to_fid, print_to_obj, or
%      print_to_fcn is specified then [default=false]
%   print_to_fid:
%      An attempt is made to also use this parameter for warnings or errors during input parsing.
%      The file identifier where console output will be printed. Errors and warnings will be
%      printed including the call stack. You can provide the fid for the command window (fid=1) to
%      print warnings as text. Errors will be printed to the specified file before being actually
%      thrown. [default=[];]
%      If print_to_fid, print_to_obj, and print_to_fcn are all empty, this will have the effect of
%      suppressing every output except errors. 
%      Array inputs are allowed.
%   print_to_obj:
%      An attempt is made to also use this parameter for warnings or errors during input parsing.
%      The handle to an object with a String property, e.g. an edit field in a GUI where console
%      output will be printed. Messages with newline characters (ignoring trailing newlines) will
%      be returned as a cell array. This includes warnings and errors, which will be printed
%      without the call stack. Errors will be written to the object before the error is actually
%      thrown. [default=[];]
%      If print_to_fid, print_to_obj, and print_to_fcn are all empty, this will have the effect of
%      suppressing every output except errors.
%      Array inputs are allowed.
%   print_to_fcn:
%      An attempt is made to also use this parameter for warnings or errors during input parsing.
%      A struct with a function handle, anonymous function or inline function in the 'h' field and
%      optionally additional data in the 'data' field. The function should accept three inputs: a
%      char array (either 'warning' or 'error'), a struct with the message, id, and stack, and the
%      optional additional data. The function(s) will be run before the error is actually thrown.
%      [default=[];]
%      If print_to_fid, print_to_obj, and print_to_fcn are all empty, this will have the effect of
%      suppressing every output except errors.
%      Array inputs are allowed.
%
% This function is aimed at providing a reliable method of reading a file. The backbone of this
% function is fread, supplemented by the fileread function. These work in slightly different ways
% and can be used under different circumstances. An attempt is made to detect the encoding (UTF-8
% or ANSI), apply the transcoding and returning the file as an n-by-1 cell array for files with
% n lines.
% You can redirect all outputs (errors only partially) to a file or a graphics object, or run a
% function based on the errors/warnings so you can more easily use this function in a GUI or allow
% it to write to a log file.
%
% Some input parameters can be used to mimic the readlines function, which was introduced in R2020a
% and returns a string vector instead of a cell array of character vectors.
%
% The test for being UTF-8 can fail. For files with chars in the 128:255 range, the test will often
% determine the encoding correctly, but it might fail, especially for files with encoding errors.
% Online files are much more limited than offline files. To avoid this the files are downloaded to
% tempdir() and deleted after reading. To avoid this the files are downloaded to tempdir() and
% deleted after reading. An additional fallback reads online files with webread/urlread, although
% this will often result in an incorrect output. This should only be relevant if there is no write
% access to the tempdir().
%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 4.0.0                                                         |%
%|  Date:    2021-05-19                                                    |%
%|  Author:  H.J. Wisselink                                                |%
%|  Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 ) |%
%|  Email = 'h_j_wisselink*alumnus_utwente_nl';                            |%
%|  Real_email = regexprep(Email,{'*','_'},{'@','.'})                      |%
%|                                                                         |%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%
% Compatibility considerations:
% Tested on several versions of Matlab (ML 6.5 and onward) and Octave (4.4.1 and onward), and on
% multiple operating systems (Windows/Ubuntu/MacOS). For the full test matrix, see the HTML doc.
% - The size of the char arrays may be different between Octave and Matlab. This is because Matlab
%   encodes characters internally with UTF-16, which means all 'normal' characters only take up a
%   single 16 bit value. A doc page seems to suggest Matlab uses UTF-8 to encode chars, but appears
%   to only be true for file interactions. If you want to include higher Unicode code points (e.g.
%   most emoji), some characters will require 2 elements in a char array. Octave use UTF-8 to
%   encode chars, but chars with values 128-255 are supported 'by accident'. This might change at
%   some point, but switching Octave to UTF-16 would require a lot of work, with the only
%   fundamental benefit being that size functions will return the same results between Matlab and
%   Octave. Judging by a discussion in the Octave bug tracker, I doubt this change will ever
%   happen.
% - It is therefore important to remember that a scalar char is not guaranteed to be a single
%   Unicode character, and that a single Unicode character is not guaranteed to be a single glyph. 
% - The readlines function was introduced in R2020b. It doesn't read to a cell of chars, but to a
%   string vector. The documentation implies that these two functions are functionally equivalent
%   (apart from that difference), but it seems to fail for characters beyond the BMP (Basic
%   Multilingual Plane). That means most emoji will fail. A future version of readlines might
%   correct this. When this bug is corrected
%   isequal(cellstr(readlines(filename)),readfile(filename)) should return true for all files.
%   Since R2021a readlines also supports reading online files.
% - Incorrect reading of files should only occur if the download to a temporary location fails.
%   (NB: this should be a rare occurence) Modern releases of Matlab (>=R2015a) are expected to read
%   every file correctly, except for ANSI files containing special characters. GNU Octave has
%   trouble with many ANSI files. Older releases of Matlab have the same results as Octave for ANSI
%   files, but also have issues with some UTF-8 files. Interestingly, R13 (v6.5) performs better on
%   ANSI files, but worse on UTF-8.

% Tested with 4 files with the following chars:
% list_of_chars_file1=[...
%     0032:0035 0037 0039:0042 0044:0059 0061 0063 0065:0091 0093 0096:0122 0160 0171 0173 0183 ...
%     0187:0189 0191:0193 0196 0200:0203 0205 0207 0209 0211 0212 0218 0224:0226 0228 0230:0235 ...
%     0237:0239 0241:0244 0246 0249:0253 8211 8212 8216:8218 8220:8222 8226 8230];
% list_of_chars_file2=[32:126 160:255 32 32 32];
% list_of_chars_file3=[...
%     0032:0126 0161:0163 0165 0167:0172 0174:0187 0191:0214 0216:0275 0278:0289 0292 0293 0295 ...
%     0298 0299 0304 0305 0308 0309 0313 0314 0317 0318 0321:0324 0327 0328 0336:0341 0344:0357 ...
%     0362:0369 0376:0382 0913:0929 0931:0974 0977 0984:0989 0991:0993 8211 8212 8216:8222 8224 ...
%     8225 8226 8230 8240 8249 8250 8260 8353 8356 8358 8361 8363 8364 8370 8482];
% list_of_chars4=[...
%    008986,009785,010084,128025,128512,128512,128513,128522,128550,128551,128552,128553,128555,...
%    128561,128578,128583,129343];
if nargin<1
    error('HJW:readfile:nargin','Incorrect number of input argument.')
end
if ~(nargout==0 || nargout==1) %might trigger 'MATLAB:TooManyOutputs' instead
    error('HJW:readfile:nargout','Incorrect number of output argument.')
end
[success,opts,ME]=readfile_parse_inputs(filename,varargin{:});
if ~success
    rethrow(ME)
else
    [filename,print_2,legacy,UseURLread,err_on_ANSI,EmptyLineRule,Whitespace,LineEnding,...
        FailMultiword_UTF16,WhitespaceRule]=...
        deal(opts.filename,opts.print_2__options,opts.legacy,opts.UseURLread,opts.err_on_ANSI,...
        opts.EmptyLineRule,opts.Whitespace,opts.LineEnding,opts.FailMultiword_UTF16,...
        opts.WhitespaceRule);
end

if opts.OfflineFile
    data=readfile_from_file(filename,LineEnding,print_2,err_on_ANSI);
else
    if ~legacy.allows_https && strcmpi(filename(1:8),'https://')
        warning_(print_2,'HJW:readfile:httpsNotSupported',...
            ['This implementation of urlread probably doesn''t allow https requests.',char(10),...
            'The next lines of code will probably result in an error.']) %#ok<CHARTEN>
    end
    str=readfile_from_URL(filename,UseURLread,print_2,LineEnding,err_on_ANSI);
    if isa(str,'cell') %file was read from temporary downloaded version
        data=str;
    else
        %This means the download failed. Some files will not work.
        invert=true;
        str=convert_from_codepage(str,invert);
        try ME=[]; %#ok<NASGU>
            [ii,isUTF8,converted]=UTF8_to_unicode(str); %#ok<ASGLU>
        catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
            if strcmp(ME.identifier,'HJW:UTF8_to_unicode:notUTF8')
                isUTF8=false;
            else
                error_(print_2,ME)
            end
        end
        if isUTF8
            str=unicode_to_char(converted);
        end
        if isa(LineEnding,'double') && isempty(LineEnding)
            data=char2cellstr(str);
        else
            data=char2cellstr(str,LineEnding);
        end
    end
end

%Determine the location of whitespace, but only if relevant.
if ~strcmp(EmptyLineRule,'read') || ~strcmp(WhitespaceRule,'preserve')
    L=cellfun('isempty',data);
    for n=find(~L).'
        % The cellfun call will only find completely empty lines, while readlines implicitly
        % considers lines with only whitespace empty.
        tmp=ismember(data{n},Whitespace);
        L(n)=all(tmp);
        if ~strcmp(WhitespaceRule,'preserve')
            % If there is only whitespace, take a shortcut by wiping the line now.
            if L(n),data{n}='';continue,end
            % If the first and last chars are whitespace, triming will have no effect.
            if ~tmp(1) && ~tmp(end),continue,end
            %Find the indices of non-whitespace.
            switch WhitespaceRule
                case 'trim'
                    inds=find(~tmp);inds=inds([1 end]);
                case 'trimleading'
                    %Use findND to extend the syntax for old Matlab releases.
                    inds=[findND(~tmp,1) numel(tmp)];
                case 'trimtrailing'
                    %Use findND to extend the syntax for old Matlab releases.
                    inds=[1 findND(~tmp,1,'last')];
            end
            data{n}=data{n}(inds(1):inds(2));
        end
    end
end
if ~strcmp(EmptyLineRule,'read')
    switch EmptyLineRule
        %To allow the expanded syntax for find(), the findND() function is used instead, as that
        %extends the syntax for find() on old releases of Matlab.
        case 'skip'
            data(L)=[];
        case 'error'
            if any(L)
                error_(print_2,'HJW:readfile:EmptyLinesRuleError',...
                    'Unexpected empty line detected on row %d',findND(L,1))
            end
        case 'skipleading'
            if L(1)
                ind=1:(findND(~L,1,'first')-1);
                data(ind)=[];
            end
        case 'skiptrailing'
            if L(end)
                ind=(1+findND(~L,1,'last')):numel(L);
                data(ind)=[];
            end
    end
end
persistent isOctave,if isempty(isOctave),isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;end
if FailMultiword_UTF16
    %R2020b and R2021a fail for multiword UTF16 characters, rendering them as char(26). To keep
    %complete equivalence, that behavior is replicated here.
    %The bit-pattern is 110110xx_xxxxxxxx 110111xx_xxxxxxxx, so we can simply detect any value
    %between 55296 and 56319. For Octave we can check if there are 4-byte characters.
    for n=1:numel(data)
        if isOctave
            if any(data{n}>=240)
                %Now we need to properly convert to UTF-32, replace by 26 and convert back.
                data{n}=replace_multiword_UTF16_by_26(data{n});
            end
        else
            if any(data{n}>=55296 & data{n}<=56319)
                %Now we need to properly convert to UTF-32, replace by 26 and convert back.
                data{n}=replace_multiword_UTF16_by_26(data{n});
            end
        end
    end
end
end
function out=replace_multiword_UTF16_by_26(in)
%Replace all multiword UTF-16 (i.e. U+10000 to U+10FFFF) with char(26)
persistent isOctave,if isempty(isOctave),isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;end
if isOctave
    unicode=UTF8_to_unicode(in);
else
    unicode=UTF16_to_unicode(in);
end
%Perform replacement.
unicode(unicode>=65536)=26;
%Convert back to proper encoding
if isOctave
    out=char(unicode_to_UTF8(unicode));
else
    out=char(unicode_to_UTF16(unicode));
end
end
function str=readfile_from_file(filename,LineEnding,print_2,err_on_ANSI)
persistent isOctave,if isempty(isOctave),isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;end
persistent ME_file_access_FormatSpec
if isempty(ME_file_access_FormatSpec)
    if isOctave,runtime='Octave';else,runtime='Matlab';end
    ME_file_access_FormatSpec=sprintf(['%s could not read the file %%s.\n',...
        'The file doesn''t exist or is not readable.\n',...
        '(Note that for online files, only http and https is supported.)'],runtime);
end
ME_file_access=struct('identifier','HJW:readfile:ReadFail','message',...
        sprintf(ME_file_access_FormatSpec,filename));

fid=fopen(filename,'rb');
if fid<0,error_(print_2,ME_file_access),end
data=fread(fid,'uint8');
fclose(fid);
data=data.';
try ME=[]; %#ok<NASGU>
    isUTF8=true;
    converted=UTF8_to_unicode(data);
catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
    if strcmp(ME.identifier,'HJW:UTF8_to_unicode:notUTF8')
        isUTF8=false;
        if err_on_ANSI
            error_(print_2,'HJW:readfile:notUTF8',...
                'The provided file "%s" is not a correctly encoded UTF-8 file.',filename)
        end
    else
        error_(print_2,ME)
    end
end

if isOctave
    if isUTF8
        str=converted;
    else
        try str=fileread(filename);catch,error_(print_2,ME_file_access),end
        str=convert_from_codepage(str);
    end
else
    if ispc
        if isUTF8
            str=converted;
        else
            if ifversion('<',7)
                try str=fileread(filename);catch,error_(print_2,ME_file_access),end
                str=convert_from_codepage(str);
            else
                try str=fileread(filename);catch,error_(print_2,ME_file_access),end
            end
        end
    else %This assumes Mac will work the same as Ubuntu.
        if isUTF8
            str=converted;
        else
            str=convert_from_codepage(data);
        end
    end
end

%Remove UTF BOM (U+FEFF) from text.
if numel(str)>=1 && double(str(1))==65279,str(1)=[];end
%Convert back to a char and split to a cellstr.
str=unicode_to_char(str);
if isa(LineEnding,'double') && isempty(LineEnding)
    str=char2cellstr(str);
else
    str=char2cellstr(str,LineEnding);
end
end
function str=readfile_from_URL(url,UseURLread,print_to,LineEnding,err_on_ANSI)
%Read the contents of a file to a char array.
%
%Attempt to download to the temp folder, read the file, then delete it.
%If that fails, read to a char array with urlread/webread.
try
    RevertToUrlread=false;%in case the saving+reading fails
    %Generate a random file name in the temp folder.
    fn=tmpname('readfile_from_URL_tmp_','.txt');
    try
        %Try to download (with 'raw' as ContentType to prevent parsing of XML/JSON/etc by websave).
        if UseURLread,fn=urlwrite(url,fn); %#ok<URLWR>
        else,         fn= websave(fn,url,weboptions('ContentType','raw'));end
        
        %Try to read
        str=readfile_from_file(fn,LineEnding,print_to,err_on_ANSI);
    catch
        RevertToUrlread=true;
    end
    
    %Delete the temp file
    try if exist(fn,'file'),delete(fn);end,catch,end
    
    if RevertToUrlread,error('revert to urlread'),end
catch
    %Read to a char array and let these functions throw an error in case of HTML errors and/or
    %missing connectivity.
    try ME=[]; %#ok<NASGU>
        %Use 'raw' as ContentType to prevent parsing of XML/JSON/etc by webread.
        if UseURLread,str=urlread(url);else%#ok<URLRD>
            str=webread(url,weboptions('ContentType','raw'));end 
    catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
        error_(print_to,ME)
    end
end
end
function [success,options,ME]=readfile_parse_inputs(filename,varargin)
%Parse the inputs of the readfile function
% It returns a success flag, the parsed options, and an ME struct.
% As input, the options should either be entered as a struct or as Name,Value pairs. Missing fields
% are filled from the default.

%Pre-assign outputs.
success=false;
options=struct;
ME=struct('identifier','','message','');

%Test the required input.
[valid,filename]=filename_is_valid(filename);%This will covert string to char.
try
    OnlineFile= ...
        (  strcmpi(filename(1:7),'http://') ...
        || strcmpi(filename(1:8),'https://'));
    if OnlineFile
        %Test if it is long enough to be a proper URL.
        if numel(filename)<10,error('trigger'),end
    else
        %Offline files must adhere to the standards of this check.
        if ~valid,error('trigger'),end
    end
catch
    ME.identifier='HJW:readfile:IncorrectInput';
    ME.message='The file name must be a non-empty char or a scalar string.';
    return
end

persistent default print_2__default_options
if isempty(default)
    %The regexp split option was introduced in R2007b.
    legacy.split = ifversion('<','R2007b','Octave','<',4);
    %Change this line when Octave does support https.
    legacy.allows_https=ifversion('>',0,'Octave','<',0);
    default.legacy=legacy;
    %Test if either webread(), websave(), or weboptions() are missing.
    try no_webread=isempty(which(func2str(@webread   )));catch,no_webread=true;end
    try no_websave=isempty(which(func2str(@websave   )));catch,no_websave=true;end
    try no_webopts=isempty(which(func2str(@weboptions)));catch,no_webopts=true;end
    default.UseURLread= no_webread || no_websave || no_webopts;
    
    print_2__default_options=struct;
    print_2__default_options.print_to_con=true;
    print_2__default_options.print_to_fid=[];
    print_2__default_options.print_to_obj=[];
    print_2__default_options.print_to_fcn=[];
    default.print_2__options=validate_print_to__options(print_2__default_options,ME);
    
    default.err_on_ANSI=false;
    %R2020b has a bug where readlines fails for chars outside the BMP (e.g. most emoji).
    default.FailMultiword_UTF16=false;
    default.EmptyLineRule='read';
    %The Whitespace parameter contains most characters reported by isspace, plus delete characters
    %and no break spaces. This is different from the default for readlines.
    %To make sure all these code points are encoded in char correctly, we need to use
    %unicode_to_char. The reason for this is that Octave uses UTF-8.
    default.Whitespace=unicode_to_char([8 9 28:32 160 5760 8192:8202 8239 8287 12288]);
    default.DefaultLineEnding=true;
    default.LineEnding=[];%(equivalent to {'\r\n','\n','\r'}, the order matters for char2cellstr)
    default.WhitespaceRule='preserve';
    
    default.UseReadlinesDefaults=false;
    default.ReadlinesDefaults.FailMultiword_UTF16=true;
    default.ReadlinesDefaults.Whitespace=sprintf(' \b\t');
    %Fill any missing fields.
    fn=fieldnames(default);
    for k=1:numel(fn)
        if ~isfield(options,fn(k)),default.ReadlinesDefaults.(fn{k})=default.(fn{k});end
    end
end
%The required inputs are checked, so now we need to return the default options if there are no
%further inputs.
if nargin==1
    options=default;
    success=true;
    
    options.OfflineFile=~OnlineFile;
    options.filename=filename;
    return
end

%Test the optional inputs.
struct_input=   nargin==2        && isa(varargin{1},'struct');
NameValue_input=mod(nargin,2)==1 && all(...
    cellfun('isclass',varargin(1:2:end),'char'  ) | ...
    cellfun('isclass',varargin(1:2:end),'string')   );
if ~( struct_input || NameValue_input )
    ME.message=['The second input (options) is expected to be either a struct, ',char(10),...
        'or consist of Name,Value pairs.']; %#ok<CHARTEN>
    ME.identifier='HJW:readfile:incorrect_input_options';
    return
end
if NameValue_input
    %Convert the Name,Value to a struct.
    for n=1:2:numel(varargin)
        try
            options.(varargin{n})=varargin{n+1};
        catch
            ME.message='Parsing of Name,Value pairs failed.';
            ME.identifier='HJW:readfile:incorrect_input_NameValue';
            return
        end
    end
else
    options=varargin{1};
end
%Do some special parsing for the print_2__options
print_2__replace_default_options=false;print_2__options=print_2__default_options;
%Create a copy of the defaults so they can be overwritten by the ReadlinesDefaults.
default__copy=default;
fn=fieldnames(options);
for k=1:numel(fn)
    curr_option=fn{k};
    item=options.(curr_option);
    ME.identifier=['HJW:readfile:incorrect_input_opt_' lower(curr_option)];
    switch curr_option
        case 'UseURLread'
            [passed,item]=test_if_scalar_logical(item);
            if ~passed
                ME.message='UseURLread should be either true or false';
                return
            end
            %Force the use of urlread/urlwrite if websave/webread are not available.
            options.UseURLread=item || default.UseURLread;
        case 'err_on_ANSI'
            [passed,item]=test_if_scalar_logical(item);
            if ~passed
                ME.message='err_on_ANSI should be either true or false';
                return
            end
            options.err_on_ANSI=item;
        case 'print_to_fid'
            print_2__replace_default_options=true;
            print_2__options.print_to_fid=item;
        case 'print_to_obj'
            print_2__replace_default_options=true;
            print_2__options.print_to_obj=item;
        case 'print_to_fcn'
            print_2__replace_default_options=true;
            print_2__options.print_to_fcn=item;
        case 'print_to_con'
            print_2__replace_default_options=true;
            print_2__options.print_to_con=item;
        case 'EmptyLineRule'
            if isa(item,'string')
                if numel(item)~=1,item=[]; % This will trigger an error.
                else,item=char(item);end   % Convert a scalar string to char.
            end
            if isa(item,'char'),item=lower(item);end
            if ~isa(item,'char') || ...
                    ~ismember(item,{'read','skip','error','skipleading','skiptrailing'})
                ME.message='EmptyLineRule must be a char or string with a specific value.';
                return
            end
            options.EmptyLineRule=item;
        case 'Whitespace'
            %Cellstr input is converted to a char array with sprintf.
            try
                switch class(item)
                    case 'string'
                        if numel(item)~=1,error('trigger error'),end
                        item=char(item);
                    case 'cell'
                        for n=1:numel(item),item{n}=sprintf(item{n});end
                        item=horzcat(item{:});
                    case 'char'
                        %Nothing to check or do here.
                    otherwise
                        error('trigger error')
                end
                options.Whitespace=item;
            catch
                ME.message=['The Whitespace parameter must be a char vector, string scalar or ',...
                    'cellstr.\nA cellstr input must be parsable by sprintf.'];
                return
            end
        case 'WhitespaceRule'
            if isa(item,'string')
                if numel(item)~=1,item=[]; % This will trigger an error.
                else,item=char(item);end   % Convert a scalar string to char.
            end
            if isa(item,'char'),item=lower(item);end
            if ~isa(item,'char') || ...
                    ~ismember(item,{'preserve','trim','trimleading','trimtrailing'})
                ME.message='WhitespaceRule must be a char or string with a specific value.';
                return
            end
            options.WhitespaceRule=item;
        case 'LineEnding'
            %character vector  - literal
            %string scalar  - literal
            %cell array of character vectors  - parse by sprintf
            %string array  - parse by sprintf
            err=false;
            if isa(item,'string')
                item=cellstr(item);
                if numel(item)==1,item=item{1};end % Convert scalar string to char.
            end
            if isa(item,'cell')
                try for n=1:numel(item),item{n}=sprintf(item{n});end,catch,err=true;end
            elseif isa(item,'char')
                item={item};% Wrap char in a cell.
            else
                err=true; % This catches [] as well, while iscellstr wouldn't.
            end
            if err || ~iscellstr(item)
                ME.message=['The LineEnding parameter must be a char vector, a string or a ',...
                    'cellstr.\nA cellstr or string vector input must be parsable by sprintf.'];
                return
            end
            if isequal(item,{char(10) char(13) char([13 10])}) %#ok<CHARTEN>
                options.LineEnding=[];
            else
                options.LineEnding=item;
            end
        case 'UseReadlinesDefaults'
            [passed,item]=test_if_scalar_logical(item);
            if ~passed
                ME.message='UseReadlinesDefaults should be either true or false';
                return
            end
            options.UseReadlinesDefaults=item;
            if options.UseReadlinesDefaults
                default__copy=default.ReadlinesDefaults;
            end
    end
end

options.OfflineFile=~OnlineFile;
options.filename=filename;

if print_2__replace_default_options
    ME.identifier='HJW:readfile:incorrect_input_opt_print_to_';
    [item,ME]=validate_print_to__options(print_2__options,ME);
    if isempty(item),return,end
    options.print_2__options=item;
else
    %If none of them were set, we can use the defaults.
    options.print_2__options=default__copy.print_2__options;
end

%Fill any missing fields.
fn=fieldnames(default__copy);
for k=1:numel(fn)
    if ~isfield(options,fn(k))
        options.(fn{k})=default__copy.(fn{k});
    end
end
success=true;ME=[];
end
function varargout=regexp_outkeys(str,expression,varargin)
%Regexp with outkeys in old releases
%
% On older versions of Matlab the regexp function did not allow you to specify the output keys.
% This function has an implementation of the 'split' and 'match' output keys, so they can be used
% on any version of Matlab or GNU Octave.
% On releases where these outkeys are available, the builtin is called.
%
% Syntax:
%   out=regexp_outkeys(str,expression,outkey);
%   [out1,...,outN]=regexp_outkeys(str,expression,outkey1,...,outkeyN);
%   [___,startIndex,endIndex]=regexp_outkeys(___);
%
%Example:
% str='lorem1 ipsum1.2 dolor3 sit amet 99 ';
% words=regexp_outkeys(str,' ','split')
% numbers=regexp_outkeys(str,'[0-9.]*','match')
% [white,end1,start,end2]=regexp_outkeys(str,' ','match','end')
%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 1.1.0                                                         |%
%|  Date:    2021-10-20                                                    |%
%|  Author:  H.J. Wisselink                                                |%
%|  Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 ) |%
%|  Email = 'h_j_wisselink*alumnus_utwente_nl';                            |%
%|  Real_email = regexprep(Email,{'*','_'},{'@','.'})                      |%
%|                                                                         |%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%
% Tested on several versions of Matlab (ML 6.5 and onward) and Octave (4.4.1 and onward), and on
% multiple operating systems (Windows/Ubuntu/MacOS). For the full test matrix, see the HTML doc.
% Compatibility considerations:
% - Only the 'match', 'split', 'start', and 'end' options are supported. The additional options
%   provided by regexp are not implemented.
% - Cell array input is not supported.

if nargin<2
    error('HJW:regexp_outkeys:SyntaxError',...
        'No supported syntax used: at least 3 inputs expected.')
    %As an undocumented feature this can also return s1,s2 without any outkeys specified.
end
if ~(ischar(str) && ischar(expression))
    %The extra parameters in varargin are checked inside the loop.
    error('HJW:regexp_outkeys:InputError',...
        'All inputs must be char vectors.')
end
if nargout>nargin
    error('HJW:regexp_outkeys:ArgCount',...
        'Incorrect number of output arguments. Check syntax.')
end

persistent legacy errorstr
if isempty(legacy)
    %The legacy struct contains the implemented options as field names. It is used in the error
    %message.
    %It is assumed that all Octave versions later than 4.0 support the expanded output, and all
    %earlier versions do not, even if it is very likely most versions will support it.

    %The switch to find matches was introduced in R14 (along with the 'tokenExtents', 'tokens' and
    %'names' output switches).
    legacy.match = ifversion('<','R14','Octave','<',4);

    %The split option was introduced in R2007b.
    legacy.split = ifversion('<','R2007b','Octave','<',4);

    fn=fieldnames(legacy);
    errorstr=['Extra regexp output type not implemented,',char(10),'only the following',...
        ' types are implemented:',char(10),sprintf('%s, ',fn{:})]; %#ok<CHARTEN>
    errorstr((end-1):end)='';%Remove trailing ', '

    legacy.any = legacy.match || legacy.split;
end

if legacy.any || nargin==2
    %Determine s1 and s2 only once for the legacy implementations.
    [s1,s2]=regexp(str,expression);
end

if nargin==2
    varargout={s1,s2};return
end

%Pre-allocate output.
varargout=cell(size(varargin));
for param=1:(nargin-2)
    if ~ischar(varargin{param})
        error('HJW:regexp_outkeys:InputError',...
            'All inputs must be char vectors.')
    end
    switch lower(varargin{param})
        case 'match'
            if legacy.match
                %Legacy implementation.
                match=cell(1,numel(s1));
                for n=1:numel(s1)
                    match{n}=str(s1(n):s2(n));
                end
            else
                [match,s1,s2]=regexp(str,expression,'match');
            end
            varargout{param}=match;
        case 'split'
            if legacy.split
                %Legacy implementation.
                split=cell(1,numel(s1)+1);
                start_index=[s1 numel(str)+1];
                stop_index=[0 s2];
                for n=1:numel(start_index)
                    split{n}=str((stop_index(n)+1):(start_index(n)-1));
                end
            else
                [split,s1,s2]= regexp(str,expression,'split');
            end
            varargout{param}=split;
        case 'start'
            varargout{param}=s1;
        case 'end'
            varargout{param}=s2;
        otherwise
            error('HJW:regexp_outkeys:NotImplemented',errorstr)
    end
end
if nargout>param
    varargout(param+[1 2])={s1,s2};
end
end
function [stmt,strct,failed,ME]=sqlite3_param_binding_parse(stmt,strct)
%Parse the SQL statement, assuming an insert statement (or equivalent).
%
%The end result is that the statement will be edited, replacing named/numbered parameters with ?
%and adapting the input struct accordingly.
%
%This function also allows duplicative references:
%  stmt ='INSERT INTO test (some_real, some_int) VALUES (?1, ?1);';
%  strct=struct('number',3);
failed=false;ME=struct;
%NB: numbered parameters are 1-indexed!

%Build the regular expression that matches the parameters.
%      ? or ?NNN        :VVV or @VVV or $VVV
p=['('  '\?\d*'  '|'  '[:@$][A-Z_][A-Z_0-9]*'  ')'];
%                          first param    subsequent params
expr=['\)\s*VALUES\s*\('  '\s*' p '\s*'  '(\s*,\s*' p '\s*)*'  '\)'];


[s1,s2]=regexp(upper(stmt),expr);
if numel(s1)~=1
    failed=true;
    ME.identifier='HJW:sqlite3:ParamaterDetectionFailed';
    ME.message='Detection of parameters for binding failed.';
    return
end
s1=s1-1+strfind(stmt(s1:s2),'(');%Account for ')VALUE('.

try
    %Match named and numbered parameters.
    ind=bind_params(stmt(s1:s2),strct);
catch
    failed=true;
    ME.identifier='HJW:sqlite3:ParamaterBindingMismatch';
    ME.message=[...
        'The number of unnamed parameters (''?'') does not match the number of ',char(10),...
        'unmatched fields (i.e. fields not referenced by '':NAME'' or ''?index'').']; %#ok<CHARTEN>
    return
end

try
    strct=reorder_fields(strct,ind);
catch
    failed=true;
    ME.identifier='HJW:sqlite3:ParameterReorderingFailed';
    ME.message='Failed to reorder struct to match parameter binding.';
    return
end

%Replace all parameters with '?'.
params=['(?' repmat(',?',1,sum(stmt(s1:s2)==',')) ')'];
stmt(s1:(s1-1+numel(params))   )=params;
stmt(   (s1  +numel(params)):s2)='';
end
function ind=bind_params(params,strct)
%Match named and numbered parameters.
params=regexp_outkeys(...
    strrep(params(2:(end-1)),' ',''),...
    ',','split').';
ind=NaN*zeros(size(params));

%Parse numbered parameters, marking '?' as 0.
numbered_params=str2double(strrep(params,'?','0'));
ind(~isnan(numbered_params))=numbered_params(~isnan(numbered_params));

%Match named parameters.
fields=fieldnames(strct);
[L,ii]=ismember(regexprep(params,'[:@$]',''),fields);
ind(L)=ii(L);
L=ind==0;unmatched_fields=setdiff(1:numel(fields),ind);
if any(L) && numel(unmatched_fields)~=sum(L)
    %This would lead to an ambiguous match.
    error('trigger error')
end
ind(L)=unmatched_fields;
end
function strct=reorder_fields(strct,ind)
%Reorder the fields of the struct to match the indices (duplicating/removing fields as needed).

fields=fieldnames(strct);
persistent legacy
if isempty(legacy)
    try   legacy.histc=isempty(which(func2str(@histcounts)));
    catch,legacy.histc=true;end
end

%If all fields are used once we can use orderfields to return the original struct in the correct
%order. Otherwise we might need to create dummy fields.
if isequal(1:numel(fields),sort(ind))
    strct=orderfields(strct,ind);
else
    try
        %Reorder the data.
        c=struct2cell(strct);
        c_ind={ind};c_ind(2:ndims(c))={':'};
        c=c(c_ind{:});
        
        %Reorder the fields (creating dummy names for duplicates).
        fn=fields(ind);
        
        %Replace later occurrences of duplicated fields with dummy variable names.
        k=1;%Enter loop.
        while ~isempty(k)
            [ignore,ind_]=ismember(fn,unique(fn)); %#ok<ASGLU>
            if legacy.histc,k=find(histc(     ind_,0.5:numel(fn))>1); %#ok<HISTC>
            else           ,k=find(histcounts(ind_,0.5:numel(fn))>1);
            end
            for k=k
                k2=find(ind_==k);
                fn{k2(2)}=dummy_field_name(fn{k2(2)});
            end
        end
        strct=cell2struct(c,fn);
    catch
        %As a fallback, we can copy all the fields to a new struct in the correct order. The other
        %method probably has a much better performance, but with this it is easy to see how it
        %works. This method should have no edge cases.
        
        %Initialize the output struct to be the same size as the input struct.
        fn=dummy_field_name;
        out=struct(fn,cell(size(strct)));
        out=rmfield(out,fn);
        
        %Loop through the fields to create the output struct.
        for k=ind(:).'
            fn=dummy_field_name(fields{k});
            [out.(fn)]=deal(strct(:).(fields{k}));
        end
        strct=out;
    end
end
end
function fn=dummy_field_name(original)
if nargin==0,original='';end
x='ABCDEFGHIJKLMNPOQRSTUVWXYZ0123456789';
tmp=x(ceil(36*rand(1,34)));
%Maximum name length is 63, 63-(7+20+2)
fn=['dummy__' original(1:min(end,20)) '__' tmp];
end
function [success,opts,ME]=sqlite3_parse_inputs(db_or_opts,SQL_stmnt,varargin)
success=false;ME=[];
if isa(db_or_opts,'struct')
    if numel(db_or_opts)~=1
        opts=struct;ME.message='The option struct should be scalar';
        ME.identifier='HJW:sqlite3:NonScalarOptionstruct';return
    end
    opts=db_or_opts;
    opts.SQL_statement=SQL_stmnt;
else
    opts=struct('filename',db_or_opts,'SQL_statement',SQL_stmnt);
end

%Attempt to match the struct to the available options. This will return a struct with the same
%fields as the default option struct.
try [opts,replaced]=parse_NameValue(sqlite3_parse_inputs_defaults,opts);
catch ME;if isempty(ME),ME=lasterror;end,return,end %#ok<LERR>

if any(ismember({'isOctave','error_on_CLI'},replaced))
    % Revert these replacements. These are meant to be const.
    defaults=sqlite3_parse_inputs_defaults;
    opts.isOctave=defaults.isOctave;
    opts.error_on_CLI=defaults.error_on_CLI;
end

%Attempt to parse the error redirection options (this generates ME.message on fail).
[opts.print_to,ME]=validate_print_to__options(opts,ME);
if isempty(opts.print_to)
    ME.identifier='HJW:sqlite3:Incorrect_print_to_option';return
end

%Check the filename.
[valid,opts.filename]=filename_is_valid(opts.filename);
if ~valid
    ME.identifier='HJW:sqlite3:FilenameInvalid';
    ME.message='Invalid file name.';
    return
end

%Check the SQL statement.
if isa(opts.SQL_statement,'string') && numel(opts.SQL_statement)==1
    %Convert a scalar string to a char array.
    opts.SQL_statement=char(opts.SQL_statement);
end
var=opts.SQL_statement;
if ~isa(var,'char') || ~(size(var,1)==1 && size(var,2)>0)
    ME.identifier='HJW:sqlite3:IncorrectDataType';
    ME.message='The SQL statement must be a char vector or a scalar string.';
    return
end

%Check if the extra input matches the placeholders.
if nargin==3
    [stmt,strct,failed,ME]=sqlite3_param_binding_parse(opts.SQL_statement,varargin{1});
    if failed,return,end
    [opts.SQL_statement,opts.input_struct]=deal(stmt,strct);
end

% Check if the CLI is available and if it is selected.
if ismember({'use_CLI'},replaced)
    [passed,opts.use_CLI]=test_if_scalar_logical(opts.use_CLI);
    if ~passed
        ME.identifier='HJW:sqlite3:InvalidParameter__use_CLI';
        ME.message='The use_CLI parameter must be a scalar logical.';
        return
    elseif opts.use_CLI && opts.error_on_CLI
        ME.identifier='HJW:sqlite3:NoCLI';
        ME.message='There is no CLI available.';
        return
    end
end

success=true;
end
function opts=sqlite3_parse_inputs_defaults
%Create a struct with default values.
persistent opts_
if isempty(opts_)
    isOctave=ifversion('<',0,'Octave','>',0);
    opts_=struct(...
        'filename','',...
        'SQL_statement','',...
        'input_struct','',...
        'isOctave',isOctave,...
        'print_to_con',[],...
        'print_to_fid',[],...
        'print_to_obj',[],...
        'print_to_fcn',[],...
        'use_CLI',false,...
        'error_on_CLI',true);
    % If and when the CLI is re-implemented, the error_on_CLI parameter will most likely be
    % replaced by a test that confirms the CLI is working properly (the result of which will be
    % stored in a persistent).
end
opts=opts_;
end
function str=stringtrim(str)
%Extend strtrim to remove double spaces as well.
if exist('strtrim','builtin')
    str=strtrim(str);
else
    %ML6.5
    if numel(str)==0,return,end
    L=isspace(str);
    if L(end)
        %The last character is whitespace, so trim the end.
        idx=find(~L);
        if isempty(idx)
            %The char only contains whitespace.
            str='';return
        end
        str((idx(end)+1):end)='';
    end
    if isempty(str),return,end
    if L(1)
        %The first character is whitespace, so trim the start.
        idx=find(~L);
        str(1:(idx(1)-1))='';
    end
end
removed_double_spaces=inf;
while removed_double_spaces~=0
    length1=length(str);
    str=strrep(str,'  ',' ');
    length2=length(str);
    removed_double_spaces=length1-length2;
end
end
function [isLogical,val]=test_if_scalar_logical(val)
%Test if the input is a scalar logical or convertible to it.
%The char and string test are not case sensitive.
%(use the first output to trigger an input error, use the second as the parsed input)
%
% Allowed values:
%- true or false
%- 1 or 0
%- 'on' or 'off'
%- matlab.lang.OnOffSwitchState.on or matlab.lang.OnOffSwitchState.off
%- 'enable' or 'disable'
%- 'enabled' or 'disabled'
persistent states
if isempty(states)
    states={true,false;...
        1,0;...
        'on','off';...
        'enable','disable';...
        'enabled','disabled'};
    try
        states(end+1,:)=eval('{"on","off"}'); %#ok<EVLCS> 
    catch
    end
end
isLogical=true;
try
    if isa(val,'char') || isa(val,'string')
        try val=lower(val);catch,end
    end
    for n=1:size(states,1)
        for m=1:2
            if isequal(val,states{n,m})
                val=states{1,m};return
            end
        end
    end
    if isa(val,'matlab.lang.OnOffSwitchState')
        val=logical(val);return
    end
catch
end
isLogical=false;
end
function tf=TestFolderWritePermission(f)
%Returns true if the folder exists and allows Matlab to write files.
%An empty input will generally test the pwd.
%
%examples:
%  fn='foo.txt';if ~TestFolderWritePermission(fileparts(fn)),error('can''t write!'),end

if ~( isempty(f) || exist(f,'dir') )
    tf=false;return
end

fn='';
while isempty(fn) || exist(fn,'file')
    %Generate a random file name, making sure not to overwrite any existing file.
    %This will try to create a file without an extension.
    [ignore,fn]=fileparts(tmpname('write_permission_test_','.txt')); %#ok<ASGLU>
    fn=fullfile(f,fn);
end
try
    %Test write permission.
    fid=fopen(fn,'w');fprintf(fid,'test');fclose(fid);
    delete(fn);
    tf=true;
catch
    %Attempt to clean up.
    if exist(fn,'file'),try delete(fn);catch,end,end
    tf=false;
end
end
function str=tmpname(StartFilenameWith,ext)
%Inject a string in the file name part returned by the tempname function.
if nargin<1,StartFilenameWith='';end
if ~isempty(StartFilenameWith),StartFilenameWith=[StartFilenameWith '_'];end
if nargin<2,ext='';else,if ~strcmp(ext(1),'.'),ext=['.' ext];end,end
str=tempname;
[p,f]=fileparts(str);
str=fullfile(p,[StartFilenameWith f ext]);
end
function str=unicode_to_char(unicode,encode_as_UTF16)
%Encode Unicode code points with UTF-16 on Matlab and UTF-8 on Octave.
%Input is either implicitly or explicitly converted to a row-vector.

persistent isOctave,if isempty(isOctave),isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;end
if nargin==1
    encode_as_UTF16=~CharIsUTF8;
end
if encode_as_UTF16
    if all(unicode<65536)
        str=uint16(unicode);
        str=reshape(str,1,numel(str));%Convert explicitly to a row-vector.
    else
        %Encode as UTF-16.
        [char_list,ignore,positions]=unique(unicode); %#ok<ASGLU>
        str=cell(1,numel(unicode));
        for n=1:numel(char_list)
            str_element=unicode_to_UTF16(char_list(n));
            str_element=uint16(str_element);
            str(positions==n)={str_element};
        end
        str=cell2mat(str);
    end
    if ~isOctave
        str=char(str);%Conversion to char could trigger a conversion range error in Octave.
    end
else
    if all(unicode<128)
        str=char(unicode);
        str=reshape(str,1,numel(str));%Convert explicitly to a row-vector.
    else
        %Encode as UTF-8
        [char_list,ignore,positions]=unique(unicode); %#ok<ASGLU>
        str=cell(1,numel(unicode));
        for n=1:numel(char_list)
            str_element=unicode_to_UTF8(char_list(n));
            str_element=uint8(str_element);
            str(positions==n)={str_element};
        end
        str=cell2mat(str);
        str=char(str);
    end
end
end
function str=unicode_to_UTF16(unicode)
%Convert a single character to UTF-16 bytes.
%
%The value of the input is converted to binary and padded with 0 bits at the front of the string to
%fill all 'x' positions in the scheme.
%See https://en.wikipedia.org/wiki/UTF-16
%
% 1 word (U+0000 to U+D7FF and U+E000 to U+FFFF):
%  xxxxxxxx_xxxxxxxx
% 2 words (U+10000 to U+10FFFF):
%  110110xx_xxxxxxxx 110111xx_xxxxxxxx
if unicode<65536
    str=unicode;return
end
U=double(unicode)-65536;%Convert to double for ML6.5.
U=dec2bin(U,20);
str=bin2dec(['110110' U(1:10);'110111' U(11:20)]).';
end
function str=unicode_to_UTF8(unicode)
%Convert a single character to UTF-8 bytes.
%
%The value of the input is converted to binary and padded with 0 bits at the front of the string to
%fill all 'x' positions in the scheme.
%See https://en.wikipedia.org/wiki/UTF-8
if numel(unicode)>1,error('this should only be used for single characters'),end
if unicode<128
    str=unicode;return
end
persistent pers
if isempty(pers)
    pers=struct;
    pers.limits.lower=hex2dec({'0000','0080','0800', '10000'});
    pers.limits.upper=hex2dec({'007F','07FF','FFFF','10FFFF'});
    pers.scheme{2}='110xxxxx10xxxxxx';
    pers.scheme{2}=reshape(pers.scheme{2}.',8,2);
    pers.scheme{3}='1110xxxx10xxxxxx10xxxxxx';
    pers.scheme{3}=reshape(pers.scheme{3}.',8,3);
    pers.scheme{4}='11110xxx10xxxxxx10xxxxxx10xxxxxx';
    pers.scheme{4}=reshape(pers.scheme{4}.',8,4);
    for b=2:4
        pers.scheme_pos{b}=find(pers.scheme{b}=='x');
        pers.bits(b)=numel(pers.scheme_pos{b});
    end
end
bytes=find(pers.limits.lower<=unicode & unicode<=pers.limits.upper);
str=pers.scheme{bytes};
scheme_pos=pers.scheme_pos{bytes};
b=dec2bin(unicode,pers.bits(bytes));
str(scheme_pos)=b;
str=bin2dec(str.').';
end
function unicode=UTF16_to_unicode(UTF16)
%Convert UTF-16 to the code points stored as uint32
%
%See https://en.wikipedia.org/wiki/UTF-16
%
% 1 word (U+0000 to U+D7FF and U+E000 to U+FFFF):
%  xxxxxxxx_xxxxxxxx
% 2 words (U+10000 to U+10FFFF):
%  110110xx_xxxxxxxx 110111xx_xxxxxxxx

persistent isOctave,if isempty(isOctave),isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;end
UTF16=uint32(UTF16);

multiword= UTF16>55295 & UTF16<57344; %0xD7FF and 0xE000
if ~any(multiword)
    unicode=UTF16;return
end

word1= find( UTF16>=55296 & UTF16<=56319 );
word2= find( UTF16>=56320 & UTF16<=57343 );
try
    d=word2-word1;
    if any(d~=1) || isempty(d)
        error('trigger error')
    end
catch
    error('input is not valid UTF-16 encoded')
end

%Binary header:
% 110110xx_xxxxxxxx 110111xx_xxxxxxxx
% 00000000 01111111 11122222 22222333
% 12345678 90123456 78901234 56789012
header_bits='110110110111';header_locs=[1:6 17:22];
multiword=UTF16([word1.' word2.']);
multiword=unique(multiword,'rows');
S2=mat2cell(multiword,ones(size(multiword,1),1),2);
unicode=UTF16;
for n=1:numel(S2)
    bin=dec2bin(double(S2{n}))';
    
    if ~strcmp(header_bits,bin(header_locs))
        error('input is not valid UTF-16 encoded')
    end
    bin(header_locs)='';
    if ~isOctave
        S3=uint32(bin2dec(bin  ));
    else
        S3=uint32(bin2dec(bin.'));%Octave needs an extra transpose.
    end
    S3=S3+65536;% 0x10000
    %Perform actual replacement.
    unicode=PatternReplace(unicode,S2{n},S3);
end
end
function [unicode,isUTF8,assumed_UTF8]=UTF8_to_unicode(UTF8,print_to)
%Convert UTF-8 to the code points stored as uint32
% Plane 16 goes up to 10FFFF, so anything larger than uint16 will be able to hold every code point.
%
% If there a second output argument, this function will not return an error if there are encoding
% error. The second output will contain the attempted conversion, while the first output will
% contain the original input converted to uint32.
%
% The second input can be used to also print the error to a GUI element or to a text file.
if nargin<2,print_to=[];end
return_on_error= nargout==1 ;

UTF8=uint32(reshape(UTF8,1,[]));% Force row vector.
[assumed_UTF8,flag,ME]=UTF8_to_unicode_internal(UTF8,return_on_error);
if strcmp(flag,'success')
    isUTF8=true;
    unicode=assumed_UTF8;
elseif strcmp(flag,'error')
    isUTF8=false;
    if return_on_error
        error_(print_to,ME)
    end
    unicode=UTF8; % Return input unchanged (apart from casting to uint32).
end
end
function [UTF8,flag,ME]=UTF8_to_unicode_internal(UTF8,return_on_error)
flag='success';
ME=struct('identifier','HJW:UTF8_to_unicode:notUTF8','message','Input is not UTF-8.');

persistent isOctave,if isempty(isOctave),isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;end

if any(UTF8>255)
    flag='error';
    if return_on_error,return,end
elseif all(UTF8<128)
    return
end

for bytes=4:-1:2
    val=bin2dec([repmat('1',1,bytes) repmat('0',1,8-bytes)]);
    multibyte=UTF8>=val & UTF8<256;%Exclude the already converted chars.
    if any(multibyte)
        multibyte=find(multibyte);multibyte=multibyte(:).';
        if numel(UTF8)<(max(multibyte)+bytes-1)
            flag='error';
            if return_on_error,return,end
            multibyte( (multibyte+bytes-1)>numel(UTF8) )=[];
        end
        if ~isempty(multibyte)
            idx=bsxfun_plus(multibyte , (0:(bytes-1)).' );
            idx=idx.';
            multibyte=UTF8(idx);
        end
    else
        multibyte=[];
    end
    header_bits=[repmat('1',1,bytes-1) repmat('10',1,bytes)];
    header_locs=unique([1:(bytes+1) 1:8:(8*bytes) 2:8:(8*bytes)]);
    if numel(multibyte)>0
        multibyte=unique(multibyte,'rows');
        S2=mat2cell(multibyte,ones(size(multibyte,1),1),bytes);
        for n=1:numel(S2)
            bin=dec2bin(double(S2{n}))';
            % To view the binary data, you can use this: bin=bin(:)';
            % Remove binary header (3 byte example):
            % 1110xxxx10xxxxxx10xxxxxx
            %     xxxx  xxxxxx  xxxxxx
            if ~strcmp(header_bits,bin(header_locs))
                % Check if the byte headers match the UTF-8 standard.
                flag='error';
                if return_on_error,return,end
                continue %leave unencoded
            end
            bin(header_locs)='';
            if ~isOctave
                S3=uint32(bin2dec(bin  ));
            else
                S3=uint32(bin2dec(bin.'));% Octave needs an extra transpose.
            end
            % Perform actual replacement.
            UTF8=PatternReplace(UTF8,S2{n},S3);
        end
    end
end
end
function [opts,ME]=validate_print_to__options(opts_in,ME)
% If any input is invalid, this returns an empty array and sets ME.message.
%
% Input struct:
% options.print_to_con=true;   % or false
% options.print_to_fid=fid;    % or []
% options.print_to_obj=h_obj;  % or []
% options.print_to_fcn=struct; % or []
%
% Output struct:
% options.fid
% options.obj
% options.fcn.h
% options.fcn.data
% options.boolean.con
% options.boolean.fid
% options.boolean.obj
% options.boolean.fcn

% Set defaults.
if nargin<2,ME=struct;end
if ~isfield(opts_in,'print_to_con'),opts_in.print_to_con=[];end
if ~isfield(opts_in,'print_to_fid'),opts_in.print_to_fid=[];end
if ~isfield(opts_in,'print_to_obj'),opts_in.print_to_obj=[];end
if ~isfield(opts_in,'print_to_fcn'),opts_in.print_to_fcn=[];end
print_to_con_default=true; % Unless a valid fid, obj, or fcn is specified.

% Initialize output.
opts=struct;

% Parse the fid. We can use ftell to determine if fprintf is going to fail.
item=opts_in.print_to_fid;
if isempty(item)
    opts.boolean.fid=false;
    opts.fid=[];
else
    print_to_con_default=false;
    opts.boolean.fid=true;
    opts.fid=item;
    for n=1:numel(item)
        try position=ftell(item(n));catch,position=-1;end
        if item(n)~=1 && position==-1
            ME.message=['Invalid print_to_fid parameter:',char(10),...
                'should be a valid file identifier or 1.']; %#ok<CHARTEN>
            opts=[];return
        end
    end
end

% Parse the object handle. Retrieving from multiple objects at once works, but writing that output
% back to multiple objects doesn't work if Strings are dissimilar.
item=opts_in.print_to_obj;
if isempty(item)
    opts.boolean.obj=false;
    opts.obj=[];
else
    print_to_con_default=false;
    opts.boolean.obj=true;
    opts.obj=item;
    for n=1:numel(item)
        try
            txt=get(item(n),'String'    ); % See if this triggers an error.
            set(    item(n),'String','' ); % Test if property is writable.
            set(    item(n),'String',txt); % Restore original content.
        catch
            ME.message=['Invalid print_to_obj parameter:',char(10),...
                'should be a handle to an object with a writeable String property.']; %#ok<CHARTEN>
            opts=[];return
        end
    end
end

% Parse the function handles.
item=opts_in.print_to_fcn;
if isempty(item)
    opts.boolean.fcn=false;
    opts.fcn=[];
else
    print_to_con_default=false;
    try
        for n=1:numel(item)
            if ~ismember(class(item(n).h),{'function_handle','inline'}) || numel(item(n).h)~=1
                error('trigger error')
            end
        end
    catch
        ME.message=['Invalid print_to_fcn parameter:',char(10),...
            'should be a struct with the h field containing a function handle,',char(10),...
            'anonymous function or inline function.']; %#ok<CHARTEN>
        opts=[];return
    end
end

% Parse the logical that determines if a warning will be printed to the command window.
% This is true by default, unless an fid, obj, or fcn is specified.
item=opts_in.print_to_con;
if isempty(item)
    opts.boolean.con=print_to_con_default;
else
    [passed,opts.boolean.con]=test_if_scalar_logical(item);
    if ~passed
        ME.message=['Invalid print_to_con parameter:',char(10),...
            'should be a scalar logical.']; %#ok<CHARTEN>
        opts=[];return
    end
end
end
function varargout=var2str(varargin)
%Analogous to func2str, return the variable names as char arrays, as detected by inputname.
%This returns an error for invalid inputs and if nargin~=max(1,nargout).
%
%You can use comma separated lists to create a cell array:
% out=cell(1,2);
% foo=1;bar=2;
% [out{:}]=var2str(foo,bar);

%One-line alternative: function out=var2str(varargin),out=inputname(1);end
err_flag= nargin~=max(1,nargout) ;
if ~err_flag
    varargout=cell(nargin,1);
    for n=1:nargin
        try varargout{n}=inputname(n);catch,varargout{n}='';end
        if isempty(varargout{n}),err_flag=true;break,end
    end
end
if err_flag
    error('Invalid input and/or output.')
end
end
function warning_(options,varargin)
%Print a warning to the command window, a file and/or the String property of an object
% The lastwarn state will be set if the warning isn't thrown with warning().
% The printed call trace omits this function, but the warning() call does not.
%
% You can also provide a struct (scalar or array) with two fields: 'h' with a function handle, and
% 'data' with arbitrary data passed as third input. These functions will be run with 'warning' as
% first input. The second input is a struct with identifier, message, and stack as fields. This
% function will be run with feval (meaning the function handles can be replaced with inline
% functions or anonymous functions).
%
% The intention is to allow replacement of most warning(___) call with warning_(options,___). This
% does not apply to calls that query or set the warning state.
%
%options.boolean.con: if true print warning to command window with warning()
%options.fid:         file identifier for fprintf (array input will be indexed)
%options.boolean.fid: if true print warning to file (options.fid)
%options.obj:         handle to object with String property (array input will be indexed)
%options.boolean.obj: if true print warning to object (options.obj)
%options.fcn          struct (array input will be indexed)
%options.fcn.h:       handle of function to be run
%options.fcn.data:    data passed as third input to function to be run (optional)
%options.boolean.fnc: if true the function(s) will be run
%
% Syntax:
%   warning_(options,msg)
%   warning_(options,msg,A1,...,An)
%   warning_(options,id,msg)
%   warning_(options,id,msg,A1,...,An)
%   warning_(options,ME)               %rethrow error as warning
%
%examples options struct:
%  % Write to a log file:
%  opts=struct;opts.fid=fopen('log.txt','wt');
%  % Display to a status window and bypass the command window:
%  opts=struct;opts.boolean.con=false;opts.obj=uicontrol_object_handle;
%  % Write to 2 log files:
%  opts=struct;opts.fid=[fopen('log2.txt','wt') fopen('log.txt','wt')];

persistent this_fun
if isempty(this_fun),this_fun=func2str(@warning_);end

[id,msg,stack,trace,no_op]=parse_warning_error_redirect_inputs( varargin{:});
if no_op,return,end % Don't waste time parsing the options in case of a no-op.
ME=struct('identifier',id,'message',msg,'stack',stack);

% Parse options struct, allowing an empty input to revert to default.
if isempty(options),options=validate_print_to__options(struct);end
options                   =parse_warning_error_redirect_options(  options  );

% Check if the warning is turned off and exit the function if this is the case.
w=warning;if any(ismember({w(ismember({w.identifier},{id,'all'})).state},'off')),return,end

if options.boolean.con
    if ~isempty(id),warning(id,'%s',msg),else,warning(msg), end
else
    if ~isempty(id),lastwarn(msg,id);    else,lastwarn(msg),end
end

if options.boolean.obj
    msg_=msg;while msg_(end)==10,msg_(end)=[];end%Crop trailing newline.
    if any(msg_==10)  % Parse to cellstr and prepend warning.
        msg_=char2cellstr(['Warning: ' msg_]);
    else              % Only prepend warning.
        msg_=['Warning: ' msg_];
    end
    set(options.obj,'String',msg_)
    for OBJ=options.obj(:).'
        try set(OBJ,'String',msg_);catch,end
    end
end

if options.boolean.fid
    for FID=options.fid(:).'
        try fprintf(FID,'Warning: %s\n%s',msg,trace);catch,end
    end
end

if options.boolean.fcn
    if ismember(this_fun,{stack.name})
        % To prevent an infinite loop, trigger an error.
        error('prevent recursion')
    end
    for FCN=options.fcn(:).'
        if isfield(FCN,'data')
            try feval(FCN.h,'warning',ME,FCN.data);catch,end
        else
            try feval(FCN.h,'warning',ME);catch,end
        end
    end
end
end
function outfilename=WBM(filename,url_part,varargin)
%This functions acts as an API for the Wayback Machine (web.archive.org)
%
% With this function you can download captures to the internet archive that matches a date pattern.
% If the current time matches the pattern and there is no valid capture, a capture will be
% generated. The WBM time stamps are in UTC, so a switch allows you to provide the date-time
% pattern in local time, whose upper and lower bounds will be shifted to UTC internally.
%
% This code enables you to use a specific web page in your data processing, without the need to
% check if the page has changed its structure or is not available at all.
% You can also redirect all outputs (errors only partially) to a file or a graphics object, so you
% can more easily use this function in a GUI or allow it to write to a log file.
%
% Usage instruction about the syntax of the WBM interface are derived from a Wikipedia help page:
% https://en.wikipedia.org/wiki/Help:Using_the_Wayback_Machine
% Fuzzy date matching behavior is based on https://archive.org/web/web-advancedsearch.php
%
% Syntax:
%   outfilename=WBM(filename,url_part)
%   outfilename=WBM(___,Name,Value)
%   outfilename=WBM(___,options)
%
% Input/output arguments:
% outfilename:
%   Full path of the output file, the variable is empty if the download failed.
% filename:
%   The target filename in any format that websave (or urlwrite) accepts. If this file already
%   exists, it will be overwritten in most cases.
% url_part:
%   This URL will be searched for on the WBM. The URL might be changed (e.g. ':80' is often added).
% options:
%   A struct with Name,Value parameters. Missing parameters are filled with the defaults listed
%   below. Using incomplete parameter names or incorrect capitalization is allowed, as long as
%   there is a unique match.
%   Parameters related to warning/error redirection will be parsed first.
%
% Name,Value parameters:
%   date_part:
%      A string with the date of the capture. It must be in the yyyymmddHHMMSS format, but doesn't
%      have to be complete. Note that this is represented in UTC. If incomplete, the Wayback
%      Machine will return a capture that is as close to the midpoint of the matching range as
%      possible. So for date_part='2' the range is 2000-01-01 00:00 to 2999-12-31 23:59:59, meaning
%      the WBM will attempt to return the capture closest to 2499-12-31 23:59:59. [default='2';]
%   UseLocalTime:
%      A scalar logical. Interpret the date_part in local time instead of UTC. This has the
%      practical effect of the upper and lower bounds of the matching date being shifted by the
%      timezone offset. [default=false;]
%   tries:
%      A 1x3 vector. The first value is the total number of times an attempt to load the page is
%      made, the second value is the number of save attempts and the last value is the number of
%      timeouts allowed. [default=[5 4 4];]
%   verbose:
%      A scalar denoting the verbosity. Level 0 will hide all errors that are caught. Level 1 will
%      enable only warnings about the internet connection being down. Level 2 includes errors NOT
%      matching the usual pattern as well and level 3 includes all other errors that get rethrown
%      as warning. Level 4 adds a warning on retrying to retrieve the timestamp (this will always
%      happen for pages without a capture, and may happen every now and then for no apparent
%      reason).[default=3;]
%      Octave uses libcurl, making error catching is bit more difficult. This will result in more
%      HTML errors being rethrown as warnings under Octave than Matlab. There is no check in place
%      to restrict this to level 0-4 (even inf is allowed), so you can anticipate on higher levels
%      in future versions of this function.
%   if_UTC_failed:
%      This is a char array with the intended behavior for when this function is unable to
%      determine the UTC. The options are 'error', 'warn_0', 'warn_1', 'warn_2', 'warn_3', and
%      'ignore'. For the options starting with warn_, a warning will be triggered if the 'verbose'
%      parameter is set to this level or higher (so 'warn_0' will trigger a warning if 'verbose' is
%      set to 0).
%      If this parameter is not set to 'error', the valid time range is expanded by -12 and +14
%      hours to account for all possible time zones, and the midpoint is shifted accordingly.
%      [default='warn_3']
%   m_date_r:
%      A string describing the response to the date missing in the downloaded web page. Usually,
%      either the top bar will be present (which contains links), or the page itself will contain
%      links, so this situation may indicate a problem with the save to the WBM. Allowed values are
%      'ignore', 'warning' and 'error'. Be aware that non-page content (such as images) will set
%      off this response. Flags other than '*' will also set off this response.
%      [default='warning';] if flags~='*' then [default='ignore']
%   response:
%      The response variable is a cell array, where each row encodes one sequence of HMTL errors
%      and the appropriate next action. The syntax of each row is as follows:
%      #1 If there is a sequence of failure that fit the first cell,
%      #2 and the HTML error codes of the sequence are equal to the second cell,
%      #3 then respond as per the third cell.
%      The sequence of failures are encoded like this:
%      t1: failed attempt to load, t2: failed attempt to save, tx: either failed to load, or failed
%      to save.
%      The error code list must be HTML status codes. The Matlab timeout error is encoded with 4080
%      (analogous to the HTTP 408 timeout error code). The  error is extracted from the identifier,
%      which is not always possible, especially in the case of Octave.
%      The response in the third cell is either 'load', 'save', 'exit', or 'pause_retry'.
%      Load and save set the preferred type. If a response is not allowed (i.e. if the
%      corresponding element of 'tries' is 0), the other response (save or load) is tried, until
%      sum(tries(1:2))==0. If the response is set to exit, or there is still no successful download
%      after tries has been exhausted, the output file will be deleted and the script will exit.
%      The pause_retry is intended for use with an error 429. See the err429 parameter for more
%      options. [default={'tx',404,'load';'txtx',[404 404],'save';'tx',403,'save';
%      't2t2',[403 403],'exit';'tx',429,'pause_retry'};]
%   err429:
%      Sometimes the webserver will return a 429 status code. This should trigger a waiting period
%      of a few seconds. This parameter controls the behavior of this function in case of a 429
%      status code. It is a struct with the following fields.
%      The CountsAsTry field (logical) describes if the attempt should decrease the tries counter.
%      The TimeToWait field (double) contains the time in seconds to wait before retrying.
%      The PrintAtVerbosityLevel field (double) contains the verbosity level at which a text should
%      be printed, showing the user the function did not freeze.
%      [default=struct('CountsAsTry',false,'TimeToWait',15, 'PrintAtVerbosityLevel',3);]
%   ignore:
%      The ignore variable is vector with the same type of error codes as in the response variable.
%      Ignored errors will only be ignored for the purposes of the response, they will not prevent
%      the tries vector from decreasing. [default=4080;]
%   flag:
%      The flags can be used to specify an explicit version of the archived page. The options are
%      'id' (identical), 'js' (Javascript), 'cs' (CSS), 'im' (image), 'fw'/'if' (iFrame), or *
%      (explicitly expand date, shows calendar in browser mode). With the 'id' flag the page is
%      show as captured (i.e. without the WBM banner, making it ideal for e.g. exe files). With the
%      'id' and '*' flags the date check will fail, so the missing date response (m_date_r) will be
%      invoked. For the 'im' flag you can circumvent this by first loading in the normal mode
%      ('*'), and then extracting the image link from that page. That way you can enforce a date
%      pattern and still get the image. The Wikipedia page suggest that a flag syntax requires a
%      full date, but this seems not to be the case, as the date can still auto-expand.
%      [default='*';]
%   waittime:
%      This value controls the maximum time that is spent waiting on the internet connection for
%      each call of this function. This does not include the time waiting as a result of a 429
%      error. The input must be convertible to a scalar double. This is the time in seconds.
%      NB: Setting this to inf will cause an infite loop if the internet connection is lost.
%      [default=60]
%   timeout:
%      This value is the allowed timeout in seconds. It is ignored if it isn't supported. The input
%      must be convertible to a scalar double. [default=10]
%   print_to_con:
%      An attempt is made to also use this parameter for warnings or errors during input parsing.
%      A logical that controls whether warnings and other output will be printed to the command
%      window. Errors can't be turned off. [default=true;] if either print_to_fid, print_to_obj, or
%      print_to_fcn is specified then [default=false]
%   print_to_fid:
%      An attempt is made to also use this parameter for warnings or errors during input parsing.
%      The file identifier where console output will be printed. Errors and warnings will be
%      printed including the call stack. You can provide the fid for the command window (fid=1) to
%      print warnings as text. Errors will be printed to the specified file before being actually
%      thrown. [default=[];]
%      If print_to_fid, print_to_obj, and print_to_fcn are all empty, this will have the effect of
%      suppressing every output except errors. 
%      Array inputs are allowed.
%   print_to_obj:
%      An attempt is made to also use this parameter for warnings or errors during input parsing.
%      The handle to an object with a String property, e.g. an edit field in a GUI where console
%      output will be printed. Messages with newline characters (ignoring trailing newlines) will
%      be returned as a cell array. This includes warnings and errors, which will be printed
%      without the call stack. Errors will be written to the object before the error is actually
%      thrown. [default=[];]
%      If print_to_fid, print_to_obj, and print_to_fcn are all empty, this will have the effect of
%      suppressing every output except errors.
%      Array inputs are allowed.
%   print_to_fcn:
%      An attempt is made to also use this parameter for warnings or errors during input parsing.
%      A struct with a function handle, anonymous function or inline function in the 'h' field and
%      optionally additional data in the 'data' field. The function should accept three inputs: a
%      char array (either 'warning' or 'error'), a struct with the message, id, and stack, and the
%      optional additional data. The function(s) will be run before the error is actually thrown.
%      [default=[];]
%      If print_to_fid, print_to_obj, and print_to_fcn are all empty, this will have the effect of
%      suppressing every output except errors.
%      Array inputs are allowed.
%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 3.1.0                                                         |%
%|  Date:    2022-01-22                                                    |%
%|  Author:  H.J. Wisselink                                                |%
%|  Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 ) |%
%|  Email = 'h_j_wisselink*alumnus_utwente_nl';                            |%
%|  Real_email = regexprep(Email,{'*','_'},{'@','.'})                      |%
%|                                                                         |%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%
% Tested on several versions of Matlab (ML 6.5 and onward) and Octave (4.4.1 and onward), and on
% multiple operating systems (Windows/Ubuntu/MacOS). For the full test matrix, see the HTML doc.
% Compatibility considerations:
% - HTML error codes are harder to catch on Octave. Depending on the selected verbosity level that
%   means the number of warnings will be larger.
% - The duration of a timeout can only be set with websave. This means that for larger files or
%   less stable internet connections, a timeout error will be more likely when using older releases
%   or Octave.

if nargin<2
    error('HJW:WBM:nargin','Incorrect number of input argument.')
end
if ~(nargout==0 || nargout==1) %might trigger 'MATLAB:TooManyOutputs' instead
    error('HJW:WBM:nargout','Incorrect number of output argument.')
end
[success,opts,ME]=WBM_parse_inputs(filename,url_part,varargin{:});
if ~success
    %If the parsing of print_to failed (which is tried first), the default will be used.
    error_(opts.print_to,ME)
else
    [        date_part,     tries,     response,     ignore,     verbose,     UseURLwrite ...
        ,     flag,     err429,     UTC_fail_response,     timeout,     print_to]=...
        deal(...
        opts.date_part,opts.tries,opts.response,opts.ignore,opts.verbose,opts.UseURLwrite ...
        ,opts.flag,opts.err429,opts.UTC_fail_response,opts.timeout,opts.print_to);
    SavesAllowed=tries(2)>0;
    %Cap the time we wait for an internet connection to avoid infinite loops (time is in seconds).
    waittime=struct('cap',opts.waittime,'total',0);
end
if ~TestFolderWritePermission(fileparts(filename))
    error_(print_to,'HJW:WBM:NoWriteFolder',...
        'The target folder doesn''t exist or Matlab doesn''t have write permission for it.')
end
if ~WBM_UTC_test % No UTC determination is available.
    if strcmp(UTC_fail_response,'error')
        error_(print_to,'HJW:WBM:UTC_missing_error','Retrieval of UTC failed.')
    elseif strcmp(UTC_fail_response,'warn')
        warning_(print_to,'HJW:WBM:UTC_missing_warning',...
            'Retrieval of UTC failed, continuing with wider date/time range.')
    end
end

%Order responses based on pattern length to match. There is no need to do this in the parser, since
%this is fast anyway, and response_lengths is needed in the loop below as well.
response_lengths=cellfun('length',response(:,2));
[response_lengths,order]=sort(response_lengths);
order=order(end:-1:1);%sort(__,'descend'); is not supported in ML6.5
response=response(order,:);
response_lengths=response_lengths(end:-1:1);

prefer_type=1;%prefer loading
success=false;%start loop
response_list_vector=[];%initialize response list
type_list=[];%initialize type list
connection_down_wait_factor=0;%initialize
while ~success && ...         %no successful download yet?
        sum(tries(1:2))>0 ... %any save or load tries left?
        && tries(3)>=0        %timeout limit reached?
    if tries(prefer_type)<=0%no tries left for the preferred type
        prefer_type=3-prefer_type;%switch 1 to 2 and 2 to 1
    end
    type=prefer_type;
    try ME=[]; %#ok<NASGU>
        if type==1 %load
            SaveAttempt=false;
            tries(type)=tries(type)-1;
            
            [status,waittime]=confirm_capture_is_available(url_part,date_part,waittime,opts);
            if status<200 || status>=300
                %No capture is available for this URL. Report this to the user (if the verbosity is
                %set high enough). Then return to the start of the loop to try saving (if there are
                %save attempts left).
                if verbose>=3
                    txt=sprintf('No capture found for this URL. (download of %s)',filename);
                    if print_to.boolean.con
                        fprintf('%s\n',txt)
                    end
                    if print_to.boolean.fid
                        for FID=print_to.fid(:).',fprintf(FID,'%s\n',txt);end
                    end
                    if print_to.boolean.obj
                        set(print_to.obj,'String',txt)
                    end
                    drawnow
                end
                %Try saving or quit.
                prefer_type=2;
                if tries(2)>0,continue
                else         ,break
                end
            end
            
            %If the execution reaches this point, there is a capture available, 
            if UseURLwrite
                outfilename=urlwrite(...
                    ['http://web.archive.org/web/' date_part flag '_/' url_part],...
                    filename);%#ok<URLWR>
                outfilename=check_filename(filename,outfilename);
            else
                outfilename=websave(filename,...
                    ['https://web.archive.org/web/' date_part flag '_/' url_part],...
                    weboptions('Timeout',timeout));
            end
        elseif type==2 %save
            SaveAttempt=true;
            tries(type)=tries(type)-1;
            if UseURLwrite
                outfilename=urlwrite(...
                    ['http://web.archive.org/save/' url_part],...
                    filename); %#ok<URLWR>
                outfilename=check_filename(filename,outfilename);
            else
                outfilename=websave(filename,...
                    ['https://web.archive.org/save/' url_part],...
                    weboptions('Timeout',timeout));
            end
        end
        success=true;connection_down_wait_factor=0;
        if SavesAllowed && ~check_date(outfilename,opts,SaveAttempt)
            %Incorrect date or live page loaded, so try saving.
            success=false;prefer_type=2;
        end
    catch ME;if isempty(ME),ME=lasterror;end%#ok<LERR>
        success=false;
        if ~isnetavl
            %If the connection is down, retry in intervals
            while ~isnetavl
                if waittime.total>waittime.cap
                    %Total wait time exceeded, return the error condition.
                    if verbose>=1
                        warning_(print_to,'Maximum waiting time for internet connection exceeded.')
                    end
                    tries=zeros(size(tries));continue
                end
                curr_time=datestr(now,'HH:MM:SS');
                if verbose>=1
                    warning_(print_to,'Internet connection down, retrying in %d seconds (@%s).',...
                        2^connection_down_wait_factor,curr_time)
                end
                pause(2^connection_down_wait_factor)
                waittime.total=waittime.total+2^connection_down_wait_factor;
                %Increment, but cap to a reasonable interval.
                connection_down_wait_factor=min(1+connection_down_wait_factor,6);
            end
            %Skip the rest of the error processing and retry without reducing points.
            continue
        end
        connection_down_wait_factor=0;
        ME_id=ME.identifier;
        ME_id=strrep(ME_id,':urlwrite:',':webservices:');
        if strcmp(ME_id,'MATLAB:webservices:Timeout')
            code=4080;
            tries(3)=tries(3)-1;
        else
            %equivalent to raw_code=textscan(ME_id,'MATLAB:webservices:HTTP%dStatusCodeError');
            raw_code=strrep(ME_id,'MATLAB:webservices:HTTP','');
            raw_code=strrep(raw_code,'StatusCodeError','');
            raw_code=str2double(raw_code);
            if isnan(raw_code)
                %Some other error occurred, set a code and throw a warning. As Octave does not
                %report an HTML error code, this will happen almost every error. To reduce command
                %window clutter, consider lowering the verbosity level.
                code=-1;
                if verbose>=2
                    warning_(print_to,ME)
                end
            else
                %Octave doesn't really returns a identifier for urlwrite, nor do very old releases
                %of Matlab.
                switch ME.message
                    case 'urlwrite: Couldn''t resolve host name'
                        code=404;
                    case ['urlwrite: Peer certificate cannot be ',...
                            'authenticated with given CA certificates']
                        %It's not really a 403, but the result in this context is similar.
                        code=403;
                    otherwise
                        code=raw_code;
                end
            end
        end
        
        if verbose>=3
            txt=sprintf('Error %d tries(%d,%d,%d) (download of %s).',...
                double(code),tries(1),tries(2),tries(3),filename);
            if print_to.boolean.con
                fprintf('%s\n',txt)
            end
            if print_to.boolean.fid
                for FID=print_to.fid(:).',fprintf(FID,'%s\n',txt);end
            end
            if print_to.boolean.obj
                set(print_to.obj,'String',txt)
            end
            drawnow
        end
        if ~any(code==ignore)
            response_list_vector(end+1)=code; %#ok<AGROW>
            type_list(end+1)=type; %#ok<AGROW>
            for n_response_pattern=1:size(response,1)
                if length(response_list_vector) < response_lengths(n_response_pattern)
                    %Not enough failed attempts (yet) to match against the current pattern.
                    continue
                end
                last_part_of_response_list=response_list_vector(...
                    (end-response_lengths(n_response_pattern)+1):end);
                last_part_of_type_list=type_list(...
                    (end-response_lengths(n_response_pattern)+1):end);
                
                %Compare the last types to the type patterns.
                temp_type_pattern=response{n_response_pattern,1}(2:2:end);
                temp_type_pattern=strrep(temp_type_pattern,'x',num2str(type));
                type_fits=strcmp(temp_type_pattern,sprintf('%d',last_part_of_type_list));
                if isequal(...
                        response{n_response_pattern,2},...
                        last_part_of_response_list)...
                        && type_fits
                    %If the last part of the response list matches with the response pattern in the
                    %current element of 'response', set prefer_type to 1 for load, and to 2 for
                    %save.
                    switch response{n_response_pattern,3}
                        %otherwise will not occur: should be caught in the input parser
                        case 'load'
                            prefer_type=1;
                        case 'save'
                            prefer_type=2;
                        case 'exit'
                            %Cause a break in the while loop.
                            tries=[0 0 -1];
                        case 'pause_retry'
                            if ~err429.CountsAsTry
                                %Increment the counter, which has the effect of not counting this
                                %as a try.
                                tries(prefer_type)=tries(prefer_type)+1;
                            end
                            if verbose>=err429.PrintAtVerbosityLevel
                                N=10;
                                s='Waiting a while until the server won''t block us anymore';
                                if print_to.boolean.con
                                    fprintf(s);drawnow
                                end
                                if print_to.boolean.fid
                                    for FID=print_to.fid(:).',fprintf(FID,s);end,drawnow
                                end
                                for n=1:N
                                    pause(err429.TimeToWait/N)
                                    if print_to.boolean.con
                                        fprintf('.');drawnow
                                    end
                                    if print_to.boolean.fid
                                        for FID=print_to.fid(:).',fprintf(FID,'.');end,drawnow
                                    end
                                    if print_to.boolean.obj
                                        s=[s '.']; %#ok<AGROW>
                                        set(print_to.obj,s);drawnow
                                    end
                                end
                                if print_to.boolean.con
                                    fprintf('\nContinuing\n');drawnow
                                end
                                if print_to.boolean.fid
                                    for FID=print_to.fid(:).',fprintf(FID,'\nContinuing\n');end
                                    drawnow
                                end
                                if print_to.boolean.obj
                                    set(print_to.obj,'Continuing');drawnow
                                end
                            else
                                pause(err429.TimeToWait)
                            end
                    end
                    break
                end
            end
        end
    end
end

if ~success || ( ~SavesAllowed && ~check_date(outfilename,opts,SaveAttempt) )
    %If saving isn't allowed and the date doesn't match the date_part, or no successful download
    %was reached within the allowed tries, delete the output file (as it will be either the
    %incorrect date, or 0 bytes).
    if exist(filename,'file'),try delete(filename);catch,end,end
    outfilename=[];
end
filename2=[filename '.html'];
if exist(filename2,'file')
    a=dir(filename2);
    if numel(a)==1 && a.bytes==0 && abs(datenum(a.date)-now)<=(1/24) %newly created
        try delete(filename2);catch,end %assume a 0 byte is never correct (although it might be)
    end
end
if nargout==0
    clear(var2str(outfilename));
end
end
function [status,waittime]=confirm_capture_is_available(url,date,waittime,opts)
%Retrieve the time stamp closest to the requested date and compare it to date_part.
[t,status,waittime]=WBM_retrieve_timestamp(url,date,opts,waittime);
if isempty(t)
    %Sometimes this will just fail for some reason. We should probably retry only once to avoid an
    %unnecessary save. See WBM_retrieve_timestamp for the conditions under which this might occur.
    if opts.verbose>3 %Print a warning to alert the user.
        msg=sprintf('Retrying to retrieve the timestamp for %s',url);
        warning_(opts.print_to,'HJW:WBM:RetryTimestamp',msg)
    end
    [t,status,waittime]=WBM_retrieve_timestamp(url,date,opts,waittime);
end
if status<200 || status>=300
    %There might be a timestamp available, but since the status is not 2xx urlwrite and websave
    %will probably error anyway.
    return
end
date_as_double=str2double(t);
if opts.date_bounds.double(1)>date_as_double || date_as_double>opts.date_bounds.double(2)
    %As a capture with an appropriate timestamp could not be found, a 404 status code seems
    %appropriate. This should normally trigger a save.
    status=404;
    return
end
end
function atomTime=WBM_getUTC_local
persistent pref_order % Some checks take a relatively long time, so store this in a persistent.
if isempty(pref_order)
    tf=CheckMexCompilerExistence;
    [ignore,status]=GetWritableFolder('ErrorOnNotFound',false); %#ok<ASGLU>
    if status==0,status=inf;end
    if tf && status<3
        % There is a compiler, and there is a writable folder other than the pwd, so the c
        % implementation should be tried before the web implementation.
        % Since the c implementation is probably faster than the OS call, try that first.
        pref_order=[1 3 2];
    else
        % Either a compiler is missing, or only the pwd is writable. Only use c if web fails.
        pref_order=[3 2 1];
    end
end
for n=1:numel(pref_order)
    atomTime=max([0 getUTC(pref_order(n))]); % Ensure it is 0 if getUTC fails
    if atomTime~=0,break,end
end
end
function [bounds,midpoint]=WBM_parse_date_to_range(date_part)
%Match the partial date to the range of allowed dates.
%This throws an error if the date_part doesn't fit a valid date pattern.
shift=[...
    1000, 0, 0, 0, 0, 0;...
    100 , 0, 0, 0, 0, 0;...
    10  , 0, 0, 0, 0, 0;...
    1   , 0, 0, 0, 0, 0;...
    0   ,10, 0, 0, 0, 0;...
    0   , 1, 0, 0, 0, 0;...
    0   , 0,10, 0, 0, 0;...
    0   , 0, 1, 0, 0, 0;...
    0   , 0, 0,10, 0, 0;...
    0   , 0, 0, 1, 0, 0;...
    0   , 0, 0, 0,10, 0;...
    0   , 0, 0, 0, 1, 0;...
    0   , 0, 0, 0, 0,10;...
    0   , 0, 0, 0, 0, 1];

date=char(zeros(1,14)+'0');
date(1:numel(date_part))=date_part;
date={date(1:4),date(5:6),date(7:8),...  %date
    date(9:10),date(11:12),date(13:14)}; %time
date=str2double(date);date(1:3)=max(date(1:3),1);

lower_bound=num2cell(date);
lower_bound=datenum(lower_bound{:});

test=datestr(lower_bound,'yyyymmddTHHMMSS');
test(9)='';
if ~strcmp(test(1:numel(date_part)),date_part)
    error('incorrect date_part format')
end

upper_bound=date+shift(numel(date_part),:);
d=[31 28 31 30 31 30 31 31 30 31 30 31];
y=upper_bound(1);if rem(y,4)==0 && (rem(y,100)~=0 || rem(y,400)==0),d(2) = 29;end
d=d(min(12,upper_bound(2)));
overflow=find(upper_bound>[inf 12 d 23 59 59])+1;
if ~isempty(overflow),upper_bound(overflow:end)=inf;else,upper_bound(end)=-1;end
upper_bound=min(upper_bound,[inf 12 d 23 59 59]);
upper_bound=num2cell(upper_bound);
upper_bound=datenum(upper_bound{:});


%If the UTC time can't be determined, the valid time should be expanded by -12 and +14 hours.
%If an error or warning should be thrown, this will happen in the main function.
if ~WBM_UTC_test % Expand by -12 and +14 hours.
    lower_bound=lower_bound+datenum(0,0,0,-12,0,0);
    upper_bound=upper_bound+datenum(0,0,0,+14,0,0);
end

midpoint=(lower_bound+upper_bound)/2;
midpoint=datestr(midpoint,'yyyymmddTHHMMSS');
midpoint(9)='';
bounds=[lower_bound upper_bound];
end
function [success,opts,ME]=WBM_parse_inputs(filename,url,varargin)
%The print_to variable will be parsed first. If the parsing of print_to fails, an empty double will
%be returned.
success=false;opts=struct;
ME=struct('identifier','','message','');

%Attempt to match the struct to the available options. This will return a struct with the same
%fields as the default option struct.
try ME_=[];[opts,replaced]=parse_NameValue(WBM_parse_inputs_defaults,varargin{:}); %#ok<NASGU>
catch ME_;if isempty(ME_),ME_=lasterror;end,ME=ME_;return,end %#ok<LERR>
%Add the required inputs to the struct.
opts.filename=filename;opts.url=url;

%Attempt to parse the error redirection options (this generates ME.message on fail).
[opts.print_to,ME]=validate_print_to__options(opts,ME);
if isempty(opts.print_to)
    ME.identifier='HJW:WBM:Incorrect_print_to_option';return
end

%Test the required inputs.
[valid,opts.filename]=filename_is_valid(opts.filename);
if ~valid
    ME.message='The first input (filename) is not char and/or empty.';
    ME.identifier='HJW:WBM:incorrect_input_filename';
    return
end
if ~ischar(opts.url) || numel(opts.url)==0
    ME.message='The second input (url) is not char and/or empty.';
    ME.identifier='HJW:WBM:incorrect_input_url';
    return
end

if numel(replaced)==0,success=true;ME=[];return,end % no default values were changed

%Check optional inputs
default=WBM_parse_inputs_defaults;
for k=1:numel(replaced)
    item=opts.(replaced{k});
    switch replaced{k}
        case {'print_to_con','print_to_fid','print_to_obj','print_to_fcn'}
            %Already checked
        case 'date_part'
            valid_date=true;
            if  ~ischar(item) || numel(item)==0 || numel(item)>14 || any(item<48 & item>57)
                valid_date=false;
            end
            %Parse the date to bounds. (editing note: this code is also in the default constructor)
            try [bounds,midpoint]=WBM_parse_date_to_range(item);catch,valid_date=false;end
            if ~valid_date
                ME.message='The value of options.date_part is empty or not a valid numeric char.';
                return
            end
            opts.date_part=midpoint;
            opts.date_bounds.datenum=bounds;
            %Note: options.date_bounds.double is set later
        case 'tries'
            if ~isnumeric(item) || numel(item)~=3 || any(isnan(item))
                ME.message=['The value of options.tries has an incorrect format.',char(10),...
                    'The value should be a numeric vector with 3 integer elements.'];%#ok<CHARTEN>
                return
            end
        case 'response'
            if WBM_parse_inputs__validate_response_format(item)
                ME.message='The value of options.response has an incorrect format.';
                return
            end
        case 'ignore'
            if ~isnumeric(item) || numel(item)==0 || any(isnan(item))
                ME.message=['The value of options.ignore has an incorrect format.',char(10),...
                    'The value should be a numeric vector with HTML error codes.'];%#ok<CHARTEN>
                return
            end
        case 'verbose'
            if ~isnumeric(item) || numel(item)~=1 || double(item)~=round(double(item))
                %The integer test could cause unexpected behavior due to float rounding, but in
                %fact an error is preferred here.
                ME.message='The value of options.verbose is not an integer scalar.';
                return
            end
        case 'm_date_r'
            if ~ischar(item) || numel(item)==0
                ME.message='Options.m_date_r should be ''ignore'', ''warning'', or ''error''.';
                return
            end
            switch lower(item)
                case 'ignore'
                    item=0;
                case 'warning'
                    item=1;
                case 'error'
                    item=2;
                otherwise
                    ME.message='Options.m_date_r should be ''ignore'', ''warning'', or ''error''.';
                    return
            end
            opts.m_date_r=item;
        case 'flag'
            if ischar(item)&&numel(item)~=0&&~ismember({item},{'*','id','js','cs','im','if','fw'})
                ME.message='Invalid flag. Must be a char with *, id, js, cs, im, fw, or if.';
                return
            end
        case 'UseLocalTime'
            [passed,opts.UseLocalTime]=test_if_scalar_logical(item);
            if ~passed
                ME.message='UseLocalTime should be either true or false';
                return
            end
        case 'UseURLwrite'
            [passed,item]=test_if_scalar_logical(item);
            if ~passed
                ME.message='UseURLwrite should be either true or false';
                return
            end
            %Force the use of urlwrite if websave is not available.
            opts.UseURLwrite=item || default.UseURLwrite;
        case 'err429'
            if ~isa(item,'struct')
                ME.message='The err429 parameter should be a struct.';
                return
            end
            %Loop through the fields in the input and overwrite defaults.
            opts.err429=default.err429;fn_=fieldnames(item);
            for n=1:numel(fn_)
                item_=item.(fn_{n});
                switch lower(fn_{n})
                    case 'countsastry'%'CountsAsTry'
                        [passed,item_]=test_if_scalar_logical(item_);
                        if ~passed
                            ME.message=['Invalid field CountsAsTry in the err429 parameter:',...
                                char(10),'should be a logical scalar.']; %#ok<CHARTEN>
                            return
                        end
                        opts.err429.CountsAsTry=item_;
                    case 'timetowait'%'TimeToWait'
                        if ~isnumeric(item_) || numel(item_)~=1
                            ME.message=['Invalid field TimeToWait in the err429 parameter:',...
                                char(10),'should be a numeric scalar.']; %#ok<CHARTEN>
                            return
                        end
                        %Under some circumstances this value is divided, so it has to be converted
                        %to a float type.
                        opts.err429.TimeToWait=double(item_);
                    case 'printatverbositylevel'%'PrintAtVerbosityLevel'
                        if ~isnumeric(item_) || numel(item_)~=1 || ...
                                double(item_)~=round(double(item_))
                            ME.message=['Invalid field PrintAtVerbosityLevel in the err429 ',...
                                'parameter:',char(10),'should be a scalar double integer.']; %#ok<CHARTEN>
                            return
                        end
                        opts.err429.PrintAtVerbosityLevel=item_;
                    otherwise
                        warning_(opts.print_to,'HJW:WBM:NameValue_not_found',...
                            ['Name,Value pair not recognized during parsing of err429 ',...
                            'parameter:\n    %s'],fn_{n});
                end
            end
        case 'waittime'
            try item=double(item);catch,end
            if ~isa(item,'double') || numel(item)~=1 || item<0
                ME.message='The waittime parameter should be a scalar convertible to double.';
                return
            end
            opts.waittime=item;
        case 'if_UTC_failed'
            valid={'error','ignore','warn_0','warn_1','warn_2','warn_3'};
            if isa(item,'string'),item=char(item);end
            if ~isa(item,'char') || ~ismember(lower(item),valid)
                ME.message=['The UTC_failed parameter is invalid.',char(10),...
                    'See the documentation for valid values.']; %#ok<CHARTEN>
                return
            end
            opts.UTC_failed=lower(item);
        case 'RunAsDate__fix'
            %Undocumented. You should probably not change this default.
            [passed,item]=test_if_scalar_logical(item);
            if ~passed
                ME.message='RunAsDate__fix is expected to be a scalar logical.';
                return
            end
            opts.RunAsDate__fix=item;
        case 'timeout'
            try
                item=double(item);
                if nargin(item)~=1 || item<0
                    error('trigger error')
                end
            catch
                ME.message='The timeout parameter should be a positive numeric scalar.';
                return
            end
            opts.timeout=item;
    end
end

%Set the behavior based on the UTC_failed setting and the verbosity.
if strcmp(opts.UTC_failed,'error')
    opts.UTC_fail_response='error';
elseif strcmp(opts.UTC_failed,'ignore')
    opts.UTC_fail_response='ignore';
else
    level=str2double(strrep(opts.UTC_failed,'warn_',''));
    if opts.verbose>=level
        opts.UTC_fail_response='warn';
    else
        opts.UTC_fail_response='ignore';
    end
end

%If the requested date doesn't match the current date, saves are not allowed, even if tries would
%suggest they are, so the code below checks this and sets tries(2) to 0 if needed.
%Because the server is in UTC (which might differ from the (local) time returned by the now
%function), the bounds need to be adjusted if the local time was used to determine the bounds.
currentUTC=WBM_getUTC_local;
if currentUTC==0,currentUTC=now;end % This will be dealt with elsewhere if needed.
if opts.UseLocalTime
    %Shift the bounds (which are currently in the local time) to UTC.
    timezone_offset=currentUTC-now;
    if opts.RunAsDate__fix && abs(timezone_offset)>(14.1/24)
        %This means Matlab does not know the correct time, as timezones range from -12 to +14.
        %Round down to true timezone offset. This is generally not what you would want.
        timezone_offset=rem(timezone_offset,1);
    end
    %Round the offset to 15 minutes (the maximum resolution of timezones).
    timezone_offset=round(timezone_offset*24*4)/(24*4);
    opts.date_bounds.datenum=opts.date_bounds.datenum-timezone_offset;
    bounds=opts.date_bounds.datenum;
else
    %The date_part is in already in UTC.
    bounds=opts.date_bounds.datenum;
end
%ML6.5 needs a specific format from a short list, the closest to the needed format (yyyymmddHHMMSS)
%is ISO 8601 (yyyymmddTHHMMSS).
item={datestr(bounds(1),'yyyymmddTHHMMSS'),datestr(bounds(2),'yyyymmddTHHMMSS')};
item{1}(9)='';item{2}(9)='';%Remove the T.
opts.date_bounds.double=[str2double(item{1}) str2double(item{2})];
if ~( bounds(1)<=currentUTC && currentUTC<=bounds(2) )
    %No saves allowed, because datepart doesn't match today.
    opts.tries(2)=0;
end

%If the m_date_r is set to error and the flag is set to something other than '*' the check_date
%function will return an error, regardless of the date stamp of the file.
%If that is the case, throw an error here.
if opts.m_date_r==2 && ~strcmp(opts.flag,'*')
    ME.message=['m_date_r set to ''error'' and the flag set to something other than ''*'' will',...
        char(10),'by definition cause an error, as the downloaded pages will not contain any',...
        char(10),'dates. See the help text for a work-around for images.']; %#ok<CHARTEN>
    ME.identifier='HJW:WBM:IncompatibleInputs';
    return
end
if ~ismember('m_date_r',replaced) && ~strcmp(opts.flag,'*')
    %If the default is not changed, but the flag is set to something else than '*', then the
    %m_date_r should be set to 0 (ignore).
    opts.m_date_r=0;
end
success=true;ME=[];
end
function opts=WBM_parse_inputs_defaults
%Create a struct with default values.
persistent opts_
if isempty(opts_)
    opts_=struct;
    opts_.date_part='2';%Load the date the closest to 2499-12-31 11:59:59.
    opts_.tries=[5 4 4];%These are the [loads saves timeouts] allowed.
    opts_.response={...
        'tx',404,'load'; %A 404 may also occur after a successful save.
        'txtx',[404 404],'save';
        'tx',403,'save';
        't2t2',[403 403],'exit'; %The page likely doesn't support the WBM.
        'tx',429,'pause_retry'...%Server is overloaded, wait a while and retry.
        };
    opts_.ignore=4080;
    opts_.verbose=3;
    opts_.UTC_failed='warn_3'; %Unused in default, but used when parsing options.
    opts_.UTC_fail_response='warn';
    opts_.m_date_r=1;%This selects 'warning' as default behavior.
    opts_.flag='*';
    opts_.UseLocalTime=false;
    %The websave function was introduced in R2014b (v8.4) and isn't built into Octave (webread was
    %introduced in 6.1.0, but websave not yet). As an undocumented feature, this can be forced to
    %true, which causes urlwrite to be used, even if websave is available. A check is in place to
    %prevent the reverse.
    try no_websave=isempty(which(func2str(@websave)));catch,no_websave=true;end
	opts_.UseURLwrite=no_websave;%(allows user-implementation in subfunction)
    opts_.err429=struct('CountsAsTry',false,'TimeToWait',15,'PrintAtVerbosityLevel',3);
    opts_.waittime=60;
    opts_.timeout=10;
    opts_.print_to_con=[];
    opts_.print_to_fid=[];
    opts_.print_to_obj=[];
    opts_.print_to_fcn=[];
    
    %(note for future edits: this code also exists in two parts elsewhere in this function)
    [bounds,midpoint]=WBM_parse_date_to_range(opts_.date_part);
    opts_.date_part=midpoint;opts_.date_bounds.datenum=bounds;
    item={datestr(bounds(1),'yyyymmddTHHMMSS'),datestr(bounds(2),'yyyymmddTHHMMSS')};
    item{1}(9)='';item{2}(9)='';
    opts_.date_bounds.double=[str2double(item{1}) str2double(item{2})];
    
    opts_.RunAsDate__fix=false;
end
opts=opts_;
end
function is_invalid=WBM_parse_inputs__validate_response_format(response)
%Check if the content of the response input is in the correct format.
%See doc('WBM') for a description of the correct format.
is_invalid=false;
if ~iscell(response) || isempty(response) || size(response,2)~=3
    is_invalid=true;return
end
for row=1:size(response,1)
    %Check col 1: t1, t2, tx or combination.
    item=response{row,1};
    item_count=numel(item(2:2:end));
    if ~ischar(item) || numel(item)==0 || ~all(ismember(item(2:2:end),'12x'))
        is_invalid=true;return
    end
    %Check col 2: html codes.
    item=response{row,2};
    if ~isa(item,'double') || numel(item)~=item_count
        %The validity of a html code is not checked.
        %A timeout caught in Matlab is encoded with 4080 (due to its similarity with HTML status
        %code 408).
        is_invalid=true;return
    end
    %Check col 3: load, save, exit or pause_retry.
    item=response{row,3};
    if ~ischar(item) || ~ismember({item},{'load','save','exit','pause_retry'})
        is_invalid=true;return
    end
end
end
function [timestamp,status,waittime]=WBM_retrieve_timestamp(url,date,opts,waittime)
% Attempt to retrieve the exact timestamp of a capture close to date_part.
% Usage is based on https://archive.org/help/wayback_api.php
%
% The timestamp is a char array and may be empty, in which case the status will be 404.
% The status is a double with either the HTTP status code returned by archive.org, or 404.

%Avoid problems caused by any & symbol in the url.
url=strrep(url,'&','%26');

avl_url=['http://archive.org/wayback/available?url=' url '&timestamp=' date];
if opts.UseURLwrite
    fn=tmpname('WBM_available_API_response','.txt');
end
for try_iterations=1:3
    try
        if opts.UseURLwrite
            fn=urlwrite(avl_url,fn); %#ok<URLWR>
            a=readfile(fn);delete(fn);
            a=JSON(a{1});
        else
            a=webread(avl_url);
        end
        break
    catch
        %Assume the connection is down, retry in intervals.
        connection_down_wait_factor=0;%Initialize.
        while ~isnetavl 
            if waittime.total<=waittime.cap
                %Total wait time exceeded, return the error condition.
                a=[];break
            end
            curr_time=datestr(now,'HH:MM:SS');
            if opts.verbose>=1
                msg=sprintf('Internet connection down, retrying in %d seconds (@%s)',...
                    2^connection_down_wait_factor,curr_time);
                warning_(opts.print_to,msg)
            end
            pause(2^connection_down_wait_factor)
            waittime.total=waittime.total+2^connection_down_wait_factor;
            %Increment, but cap to a reasonable interval.
            connection_down_wait_factor=min(1+connection_down_wait_factor,6);
        end
    end
end
try
    %This can fail for 4 reasons: the waittime has been exceeded (e.g. because there is no
    %connection in the first place), the connection is stable enough to get the ping (but not to
    %read the JSON response), an unparsed error occurred when reading the JSON, or the url simply
    %doesn't have a capture.
    timestamp=a.archived_snapshots.closest.timestamp;
    status=str2double(a.archived_snapshots.closest.status);
catch
    timestamp='';
    status=404;%A 404 would probably be the best guess.
end
end
function tf=WBM_UTC_test
%Returns a logical encoding whether the UTC determination can be expected to work.
persistent UTC_is_available
if isempty(UTC_is_available)
    UTC_is_available=WBM_getUTC_local>0;
    if ~UTC_is_available && ~isnetavl
        % The web method might have failed due to internet connection issues.
        UTC_is_available=[]; % Clear the persistent to test at the next call.
        tf=false;
    else
        tf=UTC_is_available;
    end
else
    tf=UTC_is_available;
end
end
function out=WinVer
%This returns the main Windows version number (5 for XP, 6 for Vista, etc).
%It will return an empty array for non-Windows or in case of errors.
persistent persistent_val
if ~ispc,out=[];return,end
if isempty(persistent_val)
    try
        [status,str] = system('ver'); %#ok<ASGLU>
        args={'[^0-9]*(\d*).*','$1','tokenize'};
        v=version;v(min(find(v=='.')):end)='';pre_v7=str2double(v)<7; %#ok<MXFND> 
        isOctave=exist('OCTAVE_VERSION','builtin')~=0;
        if ~pre_v7 && ~isOctave
            args(end)=[];%The 'tokenize' option became the default in R14 (v7).
        end
        persistent_val=str2double(regexprep(str,args{:}));
    catch
        try
            [status,str]=system('systeminfo'); %#ok<ASGLU>
            ind1= 1+strfind(str,':');ind1=ind1(3);
            ind2=-1+strfind(str,'.');ind2(ind2<ind1)=[];
            persistent_val=str2double(str(ind1:ind2(1)));
        catch
            persistent_val=[];
        end
    end
end
out=persistent_val;
end