clc
try
    if exist('foo.db','file')
        delete('foo.db');
    end
catch
    error( 'Unable to delete database' );
end

% creating a table
sqlite3('foo.db',...
    'CREATE TABLE test (some_text TEXT, some_int INT, some_real REAL);');

% inserting entries of a struct
x=struct;
x(1).text = 'world';
x(1).int = 1337;
x(1).dbl = 2.71828;
x(2).text = 'foobar';
x(2).int = -4131;
x(2).dbl = 9.0;
x(3).text = 'foo-bar';
x(3).int = 404;
x(3).dbl = 2*pi;
sqlite3('foo.db',...
    'INSERT INTO test (some_text, some_int, some_real) VALUES (?, ?, ?);', x)

% inserting a string query
sqlite3('foo.db',...
    'INSERT INTO test (some_text, some_int, some_real) VALUES ("hello", 42, 3.14159);')

% querying
y = sqlite3('foo.db',...
    'SELECT some_text, some_int, some_real as real FROM test;');
z=sqlite3('foo.db','SELECT * FROM test;');

y(1).some_text
y(2).some_int
y(4).real
disp(z)