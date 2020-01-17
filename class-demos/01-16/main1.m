%% In Class: Types in Matlab

clear all
close all
clc

%% Char
% character arrays
a = 'this is a char...raw'

% these are actually arrays of numbers with a mask that maps the alphabet
% to the numbers.  Take a look at what happens when we do this:
a + a

%% Strings
% these are newer to Matlab (from 2016b) and mostly look and walk like
% character arrays.  See the website for FAQ on differences between strings
% and chars
% https://www.mathworks.com/help/matlab/matlab_prog/frequently-asked-questions-about-string-arrays.html

s1 = "This is a string";

% strings act differently from character arrays, so see what happens here
s1 + s1

% That's called a concatenation.


%% Structures
% very useful named fields
b.x= 1;
b.y= 7;
b.z = rand(2,4);
b.w = 'strings';

b.b.x = b.x;

% We can also use the struct function as the constructor.  It takes in
% pairs of info "name" and "value."  The name can be either a string or a
% char array, and value is anything we want to store.
x = struct("field1", 7, 'field2', rand(3,2) )

%% Cells
% Like a list that can contain anything

c = {1, rand(4), 'string', b}

d = {'long string skjdhfksjdfksjdf',
     'short string'}

% while cells are the classic way to store text of different length, the
% modern use of strings should actually be done with arrays.
d_str = [ "long string skjdhfksjdfksjdf",
          "short string"]


