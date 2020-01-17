%% Types

clear all
close all
clc

%% Char
% character arrays
a = 'this is a char...raw'


%% Structures
% very useful named fields
b.x= 1;
b.y= 7;
b.z = rand(2,4);
b.w = 'strings';

b.b.x = b.x;


%% Cells
% Like a list that can contain anything

c = {1, rand(4), 'string'}

d = {'long string skjdhfksjdfksjdf',
     'short string'}
 
%%
