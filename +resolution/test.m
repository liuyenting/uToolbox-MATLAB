clear all; close all; %#ok<CLALL>

A = 1:20;
A = repmat(A', [1, 3]);

S = shuffle(A, 2, 5);
