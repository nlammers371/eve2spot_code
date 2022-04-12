% Execute this script to call all data processing functions for the eve 2 
% spot project
clear
close all

% add functions to path
addpath(genpath('./functions'))

% combine embryos into single master dataset
masterSet = compile_traces;