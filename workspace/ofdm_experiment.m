addpath(genpath('helper_functions'), genpath('mod_func'), genpath('data'));
clearvars; hold on; close all;

%% PARAMETERS
K = 4 ; % Number of bits per constellation
M = 2^K;
n = 256*16*K; % Number of data points



