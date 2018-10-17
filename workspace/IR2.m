addpath(genpath('helper_functions'));
clearvars; close all;

%% 2.2.1 Relation between u[k], y[k] and h
% y[k] = h(k)*u(k) (* is convolutie)
% Y[z] = H[z].U[z]

%% 2.2.2 Estimate IR h
%parameters
fs=16000;
h_order = 80; %(bepaald bij IR1.m)  
sig_time = 2;

% Generate input signal (witte ruis), play and record
% om het signaal uit de output te kunnen halen
sig =  wgn(fs*sig_time, 1, 0);
[simin, nbsecs, fs] = initparams(sig, fs);
sim('recplay'); 
rec = simout.signals.values;

start = find(rec > 0.3*max(rec), 1, 'first')-10; %geeft de eerst index die hieraan voldoet
%start is de index waar het signaal (ruis) groter is dan de ruis die er
%altijd is --> zeker dat het signaal daar start

y = rec(start :start+sig_time*fs+h_order); %geselecteerde output
u = sig; %input

r = sig(1 : h_order,1); 
U = toeplitz(u.', (flipud(r)).');
