addpath(genpath('helper_functions'));
clearvars; hold on; close all;

fs = 16000;
pulse_freq = 1500;
[sinewave, t] = sinusoid(pulse_freq, 1, 2, fs);

[simin,nbsecs,fs] = initparams(sinewave,fs);

sim('recplay');

out=simout.signals.values;

[simin,nbsecs,fs] = initparams(out,fs);

sim('recplay');


