addpath(genpath('helper_functions'));
clearvars; close all;

%% 2.2.1 Relation between u[k], y[k] and h
% y[k] = h(k)*u(k) (* is convolutie)
% Y[z] = H[z].U[z]

%% 2.2.2 Estimate IR h
%%%%%parameters
fs=16000;
h_order = 200;                                      %(bepaald bij IR1.m)  
sig_time = 2;

%%%%%Generate input signal (witte ruis), play and record
                                                    %(om het signaal uit de output te kunnen halen)
sig =  wgn(fs*sig_time, 1, 0);
[simin, nbsecs, fs] = initparams(sig, fs);
sim('recplay'); 
rec = simout.signals.values;

start = find(rec > 0.5*max(rec), 1, 'first')-20;    %geeft de eerst index die hieraan voldoet
                                                    %start is de index waar het signaal (ruis) groter is dan de ruis die er
                                                    %altijd is --> zeker dat het signaal daar start

y = rec(start+h_order :start+sig_time*fs+h_order);  %geselecteerde output
u = simin(2*fs+1+h_order:(2+sig_time)*fs+h_order+1);%input: neem de input vanaf waar je het
                                                    %signaal opstart (2*fs+1) en tot het einde van het
                                                    %signaal((2+sig_time)*fs)+ h_order want het sijpelt zoveel na
r = flipud(simin(2*fs+1 : 2*fs+h_order+1,1)).';     %moet geflipt en neen gewoon h_order owv definitie toeplitz elememten
U = toeplitz(u.', r);

%%%%%least squares 
h = mldivide(U, y);                                 %NOTE: least-square solutions in same plane -> x = A^(-1)*b 

%%%%%plot
t = (0: 1 : (length(h)-1)) /fs;                     % is een rijvector
fft_h = fft(h)';                                    % transpose want we willen een rijvector om te kunnen plotten
length_h = round((length(fft_h)-1)/2 );             % Take half the length
fft_h_trunc = abs(fft_h(1 : length_h));             % Use only first half of fft samples (=positive frequencies).
df_h = (0 : fs /(2*length_h) : (fs-1)/2);           % Frequency vector for x-axis in plots

figure;
subplot(2, 1, 1);
    plot(t, h);
    title('impulse response');
    xlabel('time (sec)');
    ylabel('magnitude');
subplot(2, 1, 2);
    plot(df_h, 20*log10(fft_h_trunc));
    title('frequency response');
    xlabel('frequency(Hz)');
    ylabel('magnitude (dB)');

%%%%%Save to file 'IRest.mat'
%     matfile = fullfile('./data', 'IRest.mat');
%     save(matfile, 'h');
    
%% 2.2.3
% Ze lijken op elkaar omdat je in hetzelfde kanaal aan het meten bent.

%% 2.2.4
% Heel gelijkaardig
% kleine verschillen owv verschillende methode van opmeten
% (kanaal vs. freq respons)

%% 2.2.5
% langere IR WANT meer interferentie, plaats van opnemen geeft
% verschillende delay 
% --> KORTOM: moeilijker op te meten

% verify experimentally impossible, don't have stereos

%% 2.2.6
% sound waves are diverted SO longer path from speaker to microphone
% + specific freq are absorbed more than other by your hand
% --> results in different TF

% experiment: toooooo much distortion, can't handle this

%% 2.3.1

% does not change a lot, the channel stays the same over all experiments 
% (not lik enviromental noise that changes all the time
% (see 2.3.4 (bandstop)))
