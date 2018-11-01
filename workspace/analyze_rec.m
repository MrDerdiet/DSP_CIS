addpath(genpath('helper_functions'));
clearvars; hold on; close all;


fs = 16000;
dftsize = 256;
pulse_freq = 1500;


%{
[sig, t] = sinusoid(pulse_freq, 1, 2, fs);



sig = wgn(fs*3,1,0); %White noise

final_sig = sinusoid(0, 1, 2, fs);
freq_array = [100,200,500,1000,1500,2000,4000,6000];
for i = 1:8
    [sig, t] = sinusoid(freq_array(i), 1, 2, fs);
    final_sig = final_sig + sig;
end
sig = final_sig;




%}

%sig = wgn(fs*3,1,0); %White noise


% WINDOWS
% windows1 = sinusoid(0, 1, 0.5, fs);
% windows2 = sinusoid(0, 1, 0.5, fs);
% windows3 = sinusoid(0, 1, 0.6, fs);
% windows4 = sinusoid(0, 1, 0.25, fs);
% windows5 = sinusoid(0, 1, 0.5, fs);
% freq_array1 = [830,1245,2489];
% for i = 1:3
%     [sig, t] = sinusoid(freq_array1(i), 1, 0.5, fs);
%     windows1 = windows1 + sig;
% end
% freq_array2 = [830,1245,1865,2489];
% for i = 1:4
%     [sig, t] = sinusoid(freq_array2(i), 1, 0.5, fs);
%     windows2 = windows2 + sig;
% end
% freq_array3 = [830,1245,1661];
% for i = 1:3
%     [sig, t] = sinusoid(freq_array3(i), 1, 0.6, fs);
%     windows3 = windows3 + sig;
% end
% freq_array4 = [2489];
% for i = 1:1
%     [sig, t] = sinusoid(freq_array4(i), 1, 0.25, fs);
%     windows4 = windows4 + sig;
% end
% freq_array5 = [933,1245,1865];
% for i = 1:3
%     [sig, t] = sinusoid(freq_array5(i), 1, 0.5, fs);
%     windows5 = windows5 + sig;
% end
% 
% sig = [windows1; windows2; windows3; windows4; windows5];
% 


fs = 16000;
t = 0:1/fs:0.3-1/fs;

%   niets do=1 re=2   mi=3   fa=4   sol=5 la=6  si=7  do=8   re=9   mi=10  fa=11 sol=12 la=13  si=14
l = [0 130.81 146.83 164.81 174.61 196.00 220 246.94 261.63 293.66 329.63 349.23 392.00 440 493.88];
h = [0 523.25 587.33 659.25 698.46 783.99 880 987.77 1046.50 1174.66 1318.50 1396.92 1567.98 1760 1975.54];
note = @(f,g) [1 1]*sin(2*pi*[l(f) h(g)]'*t);

low  = [3 5 6 6 6 6 6 7 8 8 8 8 8 9 7 7 7 7 6 5 5 6 0 0 3 5 6 6 6 6 6 7 8 8 8 8 8 9 7 7 7 7 6 5 6 6 6 6]+1;
high = [3 5 6 6 6 6 6 7 8 8 8 8 8 9 7 7 7 7 6 5 5 6 0 0 3 5 6 6 6 6 6 7 8 8 8 8 8 9 7 7 7 7 6 5 6 6 6 6]+1;

song = [];
for kj = 1:length(low)
    song = [song note(low(kj),high(kj)) zeros(1,0.01*fs)];
end
song = song/(max(abs(song))+0.1);
sig = song';


%% (1a) Play the audio signal ('sig') and record (to 'rec') using Simulink 
% Simulink needs input 'simin' to be defined in workspace

[simin,nbsecs,fs] = initparams(sig,fs);

sim('recplay');

rec=simout.signals.values;

%% (1b)  Spectrogram of 'sig' and mic recording 
% plot one figure below each other 

window = ones(dftsize, 1);          % Rectangular window
noverlap = 16;                      % N samples to overlap between segments
                                    % overlap geeft nauwkeurigheid

figure('Name', 'spectograms');      % fig1
subplot( 2, 1, 1 ); 
spectrogram(sig, window, noverlap, dftsize, fs, 'yaxis'); 
[sig_fft, df_sig, t_sig, psd_sig] = spectrogram(sig, window, noverlap, dftsize, fs);
title( 'original signal spectogram' );

subplot( 2, 1, 2 );
spectrogram(rec, window, noverlap, dftsize, fs, 'yaxis');
[rec_fft, df_rec, t_rec, psd_rec] = spectrogram(rec, window, noverlap, dftsize, fs);
title( 'recorded signal spectogram' );

%% (1c) separate figure: PSD of the transmitted and recorded signal
% Note:  PSD =~ average power spectrum over the full signal length 
% spectrogram already calculated power of instantaneous moment
psd_sig_avg = mean(psd_sig, 2)';
psd_rec_avg = mean(psd_rec, 2)';

figure('Name', 'PSD'); %fig2
subplot(2,1, 1);
plot(df_sig, 10*log10(psd_sig_avg));  title( 'original signal PSD' ); xlabel( 'frequency(Hz)' ); ylabel( 'magnitude (dB)' );
subplot(2,1,2);
plot(df_rec, 10*log10(psd_rec_avg)); title( 'recorded signal PSD' ); xlabel( 'frequency(Hz)' ); ylabel( 'magnitude (dB)' );

%% (2) Look at PSD and spectrogram sig. Expected? Other frequencies than 400Hz present? Why (not)?
% Peak @400Hz, some frequency smearing (fft), minder nauwkeurig
% (+finite precision) + harmonische frequenties
% Synthetic signal => no other frequencies in signal itself


%% (3) Influence of DFT size? Try to answer first and then verify experimentally.
% Frequency resolution ~ fs/dftsize
% Frequency resolution ~1/Time resolution
% afweging tussen nauwkeurigheid in tijdsdomein of frequentie domein:
% in tijdsdomein op elk ogenblik precies => frequentiedomein uitgesmeerd


%% (4) Compare spectogram/PSD of original signal and received signal
% spectogram : Harmonischen van basisch frequentie zijn aanwezig in
% ontvangen signaal, eveneens als ruis over heel het frequentiedomein
% PSD : origineel signaal heeft 1 duidelijk signaal, terwijl ontvangen
% signaal ook pieken heeft op resonantie frequenties en ruis


%% (5) Compare spectogram/PSD of original signal and received signal
% bij verzonden, origineel signaal is wel een DC-componentn aanwezig, maar
% bij ontvangen signaal is deze niet meer aanwezig.


%% (6) Influence of scaling before transmitting
% Clipping zorgt voor extra pieken in frequentiedomein van het het ontvangen
% signaal. (Dit komt door distortie) Om dit te demonstreren moet herschaling
% uitgezet worden in "initparams". (vb: sinusoid met amplitude 2)


%% (7) Insert automatic scaling
% Hebben we in "initparams" gedaan.  :-D


%% (8) Multiple sines added up to one signal
% Meerder pieken op de gegeven frequenties bij het origineel signaal. Bij
% het ontvangen signaal zijn deze pieken boven 1000Hz nog goed
% waarneembaar. harmonischen van lagere frequenties versterken ook
% veelvouden van de basisch frequentie. 


%% (9+10+11) White noise experiments
% Uitgezonden signaal bestaat nu uit Gaussische witte ruis. Toch zijn bij
% ontvangst paar frequenties beter ontvangen dan anderen, zelfs bij het
% bewegen van de microfoon. Het bewegen geeft enkel in tijdsdomein soms een
% sterker/zwakker signaal (als je dichtbij of ver van de luidspreker zit), maar dezelfde frequenties blijven aanwzig bij ontvangst
% (vooral 6000Hz). Dit gedrag blijft terug komen voor verschillende experimenten.
% Waarschijnlijk eigenschappen van de microfoon, die zo van 6000 Hz houdt!




