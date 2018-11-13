addpath(genpath('helper_functions'));
clearvars; hold on; close all;

fs = 16000;
dftsize = 512;
window = ones(dftsize, 1);
noverlap = 0;             


% silence is send
sig1 = zeros(2*fs, 1);
[simin,nbsecs,fs] = initparams(sig1,fs);
sim('recplay');
rec1=simout.signals.values;
[sig_fft1, df_sig1, t_sig1, psd_sig1] = spectrogram(rec1, window, noverlap, dftsize, fs);


% signal, white noise to get all frequencies
sig2 = wgn(fs*2,1,0); %White noise
[simin,nbsecs,fs] = initparams(sig2,fs);
sim('recplay');
rec2=simout.signals.values;
[sig_fft2, df_sig2, t_sig2, psd_sig2] = spectrogram(rec2, window, noverlap, dftsize, fs);

p_noise = abs(mean(psd_sig1,2))';
p_signal = abs(mean(psd_sig2-psd_sig1,2))';
%ongecorelleerde PSD, dus mag je gewoon optellen of aftrekken

figure('Name', 'PSD');
subplot(2,1, 1);
plot(df_sig1, 10*log10(p_noise));  title( 'noise PSD' ); xlabel( 'frequency(Hz)' ); ylabel( 'magnitude (dB)' );
subplot(2,1,2);
plot(df_sig2, 10*log10(p_signal)); title( 'signal PSD' ); xlabel( 'frequency(Hz)' ); ylabel( 'magnitude (dB)' );

N = dftsize/2;
ratio = p_signal./p_noise;
c_chan = sum(fs/(2*N)*log2(1+ratio));
%% 1 
% dell optiplex 3060 ... dell U989P (internal speaker)
% All frequencies audible by humans

%% 2/3
% done, see file 'plot_shannon_vs_distance' for more results, but results of
% magnitude of 10^5

%% 4
% Capacity of channel is a measure of how many information can be send over
% the channel = bits/s/Hz

%% 5 
% larger fs, means larger bandwidth, means larger Capacity of the channel
% and thus you are able to send more data if sampling frequency is higher.
% Until you reach maximum of hardware abilities.
% 2.263776353017300e+05 for 2-ish cm
