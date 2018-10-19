addpath(genpath('helper_functions'));
clearvars; hold on; close all;

fs = 16000;
dftsize = 128;
window = ones(dftsize, 1);  % Rectangular window
noverlap = 4;               % N samples to overlap between segments


% silence is send
sig1 = zeros(2*fs, 1);
[simin,nbsecs,fs] = initparams(sig1,fs);
sim('recplay');
rec1=simout.signals.values;
[sig_fft1, df_sig1, t_sig1, psd_sig1] = spectrogram(rec1, window, noverlap, dftsize, fs);


% signal, white noise to get all frequencies
sig2 = wgn(fs*2,1,0); %White noise;
[simin,nbsecs,fs] = initparams(sig2,fs);
sim('recplay');
rec2=simout.signals.values;
[sig_fft2, df_sig2, t_sig2, psd_sig2] = spectrogram(rec2, window, noverlap, dftsize, fs);

p_noise = mean(psd_sig1,2)';
p_signal = mean(psd_sig2-psd_sig1,2)';
%ongecorelleerde PSD, dus mag je gewoon optellen of aftrekken

figure('Name', 'PSD');
subplot(2,1, 1);
plot(df_sig1, 10*log10(p_noise));  title( 'noise PSD' ); xlabel( 'frequency(Hz)' ); ylabel( 'magnitude (dB)' );
subplot(2,1,2);
plot(df_sig2, 10*log10(p_signal)); title( 'signal PSD' ); xlabel( 'frequency(Hz)' ); ylabel( 'magnitude (dB)' );

c_chan = 0;
N = dftsize/2;
ratio = zeros(N,1);
for m = 1:N
    ratio(m) = p_signal(m)/p_noise(m);
end

for k = 1:N
    c_chan = c_chan + fs/(2*N)*log2(1+ratio(k));
end


%% 1 
% dell optiplex 3060 ... dell U989P (internal speaker)
% All frequencies audible by humans

%% 2/3
% done, 9.736537483390578e+04 for 2-ish cm

%% 4
% Capacity of channel is a measure of how many information can be send over
% the channel = bits/s/Hz

%% 5 
% larger fs, means larger bandwidth, means larger Capacity of the channel
% and thus you are able to sned more data if sampling frequency is higher.
% Until you reach maximum of hardware abilities.
% 2.263776353017300e+05 for 2-ish cm
