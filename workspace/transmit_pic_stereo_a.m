clearvars; hold on; close all;
addpath(genpath('help_func'), genpath('data'), genpath('mod_func'));

%% Convert BMP image to bitstream
[bitStream_in, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');

tic;

%% Parameters
fs = 16000;
N = 512; 
K = 6; 
BW = 75;
threshold = BW/100;


% h1 = [1: -1/159: 0].';
% h1 = h1/(sum(h1));
% h2 = h1;
h1 = rand(160,1);
h1 = h1/sum(h1);
h2 = h1;
% h2 = [1; 0; zeros(N-2,1)];
H = fft(h1);
freq_bins = ofdm_freq_bins(H, N, threshold);
cp_size = 176;


%% Modulate
qamStream_in = qam_mod(bitStream_in, K);



%% a) Generate two random impulse responses of used de?ned length 
%(modeling the acoustic channels from loudspeaker one and two to the 
% microphone). Are these a good model for the acoustic channels? Why (not)?

% -> miss er laten uitzien als goede h (das mooier)
% h1 = rand(160,1)/100;
% h1 = [1: -1/159: 0].';
% h1 = h1/(sum(h1));
% h2 = h1;
% h2 = rand(160,1);
h1 = [1; zeros(159,1)];
h2 = [zeros(160,1)];


%% b) -> ofdm_mod_stereo.m
% [ofdmStream1 ,ofdmStream2,seq_ofdm] = ofdm_mod_stereo(seq_qam, N, cp_size, freq_bins,a,b)
H1 = fft(h1, N);
H2 = fft(h2, N);

numerator= sqrt(H1.*conj(H1)+H2.*conj(H2));

a = conj(H1)./numerator;
b = conj(H2)./numerator;
% seq_qam = qamStream_in_pad;
[ofdmStream1, ofdmStream2] = ofdm_mod_stereo(qamStream_in, N, cp_size, freq_bins, a, b);

%% c)  Generate the received microphone signal by convolving the two OFDM 
% signals generated in 1b with the two impulse responses and adding the 
% results. Do not (yet) add any noise to the microphone signal!

ofdmStream1_Rx_no_noise = fftfilt(h1, ofdmStream1);
ofdmStream2_Rx_no_noise = fftfilt(h2, ofdmStream2);

ofdmStream1_Rx = awgn(ofdmStream1_Rx_no_noise, 20, 'measured');
ofdmStream2_Rx = awgn(ofdmStream2_Rx_no_noise, 20, 'measured');

Rx = ofdmStream1_Rx + ofdmStream2_Rx;

%% d)  ofdm_demod_stereo.m
seq_qam = ofdm_demod_stereo(Rx, N, cp_size, freq_bins, H1, H2);
qamStream_out = seq_qam(1: length(qamStream_in));

bitStream_out = qam_demod(qamStream_out, K);


%% e)  check BER You can only continue when the BER is exactly zero!

[ber_IO, error_vector] = ber(bitStream_in, bitStream_out);

%% f) (f) Now add noise to the received signal and check the BER in the 
% case of mono transmission with loudspeaker 1 (a = 1,b = 0), mono 
% transmission with loudspeaker 2 (a = 0,b = 1) and stereo transmission 
% with both loudspeaker 1 and 2 (optimal a and b). Do you see an advantage 
% of using two loudspeakers


