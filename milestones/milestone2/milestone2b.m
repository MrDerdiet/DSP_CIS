addpath(genpath('helper_functions'), genpath('mod_func'), genpath('data'));
clearvars; hold on; close all;

%% params
K = 6;
N = 512;

SNR = 5;
L = 160; %channel length
cp_size = L+16;

%% Convert BMP image to bitstream
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');

%% QAM modulation
qamStream = qam_mod(bitStream,K);

%% Channel -> with fftfilt

load('h_ir2','h');
h = [h; zeros(N -size(h, 1), 1)];

Hn_vector = fft(h);

[freq_bins] = ofdm_freq_bins(Hn_vector, N, 0.1);



%% OFDM modulation
ofdmStream = ofdm_mod_onoff(qamStream, N, cp_size,freq_bins);


%% channel
ofdmStream = fftfilt(h, ofdmStream);                % alsof het signaal over het kanaal wordt gestuurd (signaal*TF)
rxOfdmStream = awgn(ofdmStream, SNR, 'measured');   % witte ruis erop

%% OFDM demodulation
rxQamStream = ofdm_demod_onoff(rxOfdmStream, N, cp_size, Hn_vector,freq_bins);

%% QAM demodulation
rxBitStream = qam_demod(rxQamStream,K);

%% Compute BER
rxBitStream = rxBitStream(1:length(bitStream));
berTransmission = ber(bitStream,rxBitStream);


%% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);

%% Plot images
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;


