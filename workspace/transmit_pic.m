addpath(genpath('helper_functions'), genpath('mod_func'), genpath('data'));
clearvars; hold on; close all;

%% Params

K = 6;
N = 512;
fs=16000;
SNR = 60;
L = 160; %channel length
cp_size = L+16;
Lt = 1;
Ld = 6;

%% Image==Data
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');


%% trainingsblock
n_trainblock    = (N/2-1)*K;   % contains N/2-1 random QAM symbols
trainblock_bit  = randi([0 1], 1, n_trainblock);
trainblock_qam  = qam_mod(trainblock_bit, K);


%% Modulate (QAM amd OFDM)
qamStream_in = qam_mod(bitStream, K);

freq_bins = ones(N/2-1,1); 

Tx = ofdm_mod(qamStream_in, N, cp_size, freq_bins, trainblock_qam, Lt, Ld);

%% Doosturen 

%pulse  =  wgn(fs*0.1, 1, 0); % pulse = 0.1 sec of noise
[pulse, ~] = sinusoid(1000, 1, 1, fs);
[simin,nbsecs,fs,pulse ] = initparams(Tx,fs,L, pulse);

sim('recplay'); 
rec = simout.signals.values;

%% Align
[Rx] = alignIO(rec,pulse,L);
%Rx = Rx(1:length(Tx)+L+10);

%% Demod
[rxQamStream, H ] = ofdm_demod(Rx, N, cp_size,freq_bins, trainblock_qam, Lt, Ld,);


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


%% vraagjes 
% 

