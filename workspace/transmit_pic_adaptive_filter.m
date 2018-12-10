addpath(genpath('helper_functions'), genpath('mod_func'), genpath('data'));
clearvars; hold on; close all;

%% Params
alpha =0.2;
mu = 1;
K = 4;
N = 512;
fs = 16000;
SNR = 60;
L = 160; %channel length
cp_size = L+16;
Lt = 5;
Ld = 4;
BW = 75;
threshold = BW/100;

alpha_test = [0.01; 0.1; 1; 10];
mu_test = [0.01; 0.1; 1; 5];

%% Image==Data
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');


%% trainingsblock
n_trainblock    = (N/2-1)*K;   % contains N/2-1 random QAM symbols
trainblock_bit  = randi([0 1], 1, n_trainblock);
trainblock_qam  = qam_mod(trainblock_bit, K);

%% channel est

seq = randi([0 1], (N/2-1)*K, 1); % random bits
trainblock = qam_mod(seq, K); 
trainblocks = repmat(trainblock,20,1);
Tx = ofdm_mod_tb(trainblocks, N, cp_size);

%pulse  =  wgn(fs*0.1, 1, 0); % pulse = 0.1 sec of noise
[pulse, ~] = sinusoid(1000, 1, 1, fs);
[simin,nbsecs,fs,pulse ] = initparams(Tx,fs,L, pulse);

sim('recplay'); 
rec = simout.signals.values;

[Rx_2] = alignIO(rec,pulse,L);
Rx_2 = Rx_2(1:length(Tx));

[seq_qam , H] = ofdm_demod_tb(Rx_2, N, cp_size, trainblock);

receive_deqam = qam_demod(seq_qam, K); 
trainblocks_deqam = qam_demod(trainblocks, K); 

freq_bins = ofdm_freq_bins(H, N, threshold);
%freq_bins = ones((N/2-1),1);


%% Modulate (QAM amd OFDM)
qamStream_in = qam_mod(bitStream, K);

Tx = ofdm_mod_adaptive_filter(qamStream_in, N, cp_size, freq_bins, trainblock_qam, Lt);

%% Doosturen 

%pulse  =  wgn(fs*0.1, 1, 0); % pulse = 0.1 sec of noise
[pulse, ~] = sinusoid(1000, 1, 1, fs);
[simin,nbsecs,fs,pulse ] = initparams(Tx,fs,L, pulse);

sim('recplay'); 
rec = simout.signals.values;

%% Align
[Rx] = alignIO(rec,pulse,L);
Rx = Rx(1:length(Tx)+L+10);

%% Demod
[rxQamStream, H ] = ofdm_demod_adaptive_filter(Rx, N, cp_size,freq_bins, trainblock_qam, Lt, mu,K,alpha);


%% QAM demodulation
rxBitStream = qam_demod(rxQamStream,K);

%% Compute BER
berTransmission = ber(bitStream,rxBitStream);

% plot van ber per transmissie
M = sum(freq_bins);
P = ceil(size(qamStream_in, 1) / M);
bertrain = zeros(floor(P/Ld),1);
for b=0:floor(P/Ld)-1
    bertrain(b+1)=ber(bitStream(b*Ld*M*K+1:(b+1)*Ld*M*K),rxBitStream(b*Ld*M*K+1:(b+1)*Ld*M*K));
end


figure('Name', 'BER vs every transmission');
title( 'BER vs transmission' ); xlabel( 'Number of transmission(/)' ); ylabel( 'Bit Error Rate(/)' );
hold on
plot(bertrain,'o')
plot(bertrain,'r-')
plot(berTransmission*ones(length(bertrain),1), 'm--', 'LineWidth', 1)
hold off

%% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);

%% Plot images
figure
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;


%% vraagjes 
% 

