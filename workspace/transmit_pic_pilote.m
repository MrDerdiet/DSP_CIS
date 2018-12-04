addpath(genpath('helper_functions'), genpath('mod_func'), genpath('data'));
clearvars; hold on; close all;

%% Params

K = 4;
N = 1024;
fs = 16000;
SNR = 60;
L = 160; %channel length
cp_size = L+16;
Lt = 1;
Ld = 4;
BW = 80;
threshold = BW/100;

%% Image==Data
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');


%% Pilote 
n_pilote    = (N/4)*K;   % contains N/2-1 random QAM symbols
pilote_bit  = randi([0 1], 1, n_pilote);
pilote_qam  = qam_mod(pilote_bit, K);

%% freq bins bepalen
seq = randi([0 1], (N/2-1)*K, 1); % random bits
trainblock = qam_mod(seq, K); 
trainblocks = repmat(trainblock,10,1);
Tx = ofdm_mod_tb(trainblocks, N, cp_size);

[pulse, ~] = sinusoid(1000, 1, 1, fs);
[simin,nbsecs,fs,pulse ] = initparams(Tx,fs,L, pulse);
sim('recplay'); 
rec = simout.signals.values;
[Rx_2] = alignIO(rec,pulse,L);
Rx_2 = Rx_2(1:length(Tx));

[seq_qam , H] = ofdm_demod_tb(Rx_2, N, cp_size, trainblock);
%aanpassen naar enkel de even freq
pilotes = 1:2:N/2-1;
H(pilotes) = 0;
freq_bins = ofdm_freq_bins(H, N, threshold/2);


%% Modulate (QAM amd OFDM)
qamStream_in = qam_mod(bitStream, K);
Tx = ofdm_mod_pilote(qamStream_in, N, cp_size, freq_bins, pilote_qam);


%% Doosturen 
tic;
%pulse  =  wgn(fs*0.1, 1, 0); % pulse = 0.1 sec of noise
[pulse, ~] = sinusoid(1000, 1, 1, fs);
[simin,nbsecs,fs,pulse ] = initparams(Tx,fs,L, pulse);

sim('recplay'); 
rec = simout.signals.values;
toc;
%% Align
[Rx] = alignIO(rec,pulse,L);
Rx = Rx(1:length(Tx)+L+10);

%% Demod
[rxQamStream, H ] = ofdm_demod_pilote(Rx, N, cp_size,freq_bins,pilote_qam);

rxQamStream = rxQamStream(1:length(qamStream_in));

%% Channel freq response
H_fft = [zeros(1,size(H,2)) ; H ; zeros(1,size(H,2));flipud(conj(H))];
h = ifft(H_fft);
H_db = 20*log10(abs(H));

%% QAM demodulation
rxBitStream = qam_demod(rxQamStream,K);

%% Compute BER
berTransmission = ber(bitStream,rxBitStream);
fprintf(['With on-off bit loading: Ber =', num2str(berTransmission),'\n']);
% plot van ber per transmissie
M = sum(freq_bins);
P = ceil(size(qamStream_in, 1) / M);
bertrain = zeros(P,1);

for b=0:floor(P)-2
    bertrain(b+1)=ber(bitStream(b*M*K+1:(b+1)*M*K),rxBitStream(b*M*K+1:(b+1)*M*K));
end
b=b+1;
bertrain(b+1)=ber(bitStream(b*M*K+1:end),rxBitStream(b*M*K+1:end));

figure('Name', 'BER vs every transmission');
title( 'BER vs transmission' ); xlabel( 'Number of transmission(/)' ); ylabel( 'Bit Error Rate(/)' );
hold on
plot(bertrain,'o')
plot(bertrain,'r-')
plot(berTransmission*ones(length(bertrain),1), 'm--', 'LineWidth', 1)
hold off

%% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);



%% Visualisation
time = N/fs;

max_h = max(h(:));
max_H = max(H_db(:));
min_H = min(H_db(:));
figure('Name','visualisation of the demodulation')
subplot(2,2,2);colormap(colorMap); image(imageData); axis image; title('The transmitted image'); 
for i = 1:size(h,2)-1 
    subplot(2,2,1)
    plot(h(:,i))
    title('Channel in time domain')
    axis([0 160 -max_h max_h])
    
    subplot(2,2,3)
    plot((1:N/2-1)*fs/2/(N/2-1),H_db(:,i));
    title('Channel in frequency domain')
    axis([-inf, inf, min_H max_H]);
    ylabel('Magnitude')
    xlabel('Frequency')
    
    if i*sum(freq_bins)*K<length(rxBitStream)
        imageRx = bitstreamtoimage(rxBitStream(1:i*sum(freq_bins)*K), imageSize, bitsPerPixel);
    else
        imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);
    end
    subplot(2,2,4); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;
    pause(time)
end    

