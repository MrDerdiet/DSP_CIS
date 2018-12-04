addpath(genpath('helper_functions'), genpath('mod_func'), genpath('data'));
clearvars; hold on; close all;


%% Params
K = 4;
N = 512;
fs=16000;
SNR = 60;
L = 160; %channel length
cp_size = L+16;
Lt = 2;
Ld = 16;
BW = 50;
threshold = BW/100;

%% Image==Data
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');

% %%% Zonder onoff
% %% trainingsblock
n_trainblock    = (N/2-1)*K;   % contains N/2-1 random QAM symbols
trainblock_bit  = randi([0 1], 1, n_trainblock);
trainblock_qam  = qam_mod(trainblock_bit, K);

%% H-channel voor freq bins 
freq_bins = ones(N/2-1,1); 

%% Modulate (QAM amd OFDM)
qamStream_in = qam_mod(bitStream, K);
Tx = ofdm_mod(qamStream_in, N, cp_size, freq_bins, trainblock_qam, Lt, Ld);

%% Play en record
tic;
[pulse, ~] = sinusoid(1000, 1, 1, fs);
[simin,nbsecs,fs,pulse ] = initparams(Tx,fs,L, pulse);

sim('recplay'); 
rec = simout.signals.values;
toc;
%% Align

[Rx] = alignIO(rec,pulse,L);
Rx = Rx(1:length(Tx)+L+10);

%% Demodulate (QAM amd OFDM)
[rxQamStream, H ] = ofdm_demod(Rx, N, cp_size,freq_bins, trainblock_qam, Lt, Ld);
rxQamStream = rxQamStream(1:length(qamStream_in));
rxBitStream = qam_demod(rxQamStream,K);

%% Channel freq response
H_fft = [zeros(1,size(H,2)) ; H ; zeros(1,size(H,2));flipud(conj(H))];
h = ifft(H_fft);
H_db = 20*log10(abs(H))./freq_bins; % niet gebruikte freq bins op 0
 
%% BER
berTransmission = ber(bitStream,rxBitStream);
fprintf(['Without on-off bit loading: Ber =', num2str(berTransmission),'\n']);
%% Visualisation
time = N*(Lt+Ld)/fs;

max_h = max(h(:));
max_H = max(H_db(:));
min_H = min(H_db(:));
figure('Name','visualisation of the demodulation')
subplot(2,2,2);colormap(colorMap); image(imageData); axis image; title('The transmitted image'); 
for i = 1:size(h,2) 
    subplot(2,2,1)
    plot(h(:,i))
    title('Channel in time domain')
    axis([0 160 -max_h max_h])
    
    subplot(2,2,3)
    %(2:N/2-1)*fs/2/(N/2-1), 
    plot((1:N/2-1)*fs/2/(N/2-1),H_db(:,i));
    title('Channel in frequency domain')
    axis([-inf, inf, min_H max_H]);
    ylabel('Magnitude')
    xlabel('Frequency')
    
    if i*sum(freq_bins)*K*Ld<length(rxBitStream)
        imageRx = bitstreamtoimage(rxBitStream(1:i*sum(freq_bins)*K*Ld), imageSize, bitsPerPixel);
    else
        imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);
    end
    subplot(2,2,4); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;
    pause(time)
end    
waitforbuttonpress;

%%% met bitloading
%% trainingsblock
n_trainblock    = (N/2-1)*K;   % contains N/2-1 random QAM symbols
trainblock_bit  = randi([0 1], 1, n_trainblock);
trainblock_qam  = qam_mod(trainblock_bit, K);

%% H-channel voor freq bins 
seq = randi([0 1], (N/2-1)*K, 1); % random bits
trainblock = qam_mod(seq, K); 
trainblocks = repmat(trainblock,30,1);
Tx = ofdm_mod_tb(trainblocks, N, cp_size);


[pulse, ~] = sinusoid(1000, 1, 1, fs);
[simin,nbsecs,fs,pulse ] = initparams(Tx,fs,L, pulse);

sim('recplay'); 
rec = simout.signals.values;

[Rx_2] = alignIO(rec,pulse,L);
Rx_2 = Rx_2(1:length(Tx));

[seq_qam , H] = ofdm_demod_tb(Rx_2, N, cp_size, trainblock);

freq_bins = ofdm_freq_bins(H, N, threshold);
% freq_bins = ones(N/2-1,1); 

%% Modulate (QAM amd OFDM)
qamStream_in = qam_mod(bitStream, K);
Tx = ofdm_mod(qamStream_in, N, cp_size, freq_bins, trainblock_qam, Lt, Ld);

%% Play en record
tic;
[pulse, ~] = sinusoid(1000, 1, 1, fs);
[simin,nbsecs,fs,pulse ] = initparams(Tx,fs,L, pulse);

sim('recplay'); 
rec = simout.signals.values;
toc;
%% Align

[Rx] = alignIO(rec,pulse,L);
Rx = Rx(1:length(Tx)+L+10);

%% Demodulate (QAM amd OFDM)
[rxQamStream, H ] = ofdm_demod(Rx, N, cp_size,freq_bins, trainblock_qam, Lt, Ld);
rxQamStream = rxQamStream(1:length(qamStream_in));
rxBitStream = qam_demod(rxQamStream,K);

%% Channel freq response
H_fft = [zeros(1,size(H,2)) ; H ; zeros(1,size(H,2));flipud(conj(H))];
h = ifft(H_fft,N);
H_db = 20*log10(abs(H)); % niet gebruikte freq bins op 0

freq_not = ~freq_bins;
H_not = 20*log10(abs(H))./(freq_not);
%% BER
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


fprintf(['With on-off bit loading: Ber =', num2str(berTransmission),'\n']);
%% Visualisation
time = N*(Lt+Ld)/fs;

max_h = max(h(:));
max_H = max(H_db(:));
min_H = min(H_db(:));
figure('Name','visualisation of the demodulation')
subplot(2,2,2);colormap(colorMap); image(imageData); axis image; title('The transmitted image'); 
for i = 1:size(h,2) 
    subplot(2,2,1)
    plot(h(:,i))
    title('Channel in time domain')
    axis([0 160 -max_h max_h])
    
    subplot(2,2,3)
    hold on
    plot((1:N/2-1)*fs/2/(N/2-1),H_db(:,i),'b-', 'Linewidth', 0.75);
    plot((1:N/2-1)*fs/2/(N/2-1),H_not(:,i),'r-', 'Linewidth', 1);
    title('Channel in frequency domain')
    axis([-inf, inf, min_H max_H]);
    ylabel('Magnitude')
    xlabel('Frequency')
    hold off
    
    if i*sum(freq_bins)*K*Ld<length(rxBitStream)
        imageRx = bitstreamtoimage(rxBitStream(1:i*sum(freq_bins)*K*Ld), imageSize, bitsPerPixel);
    else
        imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);
    end
    subplot(2,2,4); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;
    pause(time);
end   

%% vraagjes

%%% ALIGN 
% Indien niet goed ge-alligned (meer als cp_size verschil) kan er verkeerd data
% geinterpreteerd worden en QAM-symbolen/bits verschoven zijn.
% Data kan verloren gaan indien niet mooi ergens in cp_size beginnend

%%% TOEPLITZ
% FREQ OF TIJD
% freq want computationeel interessanter 
% + geven zelfde resultaat
% + freq heeft minder geheugen nodig WANT matrix puntsgewijs delen door een
% getal (geen tweede matrix nodig)

%%% HAND VOOR MICROFOON
% ber stijgt

%%% INVLOED VAN FS
% hoe meer samples, hoe nauwkeuriger (oversampled *whispered)
% maar dan duurt het ook langer om alles te berekenen

%%% VISUALIZE DEMOD, TIME VARIATION + SHAPE OF CHANNEL
% kanaal verandert vrij snel in de tijd

%%% QUALITY AND DRAWBACK OF ONOFF BITLOADING
% betere kwaliteit aangezien enkel goede freq gebruikt worden
% helaas duurt het langer om door te sturen (minder freq gebruikt = deel van BB dus meer tijd)

%%% PILOTTONES VOORDEEL
% kanaal wordt constant afgeschat en is dus nauwkeuriger
% robuust

%%% MICROFOON BEWEGEN POLIT VS TRAINBLOCK
% pilot is beter want het nieuwe kanaal meteen afgeschat



