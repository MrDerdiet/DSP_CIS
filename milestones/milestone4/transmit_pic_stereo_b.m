clearvars; hold on; close all;
addpath(genpath('help_functions'), genpath('data'), genpath('mod_func'));

%% Convert BMP image to bitstream
[bitStream_in, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');

%% Parameters
fs = 16000;
N = 512; 
K = 4; 
Lt = 2;
Ld = 8;
BW = 50;
threshold = BW/100;
cp_size = 180;
L = 160; %channel length
sig_time = 0.5;

[pulse1, ~] = sinusoid(1000, 1, .3, fs);
[pulse2, ~] = sinusoid(1500, 1, .3, fs);

%% kanaal est
seq = randi([0 1], (N/2-1)*K, 1); % random bits
trainblock = qam_mod(seq, K); 
trainblocks = repmat(trainblock,20,1);
Tx = ofdm_mod_tb(trainblocks, N, cp_size);

% test beide kant

sig1 =  Tx;
sig2 =  zeros(length(Tx), 1);
sig_in1 = [sig1 sig2];
sig_in2 = [ sig2 sig1 ];

[simin ,nbsecs,fs,pulse1,pulse2  ] = initparams_stereo(sig_in1, sig_in2,fs,L, pulse1,pulse2);

sim('recplay'); 
rec = simout.signals.values;

[Rx_1] = alignIO(rec,pulse1,L);
Rx_1 = Rx_1(1:length(Tx));%%

[Rx_2] = alignIO(rec,pulse2,L);
Rx_2 = Rx_2(1:length(Tx));%%

[~ , H_1] = ofdm_demod_tb(Rx_1, N, cp_size, trainblock);
[~ , H_2] = ofdm_demod_tb(Rx_2, N, cp_size, trainblock);

down = min(20*log10(abs(H_1+H_2)));
up = max(20*log10(abs(H_1+H_2)));
    
figure;
    hold on
    plot( (1:N/2-1) *fs/(N/2-1)/2, 20*log10(abs(H_1)), 'b-');
    plot( (1:N/2-1) *fs/(N/2-1)/2, 20*log10(abs(H_2)), 'r-');
    plot( (1:N/2-1) *fs/(N/2-1)/2, 20*log10(abs(H_1+H_2)), 'g-');
    title('Channel in frequency domain_H 1')

    axis([-inf, inf, down, up]);
    ylabel('Magnitude')
    xlabel('Frequency')
    legend('H1','H2','H1+2');


freq_bins = ofdm_freq_bins(H_1+H_2, N, threshold);


%% bereken a en b

[a,b] = fixed_transmitter_side_beamformer(H_1,H_2);

% test voor voorwaarde
% controle = a.*conj(a)+b.*conj(b);
%% Generate trainblock, known by both transmitter and receiver
n_trainblock    = (N/2-1)*K;   % contains N/2-1 random QAM symbols
trainblock_bit  = randi([0 1], 1, n_trainblock);
trainblock_qam  = qam_mod(trainblock_bit, K);

%% Modulate
qamStream_in = qam_mod(bitStream_in, K);
Tx = ofdm_mod_stereo(qamStream_in, N, cp_size, freq_bins, trainblock_qam, a, b , Lt, Ld);


%% Play and record audio
[simin, nbsecs, fs, pulse] = initparams(Tx,fs,L, pulse1);
sim('recplay')
out=simout.signals.values;


%% Align in- and output
Rx = alignIO(out, pulse,L);
rest = mod(length(Rx),(N + cp_size));
Rx = Rx(1: end-rest);


%% Demodulate 
[qamStream_out_pad, H_out] = ofdm_demod(Rx, N, cp_size, freq_bins, trainblock_qam, Lt, Ld);
qamStream_out = qamStream_out_pad(1: length(qamStream_in));
bitStream_out = qam_demod(qamStream_out, K);

H_out_full = [zeros(1,size(H_out, 2)); H_out; zeros(1,size(H_out, 2)); flipud(conj(H_out))];

%% BER
[ber_IO, error_vector] = ber(bitStream_in, bitStream_out);


%% Output image
% imageRx = bitstreamtoimage(bitStream_out, imageSize, bitsPerPixel);
% colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;

%% Visualize
time = N*(Lt+Ld)/fs;
h = ifft(H_out_full);



figure('Name','visualization of the demodulation','units','normalized','outerposition',[0 0 1 1])
subplot(2,2,2);colormap(colorMap); image(imageData); axis image; title('The transmitted image'); 
for i = 1:size(h,2)-1
    subplot(2,2,1)
    plot(h(1:L-1,i))
    title('Channel in time domain')
    axis([0 L -0.5*max(abs(h(:))) 0.5*max(abs(h(:)))])
    
    subplot(2,2,3)
    plot( (2:N/2-1) *fs/(N/2-1)/2, 20*log10(abs(H_out_full(2:N/2-1,i))));
    title('Channel in frequency domain')
    axis([0, fs/2, -60 10]);
    ylabel('Magnitude')
    xlabel('Frequency')
    
    if i*sum(freq_bins)*K*Ld<length(bitStream_out)
        imageRx = bitstreamtoimage(bitStream_out(1:i*sum(freq_bins)*K*Ld), imageSize, bitsPerPixel);
    else
        imageRx = bitstreamtoimage(bitStream_out, imageSize, bitsPerPixel);
    end
    subplot(2,2,4); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;
    pause(time)
end    