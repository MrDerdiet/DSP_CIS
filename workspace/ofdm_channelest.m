addpath(genpath('helper_functions'), genpath('mod_func'), genpath('data'));
clearvars; hold on; close all;
% Exercise session 4: DMT-OFDM transmission scheme

K = 6;
N = 512;
fs = 16000;

SNR = 60;
L = 160; %channel length
cp_size = L+16;

%% 5.1
seq = randi([0 1], (N/2-1)*K, 1); % random bits
trainblock = qam_mod(seq, K); 
trainblocks = repmat(trainblock,100,1);
Tx = ofdm_mod_tb(trainblocks, N, cp_size);


load('h_ir2','h');
h = [zeros(20,1); h; zeros(N -20 -size(h, 1), 1)];

Rx = fftfilt(h, Tx);

[seq_qam , H] = ofdm_demod_tb(Rx, N, cp_size, trainblock);

receive_deqam = qam_demod(seq_qam, K); 
trainblocks_deqam = qam_demod(trainblocks, K); 

[bert, bert_seq] = ber(receive_deqam,trainblocks_deqam);

t = (0: 1 : (length(h)-1)) /fs;
fft_h = fft(h)';                                    % transpose want we willen een rijvector om te kunnen plotten
length_h = round((length(fft_h)-1)/2) -1;           % Take half the length
fft_h_trunc = abs(fft_h(1 : length_h));
df_h = (0 : fs /(2*length_h) : (fs-1)/2);

figure('Name', 'H accoustic');
subplot(2, 1, 1);
    plot(t, h);
    title('acoustic impulse response');
    xlabel('time (sec)');
    ylabel('magnitude');
subplot(2, 1, 2);
    plot(df_h, 20*log10(fft_h_trunc));
    title('acoustic channel frequency response');
    xlabel('frequency(Hz)');
    ylabel('magnitude (dB)');


h_est = (real(ifft(H)'))';
h_est = [h_est; zeros(N -size(h_est, 1), 1)];
figure('Name', 'H estimated');
subplot(2, 1, 1);
    plot(t, h_est);
    title('estimated impulse response');
    xlabel('time (sec)');
    ylabel('magnitude');
subplot(2, 1, 2);
    plot(df_h, real(20*log10(H)));
    title('estimated channel frequency response');
    xlabel('frequency(Hz)');
    ylabel('magnitude (dB)');
     

% h_error = h - h_est;
% h_log_error = 20*log10(fft_h_trunc) - (real(20*log10(H)))';
% figure('Name', 'Error estimate-acoustic');
% subplot(2, 1, 1);
%     plot(t, h_error);
%     title('acoustic error');
%     xlabel('time (sec)');
%     ylabel('magnitude');
% subplot(2, 1, 2);
%     plot(df_h, h_log_error);
%     title('acoustic channel frequency response');
%     xlabel('frequency(Hz)');
%     ylabel('magnitude (dB)');
    
%% 5.2 
% BER is vrij hoog, we gaan voort hopelijk komt dit doorat we niet
% optimaliseren
seq = randi([0 1], (N/2-1)*K, 1); % random bits
trainblock = qam_mod(seq, K); 
trainblocks = repmat(trainblock,100,1);
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



[bert, bert_seq] = ber(receive_deqam,trainblocks_deqam);

h_est = (real(ifft(H)'))';
h_est = [h_est; zeros(N -size(h_est, 1), 1)];
figure('Name', 'H estimated over channel');
subplot(2, 1, 1);
    plot(t, h_est);
    title('estimated impulse response');
    xlabel('time (sec)');
    ylabel('magnitude');
subplot(2, 1, 2);
    plot(df_h, real(20*log10(H)));
    title('estimated channel frequency response');
    xlabel('frequency(Hz)');
    ylabel('magnitude (dB)');
    
    
 







   
   

  