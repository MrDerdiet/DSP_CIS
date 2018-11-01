addpath(genpath('helper_functions'), genpath('mod_func'), genpath('data'));
clearvars; hold on; close all;

%% PARAMETERS
K = 4 ; % Number of bits per constellation
M = 2^K;
N = 32;
n = 256*32*K; % Number of data points
fs = 16000;
SNR = 2.5*7;

%% 3.2.1
% mirror operation: nodig voor de negatieve frequenties \\ -> opzoeken 

%% 3.2.3

cp_size = 2; %cyclic prefix size

seq = randi([0 1], n, 1); % create random sequence of n bits
%SNR = 20;   % Channel SNR, [3.2.5]

seq_qam = qam_mod(seq, K);  % 2^K-QAM mapping

seq_ofdm = ofdm_mod(seq_qam, N, cp_size);

seq_res_qam = ofdm_demod(seq_ofdm, N, cp_size, 1);

seq_res = qam_demod(seq_res_qam, K);
seq_res = seq_res(1:n); %in ofdm_mod worden er nullen toegevoegd, deze worden nu weer weggelaten -> seq_res even lang maken als oorspronkelijke freq

ber_ideal = ber(seq, seq_res); % ber alleen nemen van hoe lang de oorspronkelijke seq was (toegevoegde nullen weglaten)

%% 3.2.4
% done

%% 3.2.5
% expression
R = K * (N/2-1) * fs/2;
% bit/sec = bit/qam * qam/frame * frames/sec
% fs/2: hoe snel frames na elkaar sturen zonder overlap = bandbreedte (-> Nyquist: fs=2*B dus B=fs/2)

%% 3.2.6
seq_noise_ofdm = awgn(seq_ofdm, SNR, 'measured'); % ruis op OFDM signaal voegen, dit komt aan de input

seq_noise_res_qam = ofdm_demod(seq_noise_ofdm, N, cp_size, 1);

%scatterplot(seq_noise_res_qam(:));
seq_noise_result = qam_demod(seq_noise_res_qam, K);

seq_noise_result = seq_noise_result(1:n);
ber_noise = ber(seq, seq_noise_result);



