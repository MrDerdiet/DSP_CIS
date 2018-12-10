function [ofdmStream1 ,ofdmStream2] = ofdm_mod_stereo(seq_qam, N, cp_size, freq_bins, trainblock, a, b, Lt, Ld)

%% DESCRIPTION
%% INPUT
% NOTE: In order to work: N > size(freq_bins);
%% OUTPUT
% Bepaal de P 
M = sum(freq_bins); % M = aantal ellementen die we in 1 frame kunnen steken
P = ceil(size(seq_qam, 1) / M); % P berekenen : totale lengte van seq_qam

%pad data als nodig
seq_qam_pad = padarray(seq_qam,P*M-size(seq_qam, 1),0,'post');

% data groeperen in P x M matrix Group data in in P packets (X[1], ..., X[N/2-1] in each)
seq_qam_reshaped = reshape(seq_qam_pad, M, P);

% gewone frames maken
frames_data_temp = zeros(N/2-1,P);


data_i = find(freq_bins == 1);
frames_data_temp(data_i,:) = seq_qam_reshaped;

frames_data = [zeros(1, P); frames_data_temp; zeros(1,P); flipud(conj(frames_data_temp))];

%% trainingblokken maken
n_tb = ceil(P/Ld);
dim_P = P+n_tb*Lt;
frames = zeros(N, dim_P);


trainblocks = repmat(trainblock, 1, Lt);

trainblock_frames = [zeros(1, Lt); trainblocks; zeros(1,Lt); flipud(conj(trainblocks))];

% begin met tb en dan data => voor ld datablock, steek er een tb in
for i = 1:n_tb -1
    frames(:,1+(Lt+Ld)*(i-1):(Lt+Ld)*(i-1)+Lt) = trainblock_frames;
    frames(:,1+(Lt+Ld)*(i-1)+Lt:(Lt+Ld)*(i)) = frames_data(:, 1+(i-1)*Ld: i*Ld);
end

i = i+1 ;

% Laatste appart doen (niet altijd de Ld lengte)
frames(:,1+(Lt+Ld)*(i-1):(Lt+Ld)*(i-1)+Lt) = trainblock_frames;
frames(:,1+(Lt+Ld)*(i-1)+Lt:end) = frames_data(:, 1+(i-1)*Ld: end);

% Total power = 1
frames1 = a.*frames;
frames2 = b.*frames;

seq1_ifft = ifft(frames1);
seq2_ifft = ifft(frames2);

% Cyclic Prefix
seq1_ifft = [seq1_ifft(end-cp_size+1:end,:) ; seq1_ifft];
seq2_ifft = [seq2_ifft(end-cp_size+1:end,:) ; seq2_ifft];

% Serialize
ofdmStream1 = seq1_ifft(:);
ofdmStream2 = seq2_ifft(:);
end


