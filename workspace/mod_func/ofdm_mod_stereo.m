function [ofdmStream1 ,ofdmStream2] = ofdm_mod_stereo(seq_qam, N, cp_size, freq_bins, trainblock_qam, a, b, Lt, Ld)

%% DESCRIPTION
%% INPUT
% NOTE: In order to work: N > size(freq_bins);
%% OUTPUT

%% real_code
delta_sf = Lt+Ld;

% Group data in P packets
num_bins = length(freq_bins); 
% P = n_data_packets
P = length(seq_qam) / num_bins; 
seq_qam_reshaped = reshape(seq_qam, num_bins, P);

data_frames = zeros(N, P);
for i = 1: P
    for j = 1: num_bins
        data_frames(freq_bins(j), i) = seq_qam_reshaped(j, i);
        data_frames(N+2-freq_bins(j), i) = conj(seq_qam_reshaped(j, i));
    end
end

% Ld = delta_tb -Lt;
n_packets = ceil(P/Ld);
frames = zeros(N, P+n_packets*Lt);

% Make trainblock frame
trainblocks = repmat(trainblock_qam, 1, Lt);
trainblock_frames = [zeros(1, Lt); trainblocks; zeros(1,Lt); flipud(conj(trainblocks))];

% data_sf = data_superframe = Ld*data_frame
% tb_sf = trainblock_superframe = Lt*trainblock_frames
% sf = superframe = data_superframe + trainblock_superframe  

for i= 1: n_packets-1
   %volledige superframes
   frames(:, 1+delta_sf*(i-1): delta_sf*(i-1)+Lt) = trainblock_frames;
   frames(:, 1+delta_sf*(i-1)+Lt: delta_sf*i)  = data_frames(:, 1+(i-1)*Ld: i*Ld);
end
% laatste (evt onvolledig) superframe
frames(:, 1+delta_sf*(i): delta_sf*(i)+Lt) = trainblock_frames;
frames(:, 1+delta_sf*(i)+Lt: end)  = data_frames(:, 1+(i)*Ld: end);

% Total power = 1
frames1 = a.*frames;
frames2 = b.*frames;

ofdmStream1_ifft = ifft(frames1);
ofdmStream2_ifft = ifft(frames2);

% Cyclic Prefix
ofdmStream1_cp = [ofdmStream1_ifft(end-cp_size+1:end,:) ; ofdmStream1_ifft];
ofdmStream2_cp = [ofdmStream2_ifft(end-cp_size+1:end,:) ; ofdmStream2_ifft];

% Serialize
ofdmStream1 = ofdmStream1_cp(:);
ofdmStream2 = ofdmStream2_cp(:);
end


