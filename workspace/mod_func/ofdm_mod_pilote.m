function [seq_ofdm, dim_P] = ofdm_mod_pilote(seq_qam, N, cp_size, freq_bins, pilote, Lt, Ld)
%% updated for ex 6 

% Bepaal de P 
M = sum(freq_bins); % M = aantal ellementen die we in 1 frame kunnen steken
P = ceil(size(seq_qam, 1) / M); % P berekenen : totale lengte van seq_qam

%pad data als nodig
seq_qam_pad = padarray(seq_qam,P*M-size(seq_qam, 1),0,'post');

% data groeperen in P x M matrix Group data in in P packets (X[1], ..., X[N/2-1] in each)
seq_qam_reshaped = reshape(seq_qam_pad, M, P);

% gewone frames maken
frames_data = zeros(N,P);

i = 1;
for j = 1: length(freq_bins)
    if  freq_bins(j)
        frames_data(j+1,:)= seq_qam_reshaped(i,:); 
        frames_data(N+2-j-1, :) = conj(seq_qam_reshaped(i, :));
        i = i+1;
    end
end

%% trainingblokken maken
n_tb = ceil(P/Ld);
dim_P = P+ n_tb*Lt;
frames = zeros(N, dim_P);

trainblocks = repmat(trainblock, 1, Lt);
trainblock_frames = [zeros(1, Lt); trainblocks; zeros(1,Lt); flipud(conj(trainblocks))];

% begin met tb en dan data => voor ld datablock, steek er een tb in
for i = 1:n_tb -1
    frames(:,1+(Lt+Ld)*(i-1):(Lt+Ld)*(i-1)+Lt) = trainblock_frames;
    frames(:,1+(Lt+Ld)*(i-1)+Lt:(Lt+Ld)*(i)) = frames_data(:, 1+(i-1)*Ld: i*Ld);
end
i = i+1;
% Laatste appart doen (niet altijd de Ld lengte)
frames(:,1+(Lt+Ld)*(i-1):(Lt+Ld)*(i-1)+Lt) = trainblock_frames;
frames(:,1+(Lt+Ld)*(i-1)+Lt:end) = frames_data(:,1+(i-1)*Ld: end);

% 
seq_ifft = ifft(frames); % ifft van nemen
seq_ifft = [seq_ifft(end-cp_size+1 : end, :); seq_ifft ]; % add cyclic prefix, so total packet size = P+size(cyclic_prefix)
seq_ofdm = seq_ifft(:); % Serialize data stream


end

