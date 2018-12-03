function [seq_ofdm] = ofdm_mod_pilote(seq_qam, N, cp_size, freq_bins, pilote)
%% updated for ex 6 

% Bepaal de P 
M = sum(freq_bins); % M = aantal ellementen die we in 1 frame kunnen steken
P = ceil(size(seq_qam, 1) / M); % P berekenen : totale lengte van seq_qam

%pad data als nodig
seq_qam_pad = padarray(seq_qam,P*M-size(seq_qam, 1),0,'post');

% data groeperen in P x M matrix Group data in in P packets (X[1], ..., X[N/2-1] in each)
seq_qam_reshaped = reshape(seq_qam_pad, M, P);

% pilotes maken
pilote = repmat(pilote,1,P);
% gewone frames maken
frames_data = zeros(N/2-1,P);

% Data erin
data_i = find(freq_bins == 1);
frames_data(data_i,:) = seq_qam_reshaped;

pilotes = 1:2:N/2-1;
frames_data(pilotes,:) = pilote;
frames = [zeros(1, P); frames_data; zeros(1,P); flipud(conj(frames_data))];


% 
seq_ifft = ifft(frames); % ifft van nemen
seq_ifft = [seq_ifft(end-cp_size+1 : end, :); seq_ifft ]; % add cyclic prefix, so total packet size = P+size(cyclic_prefix)
seq_ofdm = seq_ifft(:); % Serialize data stream


end

