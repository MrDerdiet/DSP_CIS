function [seq_qam_vector, W_tot] = ofdm_demod_adaptive(Rx, N, cp_size, freq_bins, trainblock, Lt, K, alpha, mu)
%%INPUT
% trainblock = exact known trainblock (1 packet)
% n_tb       = #trainblocks per H_channel estimation
% delta_tb   = #packets between 2 H_channel estimations

%%OUTPUT

%% real_code

M = sum(freq_bins); % M = aantal ellementen die we in 1 frame kunnen steken

%% reshape
P = ceil(length(Rx) / (N + cp_size)); %Total #packets
seq_reshaped = reshape(Rx, N + cp_size, P);
seq_reshaped = seq_reshaped(cp_size+1 : end , :); % Take out CP

%% data eruit
seq_fft = fft(seq_reshaped);
seq_fft = seq_fft((2:N/2), :);

del_row = find(~freq_bins);
seq_fft(del_row,:) = []; % alle data die nul zijn (ongebruikte freq) eraf
trainblock(del_row,:) = []; 


%% channel est

trainblock_tx_frames = repmat(trainblock, 1, Lt);
trainblock_rx_frames = seq_fft(:,(1: Lt));
H = zeros(M,1);
for k= 1: M
        H(k)= mrdivide(trainblock_rx_frames(k,:), trainblock_tx_frames(k,:));
end 

% Isolate data
data_rx_frames = seq_fft(:, (Lt+1: end));

% Add zero because of padded data
qam_symbols = [0; qammod((0:2^K-1)',2^K,'bin','UnitAveragePower', true)]; % generate the different possible qam symboles

% data is zonder Dc -> alle frequencies beginnen eentje lager 



W = zeros(M, P-Lt);
seq_qam = zeros(M,P-Lt);

for j = 1: M
    W(j, 1) =  1/conj(H(j));
end

% calculate next W for each frame, frequency by frequency, extract usefull
% data, save only W and the symbole
for i= 1: P-Lt
    for j = 1: M
        Xtilde_temp = conj(W(j,i))*data_rx_frames(j, i);
        [~, index] = min(abs(qam_symbols-Xtilde_temp));
        Xcirconflex_temp = qam_symbols(index);
        seq_qam(j,i) = Xcirconflex_temp;  
        W(j,i+1) =  W(j,i) + mu/(alpha+conj(data_rx_frames(j, i))*data_rx_frames(j, i))*data_rx_frames(j, i)*conj(Xcirconflex_temp-Xtilde_temp);
    end
end
% Delete zero packets (padding) in caller func

seq_qam_vector = seq_qam(:);
W_tot = zeros(N/2-1, P-Lt+1);
insert_row = find(freq_bins);
W_tot(insert_row,:) = W;
end




