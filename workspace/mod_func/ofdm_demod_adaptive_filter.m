function [seq_qam, W_tot ] = ofdm_demod_adaptive_filter(seq_ifft, N, cp_size,freq_bins, trainblock, Lt, mu, K,alpha)

qam_symbols = qammod((0:2^K-1)',2^K,'bin','UnitAveragePower', true); % generate the different possible qam symboles
% Bepaal de P 
M = sum(freq_bins); % M = aantal ellementen die we in 1 frame kunnen steken
P = ceil(length(seq_ifft) / (N+cp_size)); % P berekenen : totale lengte van seq_qam

% pad seq_ifft
seq_ifft = padarray(seq_ifft,P*(N+cp_size)-size(seq_ifft, 1),0,'post'); % dan kan dit goed gereshaped worden

%reshape P en remove cp
seq_reshaped = reshape(seq_ifft, N + cp_size,P );
seq_reshaped = seq_reshaped(cp_size+1 : end , :);% Prefix eraf halen

% fft ervan
seq_fft = fft(seq_reshaped);
seq_fft = seq_fft((2:N/2), :); % nuttige data
del_row = find(~freq_bins);
seq_fft(del_row,:) = []; % alle data die nul zijn (ongebruikte freq) eraf
trainblock(del_row,:) = []; 
% channel est
tb = seq_fft(:, 1:Lt);
H_init = tb./repmat(trainblock,1,Lt);
H_init = mean(H_init,2);


% data eruit

data_rx_frames = seq_fft(:, (Lt+1: end)); % trainblock eraf halen 
W = zeros(M, P-Lt);
seq_qam = zeros(M, P-Lt);
W(:,1) = 1./conj(H_init);

%% calculate new W each frame

for i=1:P-Lt
    for k =  1:M 
        Xtilde_temp = conj(W(k,i))*data_rx_frames(k, i);
        [~, index] = min(abs(qam_symbols-Xtilde_temp));
        Xcirconflex_temp = qam_symbols(index);
        seq_qam(k,i) = Xcirconflex_temp;
        W(k,i+1) =  W(k,i) + mu/(alpha+conj(data_rx_frames(k, i))*data_rx_frames(k, i))*data_rx_frames(k, i)*conj(Xcirconflex_temp-Xtilde_temp);
    end
end    
W_tot = zeros(N/2-1, P-Lt+1);
insert_row = find(freq_bins);
W_tot(insert_row,:) = W;

seq_qam = seq_qam(:); % terug mooi een vector van maken 

end

