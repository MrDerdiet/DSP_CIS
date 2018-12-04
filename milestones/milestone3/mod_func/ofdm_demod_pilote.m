function [seq_qam, H ] = ofdm_demod(seq_ifft, N, cp_size,freq_bins, pilote)

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

% pilote tones eruit
pilotes = 1:2:N/2-1;
tones = seq_fft(pilotes,:);

% channelest
c_est = tones./pilote;
H_est = (c_est(1:end-1,:) + c_est(2:end,:))/2; % interpolate linear

% data eruit
seq_fft(pilotes,:) = [];
freq_bins(pilotes,:) = [];

% channel eq
seq_qam = seq_fft./H_est;

% freq bins eruit
del_row = find(~freq_bins); %neem de index van alle rijen zonder info 
seq_qam(del_row,:) = [];

seq_qam = seq_qam(:); % terug mooi een vector van maken 

% H en channelest mergen
H = zeros(N/2-1,P);
H(pilotes,:) = c_est; 
H(pilotes(2:end)-1,:) = H_est; 
end

