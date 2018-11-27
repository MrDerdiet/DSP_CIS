function [seq_qam, H ] = ofdm_demod(seq_ifft, N, cp_size,freq_bins, trainblock, Lt, Ld)

% Bepaal de P 
M = sum(freq_bins); % M = aantal ellementen die we in 1 frame kunnen steken
P = ceil(length(seq_ifft) / (N+cp_size)); % P berekenen : totale lengte van seq_qam

% bepaal lengte van data 
data_length = P -(floor(P/(Lt+Ld))+1)*(Lt); 

%reshape P en remove cp
seq_reshaped = reshape(seq_ifft, N + cp_size, P);
seq_reshaped = seq_reshaped(cp_size+1 : end , :);% Prefix eraf halen

% Nul matrix toevoegen om for loop makkelijker 

% fft ervan
seq_fft = fft(seq_reshaped);
seq_fft = seq_fft((2:N/2), :); % nuttige data

% tb maken 
trainblocks = repmat(trainblock, 1, Lt);

H = [];

H_temp = zeros(N/2-1,1);

seq_qam = []; % Nakijken

for i=1:(Ld+Lt):P
    % gegevens eruit
    train_rx = seq_fft(:,(i : i + Lt-1));
    data = seq_fft(:, (i +(Lt): i+delta_tb-1));
    % Channel estimation
    for j= 1: (N/2-1)
        H_temp(j)= mrdivide(train_rx(j,:),trainblocks(j,:));
    end
    % Channel equal
    data_temp = (H_temp.').\data_block;
    seq_qam = [seq_qam, data_temp];
    H = [H,H_temp.'];
end

% data eruit halen
del_row = find(~freq_bins); %neem de index van alle rijen zonder info 

seq_qam(del_row,:) = []; 

% Nulmatrix eraf halen 
 todo


seq_qam = seq_qam(:); % terug mooi een vector van maken 

end

