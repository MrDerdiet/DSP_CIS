function [seq_qam, H ] = ofdm_demod(seq_ifft, N, cp_size,freq_bins, trainblock, Lt, Ld)

% Bepaal de P 
M = sum(freq_bins); % M = aantal ellementen die we in 1 frame kunnen steken
P = ceil(length(seq_ifft) / (N+cp_size)); % P berekenen : totale lengte van seq_qam

% pad seq_ifft
seq_ifft = padarray(seq_ifft,P*(N+cp_size)-size(seq_ifft, 1),0,'post'); % dan kan dit goed gereshaped worden

% bepaal lengte van data 
data_length = P - ceil(P/(Lt+Ld))*Lt; % P - #trainingblocks*length_trainingblock = lengte data (de verschillende blokken opgeteld)

%reshape P en remove cp
seq_reshaped = reshape(seq_ifft, N + cp_size,P );
seq_reshaped = seq_reshaped(cp_size+1 : end , :);% Prefix eraf halen



% fft ervan
seq_fft = fft(seq_reshaped);
seq_fft = seq_fft((2:N/2), :); % nuttige data

% Nul matrix toevoegen om for loop makkelijker
padded = Lt+Ld-rem(P,Lt+Ld); % #toegevoegde kolommen
seq_fft = [seq_fft, zeros(N/2-1,padded)];

% tb maken 
trainblocks = repmat(trainblock, 1, Lt);

H = [];

H_temp = zeros(N/2-1,1);

seq_qam = []; 

for i=1:(Ld+Lt):P+1
    % gegevens eruit
    train_rx = seq_fft(:,(i : i + Lt-1));
    data = seq_fft(:, (i +(Lt): i+(Ld+Lt)-1));
    % Channel estimation
    for j= 1: (N/2-1)
        H_temp(j)= mrdivide(train_rx(j,:),trainblocks(j,:));
    end
    % Channel equal
    data_temp = data./H_temp;
    seq_qam = [seq_qam, data_temp];
    H = [H ,H_temp];
end

% data eruit halen
del_row = find(~freq_bins); %neem de index van alle rijen zonder info 

seq_qam(del_row,:) = []; 

% Nulmatrix eraf halen + trasponeren
seq_qam = seq_qam(:,1:end-padded-1);

seq_qam = seq_qam.';
seq_qam = seq_qam(:); % terug mooi een vector van maken 

end

