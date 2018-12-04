function [seq_qam, H ] = ofdm_demod(seq_ifft, N, cp_size, tb)


% bereken aantal pakketjes (P frames van N+cp lenghte)
P = ceil(size(seq_ifft, 1)/(N + cp_size));

seq_ifft = padarray(seq_ifft,P*(N + cp_size)-size(seq_ifft, 1),0,'post');

% de vector terug mooi omvormen naar een matrix met elke kollom een frame
seq_reshaped = reshape(seq_ifft, N + cp_size, P);
seq_reshaped = seq_reshaped(cp_size+1 : end , :);% Prefix eraf halen

seq_fft = fft(seq_reshaped); % FFT ervan nemen

tb = repmat(tb,1,P);

% 'Channel equalization'

trainblock_rx = seq_fft(2 : N/2,:); % enkel de sectie eruithalen met nuttige data -> complex toegevoegdes weglaten

H = zeros(N/2-1,1);

for j= 1: (N/2-1)
    H(j)= mrdivide(trainblock_rx(j,:),tb(j,:));
end 

seq_qam = trainblock_rx./H;

seq_qam = seq_qam(:); % terug mooi een vector van maken 

end

