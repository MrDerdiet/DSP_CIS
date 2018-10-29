function [seq_qam] = ofdm_demod(seq_ifft, N, cp_size, H_channel)

% bereken aantal pakketjes (P frames van N+cp lenghte)
P = size(seq_ifft, 1)/(N + cp_size);

% de vector terug mooi omvormen naar een matrix met elke collum een frame
seq_reshaped = reshape(seq_ifft, N + cp_size, P);
seq_reshaped = seq_reshaped(cp_size+1 : end , :);% Prefix eraf halen

seq_fft = fft(seq_reshaped); % FFT ervan nemen

% 'Channel equalization'
seq_qam = mldivide(H_channel, seq_fft);

seq_qam = seq_qam(2 : N/2, 1:P); % enkel de sectie eruithalen met nuttige data

seq_qam = seq_qam(:); % terug mooi een vector van maken 
end

