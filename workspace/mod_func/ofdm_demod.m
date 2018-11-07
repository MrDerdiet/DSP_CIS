function [seq_qam] = ofdm_demod(seq_ifft, N, cp_size, H_channel)

% bereken aantal pakketjes (P frames van N+cp lenghte)
P = size(seq_ifft, 1)/(N + cp_size);

% de vector terug mooi omvormen naar een matrix met elke collum een frame
seq_reshaped = reshape(seq_ifft, N + cp_size, P);
seq_reshaped = seq_reshaped(cp_size+1 : end , :);% Prefix eraf halen

seq_fft = fft(seq_reshaped); % FFT ervan nemen

% 'Channel equalization'
seq_qam = mldivide(H_channel, seq_fft); %karakterisktieken van kanaal er terug uit halen
                                        %ga ervan uit dat H_channel een diagonale matrix is zoals slide 33 les 3
                                        %dit is schalen met de inverse (delen door de elem)
%%%% ??????????? DIT KAN MISSCHIEN ELEMENTGEWIJS EN DAN H_CHANNEL ALS VECTOR IPV DIAGMATRIX ??????????????????????????

seq_qam = seq_qam(2 : N/2, 1:P); % enkel de sectie eruithalen met nuttige data -> complex toegevoegdes weglaten

seq_qam = seq_qam(:); % terug mooi een vector van maken 
end

