function [seq_out] = ofdm_demod_ada(seq_ifft, N, cp_size, H_channel,freq_bins)

% bereken aantal pakketjes (P frames van N+cp lenghte)
P = size(seq_ifft, 1)/(N + cp_size);

% de vector terug mooi omvormen naar een matrix met elke kollom een frame
seq_reshaped = reshape(seq_ifft, N + cp_size, P);
seq_reshaped = seq_reshaped(cp_size+1 : end , :);% Prefix eraf halen

seq_fft = fft(seq_reshaped); % FFT ervan nemen

% 'Channel equalization'
%seq_qam = mldivide(H_channel, seq_fft); %karakterisktieken van kanaal er terug uit halen
                                        %ga ervan uit dat H_channel een diagonale matrix is zoals slide 33 les 3
                                        %dit is schalen met de inverse (delen door de elem)
frames = seq_fft./H_channel;



seq_out = [];
% enkel de sectie eruithalen met nuttige data -> complex toegevoegdes weglaten

for j = 1: length(freq_bins)
    if  freq_bins(j)
        K = freq_bins(j);
        seq_mod = frames(j+1,:);
        seq_demod_int = qamdemod(seq_mod, 2^K, 'bin', 'UnitAveragePower', true);    % demodulate: gemoduleerde seq weer omzetten naar symbolen (= decimale getallen) 
        seq_demod_cluster = de2bi(seq_demod_int).';                                 % de symbolen terug omzetten naar bitsequenties (in rijen DUS transpose)
        seq_demod = seq_demod_cluster(:);
        seq_out = [seq_out ; seq_demod];
    end
end 


end

