function [seq_mod] = qam_mod(seq, K)
%%INPUT
%%OUTPUT


%% real_code
seq_cluster = reshape(seq, length(seq)/K, K);   % Group into sequences of N bits
seq_int = bi2de(seq_cluster);                   % Convert each goup to single symbol (=integer)
seq_mod = qammod(seq_int, 2^K, 'bin', 'UnitAveragePower', true);          % Modulate 

end


