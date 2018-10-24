function [seq_demod] = qam_demod(seq_mod, K)
%%INPUT

%%OUTPUT
    
%% Code
seq_demod_int = qamdemod(seq_mod, 2^K, 'bin', 'UnitAveragePower', true);    % demodulate: gemoduleerde seq weer omzetten naar symbolen (= decimale getallen) 
seq_demod_cluster = de2bi(seq_demod_int).';                                 % de symbolen terug omzetten naar bitsequenties (in rijen DUS transpose)
seq_demod = seq_demod_cluster(:);                                           % maak vector van matrix

end
