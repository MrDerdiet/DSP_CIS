function [seq_demod] = qam_demod(seq_mod, K)
%%INPUT

%%OUTPUT

%% Parameters
if nargin < 2
    padded = 0;
end
    

%% Code
seq_demod_int = qamdemod(seq_mod, 2^K, 'bin', 'UnitAveragePower', true);          % Modulate 
seq_demod_cluster = de2bi(seq_demod_int);
seq_demod = seq_demod_cluster(:); 

end


