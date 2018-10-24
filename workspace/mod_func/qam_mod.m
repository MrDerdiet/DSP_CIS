function [seq_mod] = qam_mod(seq, K)
%%INPUT
%%OUTPUT

%% real_code
seq_cluster = reshape(seq, K, length(seq)/K).';                     % seq omvormen naar matrix (rijen volgen elkaar op in de seq)) = Group into sequences of N bits
seq_int = bi2de(seq_cluster);                                       % Elke rij wordt een decimaal getal (single symbol (=integer))
M = 2^K;
seq_mod = qammod(seq_int, M, 'bin', 'UnitAveragePower', true);      % Elk symbool een cmloexe waarde toekennen = Modulate 
                                                                    % 'bin': natural binary-coded ordering 
                                                                    % 'UnitAveragePower': logical scalar value, scales to average power 1 mbv 'true'
                                                                   

end