function [ber, seq_error] = ber(seq_1, seq_2)

if (size(seq_1) ~= size(seq_2))
    error('seq_1 and seq_2 should have the same dimensions');
end

%% BER calculation
seq_error = (seq_1 ~= seq_2);       % equality test (fout=1, juist=0)
num_errors = sum(sum(seq_error));   % #errors optellen
ber = num_errors / length(seq_1);   % ber = #errors / #data (lenght seq_1 = seq_2)

end

