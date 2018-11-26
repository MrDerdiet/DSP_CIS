function [ber, seq_error] = ber(seq_1, seq_2)


if length(seq_1) > length(seq_2)
    seq_1 = seq_1(1:length(seq_2));
    warning('shorted seq_1');
end

if length(seq_1) < length(seq_2)
    seq_2 = seq_2(1:length(seq_1));
    warning('shorted seq_2');
end
        

%% BER calculation
seq_error = (seq_1 ~= seq_2);       % equality test (fout=1, juist=0)
num_errors = sum(sum(seq_error));   % #errors optellen
ber = num_errors / length(seq_1);   % ber = #errors / #data (lenght seq_1 = seq_2)

end

