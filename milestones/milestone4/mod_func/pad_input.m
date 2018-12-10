% function [data_pad] = pad_input(seq_in, K, num_bins)
function [data_pad, pad_bits] = pad_input(seq_in, modulus)
% modulus = K*num_bins
    rest = mod(length(seq_in), modulus);
    pad_bits = modulus-rest;
    if rest ~= 0
        data_pad = padarray(seq_in, pad_bits, 'post');
    else
        data_pad = seq_in;
    end
end
