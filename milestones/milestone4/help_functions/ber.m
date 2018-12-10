function [ber, seq_error] = ber(seq_in, seq_out)
% BER - Calcule Bit Error Rate (BER) in communication channel from input to output
%% INPUT
%   'seq_in'    - data sequence (m x n) corresponding to input 
%   'seq_out'   - data sequence (m x n) corresponding to output
% 

%% OUTPUT
%   'ber'       - Calculated BER from input to output
%   'seq_error' - Vector corresponding to error positions (=different values between 'seq_in' and 'seq_out')
%               - length

%% Parameter check
if (nargin < 2)
    return;
end
if (size(seq_in) ~= size(seq_out))
    error('seq_in and seq_out should have the same dimensions');
end

%% BER calculation
seq_error = (seq_in ~= seq_out);    % equality test
num_errors = sum(sum(seq_error));   % #errors
ber = num_errors / length(seq_in);  % ber = #errors / #data

end