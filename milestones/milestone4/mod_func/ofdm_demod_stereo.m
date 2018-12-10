function [seq_qam] = ofdm_demod_stereo(Rx, N, cp_size, freq_bins, H_1, H_2)
%%INPUT
%%OUTPUT

%% real_code
P = ceil(length(Rx) / (N + cp_size));
seq_reshaped = reshape(Rx, N + cp_size, P);
seq_reshaped = seq_reshaped(cp_size+1 : end , :); % Prefix eraf halen
seq_fft = fft(seq_reshaped);

num_bins = length(freq_bins);

seq_qam_reshaped = seq_fft ./  sqrt(H_1.*conj(H_1)+H_2.*conj(H_2));

seq_qam = zeros(num_bins, P);
for i = 1: num_bins
    seq_qam(i,:) = seq_qam_reshaped(freq_bins(i), :);
end

seq_qam = seq_qam(:);
end




