function [freq_bins] = ofdm_freq_bins(h_channel_init, N, bit_loading, fs, threshold, sort_bins)


   % FFT, throws away complex conjugates (negative frequencies) + abs

H_channel_abs = H_channel_abs(2: end);                  % Delete DC frequency
if (bit_loading)
    if (~sort_bins)
        max_H = max(H_channel_abs);                            % Find bin with largest frequency response magnitude
        freq_bins = find( H_channel_abs > threshold*max_H);    % Only keep bins within certain range from max frequency response (eg 50% = -6dB).
    else                                                       % Keep only (threshold*100)% best frequency bins (regardless of relative difference)
        H_sort = sort(H_channel_abs, 'descend');
        H_threshold = H_sort(floor(threshold*length(H_channel_abs)));
        freq_bins = find(H_channel_abs > H_threshold); 
    end
else
    freq_bins = [1: N/2-1];
end

freq_bins = freq_bins +1;       %DC already deleted
end

