function [freq_bins] = ofdm_freq_bins(Hn, N, threshold)


Hn_abs = abs(Hn(2: N/2)); % Neem de Hn van elke frequentie (negatieve weg en DC ook weg)

max_H = max(Hn_abs);
freq_bins = ( Hn_abs >= threshold*max_H).*1; % Only keep bins within certain range from max frequency response (eg 50% = -6dB).
  
end

