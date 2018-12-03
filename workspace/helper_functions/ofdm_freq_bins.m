function [freq_bins] = ofdm_freq_bins(Hn, N, threshold)


Hn_abs = abs(Hn); % Neem de Hn van elke frequentie (negatieve weg en DC ook weg)

[out,idx] = sort(Hn_abs); 

comp_H = out(ceil((1-threshold+0.0000001)*length(Hn)));
freq_bins = ( Hn_abs >= comp_H).*1; % Only keep bins within certain range from max frequency response (eg 50% = -6dB).
  
end

