function [seq_ofdm] = ofdm_mod_ada(stream_in, N, cp_size, freq_bins)

M = sum(freq_bins); % M = aantal ellementen die we in 1 frame kunnen steken

P = ceil(size(stream_in, 1) / M); % P berekenen : totale lengte van seq_qam
% gedeeld door het aamtal elementen dat we in elke frame kunnen steken. in
% elke frame passen N getallen in totaal waarvan 2 nul zijn (dc) en de
% helft van de rest de comlex toegevoegde. dus N/2-1 elementen van de seq
% in 1 frame.

% P -> moet een geheel getal zijn dus afronden naar boven en dan de seq padden met nullen

%Nullen toevoegen zodat we kunnen reshapen
seq_qam_pad = padarray(stream_in,P*M-size(stream_in, 1),0,'post');


frames = zeros(N,P);
k = 1;

for j = 1: length(freq_bins)
    if  freq_bins(j)
        K = freq_bins(j);
        seq = seq_qam_pad(k:k+P*K-1);
        k = k+P*K;
        seq_cluster = reshape(seq, K, length(seq)/K).';                     
        seq_int = bi2de(seq_cluster);                                     
        seq_mod = qammod(seq_int, 2^K, 'bin', 'UnitAveragePower', true); 
        frames(j+1,:)= seq_mod; 
        
        frames(N+2-j-1, :) = conj(seq_mod);
       
    end
end

seq_ifft = ifft(frames); % ifft van nemen
seq_ifft = [seq_ifft(end-cp_size+1 : end, :); seq_ifft ]; % add cyclic prefix, so total packet size = P+size(cyclic_prefix)
seq_ofdm = seq_ifft(:); % Serialize data stream

end

