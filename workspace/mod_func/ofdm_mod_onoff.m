function [seq_ofdm] = ofdm_mod_onoff(seq_qam, N, cp_size, freq_bins)

M = sum(freq_bins); % M = aantal ellementen die we in 1 frame kunnen steken

P = ceil(size(seq_qam, 1) / M); % P berekenen : totale lengte van seq_qam
% gedeeld door het aamtal elementen dat we in elke frame kunnen steken. in
% elke frame passen N getallen in totaal waarvan 2 nul zijn (dc) en de
% helft van de rest de comlex toegevoegde. dus N/2-1 elementen van de seq
% in 1 frame.

% P -> moet een geheel getal zijn dus afronden naar boven en dan de seq padden met nullen

%Nullen toevoegen zodat we kunnen reshapen
seq_qam_pad = padarray(seq_qam,P*M-size(seq_qam, 1),0,'post');

% data groeperen in P x (N/2-1) matrix Group data in in P packets (X[1], ..., X[N/2-1] in each)
seq_qam_reshaped = reshape(seq_qam_pad, M, P);


frames = zeros(N,P);
i = 1;
for j = 1: length(freq_bins)
    if  freq_bins(j)
        frames(j+1,:)= seq_qam_reshaped(i,:); 
        frames(N+2-j-1, :) = conj(seq_qam_reshaped(i, :));
        i = i+1;
    end
end

seq_ifft = ifft(frames); % ifft van nemen
seq_ifft = [seq_ifft(end-cp_size+1 : end, :); seq_ifft ]; % add cyclic prefix, so total packet size = P+size(cyclic_prefix)
seq_ofdm = seq_ifft(:); % Serialize data stream

end

