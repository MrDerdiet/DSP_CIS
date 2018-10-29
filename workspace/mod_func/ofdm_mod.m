function [seq_ofdm] = ofdm_mod(seq_qam, N, cp_size)


P = ceil(size(seq_qam, 1) / (N/2 -1)); % P berekenen : totale lengte van seq_qam
% gedeeld door het aamtal elementen dat we in elke frame kunnen steken. in
% elke frame passen N getallen in totaal waarvan 2 nul zijn (dc) en de
% helft van de rest de comlex toegevoegde. dus N/2-1 elementen van de seq
% in 1 frame.

% P -> moet een geheel getal zijn dus afronden naar boven en dan de seq padden met nullen

%Nullem toevoegen zodat we kunnen reshapen
seq_qam_pad = padarray(seq_qam,P*(N/2-1)-size(seq_qam, 1),0,'post');

% data groeperen in P x (N/2-1) matrix Group data in in P packets (X[1], ..., X[N/2-1] in each)
seq_qam_reshaped = reshape(seq_qam_pad, N/2-1, P);

% De frames maken, 0 vam boven, de seq, 0 rij, toegevoegd compex  
frames = [zeros(1, P); seq_qam_reshaped; zeros(1,P); flipud(conj(seq_qam_reshaped))];

% 
seq_ifft = ifft(frames); % ifft van nemen
seq_ifft = [seq_ifft(end-cp_size+1 : end, :); seq_ifft ]; % add cyclic prefix, so total packet size = P+size(cyclic_prefix)
seq_ofdm = seq_ifft(:); % Serialize data stream


end

