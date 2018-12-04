function [out_aligned] = alignIO(out,pulse,L)
%% Nodig 
% Nodig om impulse response te beperken in lengte + we willen enkel de data
% output 

%% Alternatief userdefined threshold -> nadelen
% afhankelijk van SNR (als kleine snr moet threshold hoger; maar we weten snr niet op voorhand dus trail and error)
% Niet robuust tegen over luide ruis pieken (of luide  buren tijdens het project) 

%%

l_pulse = length(pulse);

% calculate & find max of cross-correlation
[c, lags] = xcorr(out, pulse);
[~, I] = max(c);
start_seq = lags(I);

% gooi stilte eraf
out_sync = out(start_seq+1: end); 

% delete pulse + silence IR -20
out_aligned = out_sync(l_pulse-20+L: end);

end

