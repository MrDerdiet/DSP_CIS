function [simin,nbsecs,fs,pulse1,pulse2 ] = initparams_stereo(toplay1, toplay2,fs,L, pulse1,pulse2)
%%INPUT
% toplay: vector (nx1) that contains the samples of an audio signal with
% fs: samplerate



if max(abs(toplay1)) > 0
       toplay1 = toplay1/max(abs(toplay1)); % normalize input to '1' (avoid clipping/increase DR)
end 
if max(abs(toplay2)) > 0
       toplay2 = toplay2/max(abs(toplay2)); % normalize input to '1' (avoid clipping/increase DR)
end 

if max(abs(pulse1)) > 0
       pulse1 = pulse1/max(abs(pulse1)); % normalize input to '1' (avoid clipping/increase DR)
end

if max(abs(pulse2)) > 0
       pulse2 = pulse2/max(abs(pulse2)); % normalize input to '1' (avoid clipping/increase DR)
end

simin = [zeros(2*fs, 2); pulse1 pulse1 ; zeros(L,2); toplay1; zeros(2*L,2);pulse2 pulse2; toplay2; zeros(fs, 2)];

nbsecs = length(simin)/fs;
end

