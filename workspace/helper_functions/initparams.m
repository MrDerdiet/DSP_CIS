function [simin,nbsecs,fs] = initparams(toplay,fs)
%%INPUT
% toplay: vector (nx1) that contains the samples of an audio signal with
% fs: samplerate

toplay = toplay/max(abs(toplay)); % normalize input to '1' (avoid clipping/increase DR)

simin = [zeros(2*fs, 2); toplay toplay; zeros(fs, 2)];

nbsecs = length(simin)/fs;

end

