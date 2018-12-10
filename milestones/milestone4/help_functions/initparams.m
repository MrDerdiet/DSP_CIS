function [simin,nbsecs,fs,pulse ] = initparams(toplay,fs,L, pulse)
%%INPUT
% toplay: vector (nx1) that contains the samples of an audio signal with
% fs: samplerate

if (nargin < 3)
    L = 0;
    pulse = [];
end



if max(abs(toplay)) > 0
       toplay = toplay/max(abs(toplay)); % normalize input to '1' (avoid clipping/increase DR)
end 

if max(abs(pulse)) > 0
       pulse = pulse/max(abs(pulse)); % normalize input to '1' (avoid clipping/increase DR)
end
if size(toplay,2) <2
    simin = [zeros(2*fs, 2); pulse zeros(length(pulse),1) ; zeros(L,2); toplay zeros(length(toplay),1); zeros(fs, 2)];
else
    simin = [zeros(2*fs, 2); pulse pulse ; zeros(L,2); toplay ; zeros(fs, 2)];
end
nbsecs = length(simin)/fs;

end

