function [sig, t] = sinusoid(freq, A, time, fs)
if nargin < 4
    error('Input parameters');
end

ts = 1/fs;
t = (0: ts: time)';


sig = 0;
for i = 1 : length(freq)
    sig = sig + A*sin(2*pi*freq(i)*t);
end


end

