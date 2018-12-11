function [a,b] = fixed_transmitter_side_beamformer(H_1,H_2)
H_1 = [0; H_1;0;flipud(conj(H_1))];
H_2 = [0; H_2;0;flipud(conj(H_2))];
numerator= sqrt(H_1.*conj(H_1)+H_2.*conj(H_2));
a = conj(H_1)./numerator;
b = conj(H_2)./numerator;
a(1) = 0;
b(1) = 0;
N = length(a)/2+1;
a(N) = 0;
b(N) = 0;
end

