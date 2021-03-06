function [a,b] = fixed_transmitter_side_beamformer(H_1,H_2)
H_1 = [0; H_1;0;flipud(conj(H_1))];
H_2 = [0; H_2;0;flipud(conj(H_2))];
numerator= sqrt(H_1.*conj(H_1)+H_2.*conj(H_2));
a = conj(H_1)./numerator;
b = conj(H_2)./numerator;
end

