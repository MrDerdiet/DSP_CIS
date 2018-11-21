addpath(genpath('helper_functions'), genpath('mod_func'), genpath('data'));
%clearvars; hold on; close all;

%% PARAMETERS
K = 6 ; % Number of bits per constellation
M = 2^K;
n = 256*16*K; % Number of data points
E_tot = 2/3*(M-1);
plot_normalized = 1; %If 1, scatterplots will show normalized values.
SNR = 35; 


%% (3.1) 
%%% 3.3.1 %%%
seq = randi([0 1], n, 1); % Generate a pseudo random binary sequence of a user defined length.

%%% 3.1.2 %%%
seq_mod = qam_mod(seq, K);

%%% 3.1.3 %%%
scatterplot(seq_mod*(sqrt(E_tot))^(1-plot_normalized), 1, 0, 'b*'); % constellation diagram

%%% 3.1.4 %%%
%%% Normalization to unit power necessary?
%   -> Yes, if no normalization: larger QAM => larger constellation => much larger total power transmitted
%       => the constellation points need to be closed to each other, for the same transmitted power.
%       = normalizing the constellation points      
%   If yes, dif normalization factors for constellation sizes?
%   -> Normalization factor = sqrt(total energy in constellation)
%                           = sqrt( 2/3*(M-1) );
%                           ?????????????????????????????????????????????????????????????????????????????????????????????
% 2/3*(M-1) = cumsum(abs(seq_mod(i))^2) for(i= 1:n/K) 
% Veronderstelling: random data: dus elke conctellatie even waarschijnlijk.

%%% 3.1.5 %%%                                          
seq_mod_noise = awgn(seq_mod, SNR, 'measured');                             % zet ruis op de gemoduleerde sequentie (dus seq + wgn)
scatterplot(seq_mod_noise*(sqrt(E_tot))^(1-plot_normalized), 1, 0, 'b*');   % constellation diagram
% afhankelijk van SNR kan een grotere K genomen worden
% SNR groot: ruis is klein DUS waarden hebben nog groot onderscheid --> grote K kan genomen worden
% SNR klein: ruis is groot DUS waarden komen dicht bij elkaar te liggen --> K moet beperkt worden om het onderscheid te behouden

%%% 3.1.6 %%%
seq_demod = qam_demod(seq_mod_noise, K);

%%% 3.1.7 %%%
[ber, seq_error] = ber(seq, seq_demod);
% SNR omlaag => BER omhoog
% SNR omhoog => BER omlaag
% K omhoog => hogere SNR nodig voor zelfde BER
% K omlaag => lagere SNR nodig voor zelfde BER

