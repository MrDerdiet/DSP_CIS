addpath(genpath('helper_functions'), genpath('mod_func'), genpath('data'));
clearvars; hold on; close all;

%% 7.1 
K = 4;
n = 1000;
n_seq    = n*K;   % contains N/2-1 random QAM symbols
seq  = randi([0 1], 1, n_seq);
seq_qam  = qam_mod(seq, K);


H = (0.666 + 0.666j);


%% 7.2 Implement the adaptive filter

delta = (0.5+0.1j); % initiele afwijking
mu = [0.05;0.25; 0.5; 5]; %5 % bepaalt hoeveel rekening je houdt met de fout ~stapgrootte
alpha = 10;


W = zeros(n+1, length(mu)); %initialisatie
W(1,:) = 1/conj(H)+delta; % w = 1/Hk* +delta
W_error = zeros(n, length(mu)); % initialisatie
Xtilde = zeros(n, length(mu)); % initialisatie
Xcirconflex = zeros(n, length(mu)); % initialisatie
legend_names = cell.empty(length(mu), 0); % naam van de legendes in de grafieken die straks geplot worden


Y = H*seq_qam; % kanaal op seq
for i=1:length(mu)
    for k=1:n
        Xtilde(k,i) = conj(W(k,i))*Y(k);
        Xtildedeqam = qamdemod(Xtilde(k,i), 2^K, 'bin', 'UnitAveragePower', true); % Decision Device 
        Xcirconflex(k,i) = qammod(Xtildedeqam, 2^K, 'bin', 'UnitAveragePower', true);
        W(k+1,i) =  W(k,i) + mu(i)/(alpha+conj(Y(k))*Y(k))*Y(k)*conj(Xcirconflex(k,i)-Xtilde(k,i)); % NLMS formule
        W_error(k,i) = W(k,i) - 1/conj(H);
    end
    legend_names{i} = strcat('Mu = ',num2str(mu(i)));
end
error = Xcirconflex - Xtilde;


figure('Name', 'W_error signal (constant Hk)');
subplot(3,1,1);
plot(real(W_error)); title('Real')
subplot(3,1,2)
plot(imag(W_error)); title('imag');
subplot(3,1,3)
plot(abs(W_error)); title('abs');
legend(legend_names);

figure('Name', 'QAM error signal (constant Hk)');
subplot(3,1,1);
plot(real(error)); title('Real')
subplot(3,1,2)
plot(imag(error)); title('imag');
subplot(3,1,3)
plot(abs(error)); title('abs');
legend(legend_names);

%% 7.1.5

H = zeros(n+1, 1);
H(1) = (0.666 + 0.666j);
for i=1:length(mu)
    for k=1:n
        Y(k) = H(k)*seq_qam(k);
        Xtilde(k,i) = conj(W(k,i))*Y(k); 
        Xtildedeqam = qamdemod(Xtilde(k,i), 2^K, 'bin', 'UnitAveragePower', true); % Decision Device 
        Xcirconflex(k,i) = qammod(Xtildedeqam, 2^K, 'bin', 'UnitAveragePower', true);
        W(k+1,i) =  W(k,i) + mu(i)/(alpha+conj(Y(k))*Y(k))*Y(k)*conj(Xcirconflex(k,i)-Xtilde(k,i));
        W_error(k,i) = W(k,i) -1/conj(H(k));
        H(k+1) = H(k) + (randn(1) + randn(1)*1j)/100; %100 is snel varierende H, 500 is traag varierende H
    end
    legend_names{i} = strcat('Mu = ',num2str(mu(i)));
end
error = Xcirconflex - Xtilde;

figure('Name', 'W_error signal (variable Hk)');
subplot(3,1,1);
plot(real(W_error)); title('Real')
subplot(3,1,2)
plot(imag(W_error)); title('imag');
subplot(3,1,3)
plot(abs(W_error)); title('abs');
legend(legend_names);

figure('Name', 'QAM Error (variable Hk)');
subplot(3,1,1);
plot(real(error)); title('Real')
subplot(3,1,2)
plot(imag(error)); title('imag');
subplot(3,1,3)
plot(abs(error)); title('abs');
legend(legend_names);

figure('Name', 'Hk');
subplot(3,1,1);
plot(real(H)); title('Real')
subplot(3,1,2)
plot(imag(H)); title('imag');
subplot(3,1,3)
plot(abs(H)); title('abs');

%%% CONCLUSIE
% Kies stapgrootte zo dat het de H geod kan tracken, niet te snel, niet te
% traag
% snel varierende H -> grootte stapgrootte (anders zal je H nooit inhalen)
% traag varierende H -> kleine stapgrootte
% let op: te kleine stapgrootte is nooit goed









