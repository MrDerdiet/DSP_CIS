addpath(genpath('helper_functions'), genpath('mod_func'), genpath('data'));
clearvars; hold on; close all;
% Exercise session 4: DMT-OFDM transmission scheme

K = 6;
N = 512;

SNR = 35;
L = 160; %channel length
cp_size = L+16;
fs=16000;
h_order = 160; 
dftsize = N;
window = ones(dftsize, 1);
noverlap = 0;  

%% Convert BMP image to bitstream
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');


%% 1 zend ruis; bereken p(k) 

sig_time = 3;

% silence is send
sig1 = zeros(sig_time*fs, 1);
[simin,nbsecs,fs] = initparams(sig1,fs);
sim('recplay');
rec1=simout.signals.values;
[sig_fft1, df_sig1, t_sig1, psd_sig1] = spectrogram(rec1, window, noverlap, dftsize, fs);


% signal, white noise to get all frequencies
sig2 = wgn(sig_time*fs,1,0); %White noise
[simin,nbsecs,fs] = initparams(sig2,fs);
sim('recplay');
rec2=simout.signals.values;
[sig_fft2, df_sig2, t_sig2, psd_sig2] = spectrogram(rec2, window, noverlap, dftsize, fs);

p_noise = abs(mean(psd_sig1,2))';
p_signal = abs(mean(psd_sig2-psd_sig1,2))';
%ongecorelleerde PSD, dus mag je gewoon optellen of aftrekken

figure('Name', 'PSD');
subplot(2,1,1);
plot(df_sig1, 10*log10(p_noise));  title( 'noise PSD' ); xlabel( 'frequency(Hz)' ); ylabel( 'magnitude (dB)' );
subplot(2,1,2);
plot(df_sig2, 10*log10(p_signal)); title( 'signal PSD' ); xlabel( 'frequency(Hz)' ); ylabel( 'magnitude (dB)' );

N = dftsize/2;
ratio = p_signal./p_noise;
c_chan = sum(fs/(2*N)*log2(1+ratio));


%% QAM modulation
qamStream = qam_mod(bitStream,K);



%% Channel -> with fftfilt

% h = randi(10000, L, 1)./5000 - 1;                   % create random sequence of L numbers (matrix met L rijen en 1 kolom --> kolomvector)
% h = [h; zeros(N -size(h, 1), 1)];                   % aanvullen met nullen WANT zie slide 34
%                                                     %(mbv fft van aangelengde vector de diagonaalelem bepalen)

% h van IR2
%%%%%parameters
                                     %(bepaald bij IR1.m)  

%%%%%Generate input signal (witte ruis), play and record (zie vorig deel)
                                                    %(om het signaal uit de output te kunnen halen)

start = find(rec2 > 0.3*max(rec2), 1, 'first')-20;    %geeft de eerst index die hieraan voldoet
                                                    %start is de index waar het signaal (ruis) groter is dan de ruis die er
                                                    %altijd is --> zeker dat het signaal daar start

y = rec2(start+h_order :start+sig_time*fs+h_order);  %geselecteerde output
u = simin(2*fs+1+h_order:(2+sig_time)*fs+h_order+1);%input: neem de input vanaf waar je het
                                                    %signaal opstart (2*fs+1) en tot het einde van het
                                                    %signaal((2+sig_time)*fs)+ h_order want het sijpelt zoveel na
r = flipud(simin(2*fs+1 : 2*fs+h_order+1,1)).';     %moet geflipt en neen gewoon h_order owv definitie toeplitz elememten
U = toeplitz(u.', r);

%%%%%least squares 
h = mldivide(U, y);                                 %NOTE: least-square solutions in same plane -> x = A^(-1)*b 
h = [h; zeros(N -size(h, 1), 1)];

Hn_vector = fft(h);
%Hn_matrix = diag(Hn_vector);                        % mbv deze matrix kan je dan doen zoals slide 33 (fft_van_yk = Hn_elem*fft_van_xk)


%% b(k)
%[freq_bins] = ofdm_freq_bins(Hn_vector, N, 0);

p_noise = p_noise(2:N/2);
p_signal = (10*log10(p_signal));
p_signal = p_signal(2:N/2);

Hk = Hn_vector(2:N/2);
Hk2 = abs(Hk).^2;
bins = floor(log2(1+Hk2./(10*p_noise)));
bins = floor((bins-min(bins))/(max(bins)-min(bins))*K);

freq_bins = p_signal-min(p_signal);
threshold = 0.33;
freq_bins = (freq_bins >= threshold*max(freq_bins)).* freq_bins;
freq_bins = floor(freq_bins / (max(freq_bins)-min(freq_bins)) *K);

figure('Name', 'Bits_adaptive');
subplot(2,1,1);
plot(bins);title( 'Transfer function' ); xlabel( 'frequency_bin(Hz)' ); ylabel( 'bits (K)' );
subplot(2,1,2);
plot(freq_bins);title( 'Shannon measurement' ); xlabel( 'frequency_bin(Hz)' ); ylabel( 'bits (K)' );

%% OFDM modulation
ofdmStream = ofdm_mod_ada(bitStream, N, cp_size,bins);

%% channel
ofdmStream = fftfilt(h, ofdmStream);                % alsof het signaal over het kanaal wordt gestuurd (signaal*TF)
rxOfdmStream = awgn(ofdmStream, SNR, 'measured');   % witte ruis erop

%% OFDM demodulation
rxBitStream = ofdm_demod_ada(rxOfdmStream, N, cp_size, Hn_vector,bins);

%% QAM demodulation
%%rxBitStream = qam_demod(rxQamStream,K);

%% Compute BER
rxBitStream = rxBitStream(1:length(bitStream));
berTransmission = ber(bitStream,rxBitStream);


%% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);

%% Plot images
figure;
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;


%% ON-OF bitloading

% -> errors bij elkaar in de buurt: lijnen over het beeld
% bepaalde frequenties die slecht zijn gekozen


