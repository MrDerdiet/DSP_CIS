addpath(genpath('helper_functions'), genpath('mod_func'), genpath('data'));
clearvars; hold on; close all;
% Exercise session 4: DMT-OFDM transmission scheme

K = 8;
N = 512;

SNR = 60;
L = 160; %channel length
cp_size = L+16;

%% Convert BMP image to bitstream
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image_extra.bmp');

%% QAM modulation
qamStream = qam_mod(bitStream,K);



%% Channel -> with fftfilt

% h = randi(10000, L, 1)./5000 - 1;                   % create random sequence of L numbers (matrix met L rijen en 1 kolom --> kolomvector)
% h = [h; zeros(N -size(h, 1), 1)];                   % aanvullen met nullen WANT zie slide 34
%                                                     %(mbv fft van aangelengde vector de diagonaalelem bepalen)

% h van IR2
%%%%%parameters
fs=16000;
h_order = 160;                                      %(bepaald bij IR1.m)  
sig_time = 2;

%%%%%Generate input signal (witte ruis), play and record
                                                    %(om het signaal uit de output te kunnen halen)
sig =  wgn(fs*sig_time, 1, 0);
[simin, nbsecs, fs] = initparams(sig, fs);
sim('recplay'); 
rec = simout.signals.values;

start = find(rec > 0.3*max(rec), 1, 'first')-20;    %geeft de eerst index die hieraan voldoet
                                                    %start is de index waar het signaal (ruis) groter is dan de ruis die er
                                                    %altijd is --> zeker dat het signaal daar start

y = rec(start+h_order :start+sig_time*fs+h_order);  %geselecteerde output
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

[freq_bins] = ofdm_freq_bins(Hn_vector, N, 0.5);



%% OFDM modulation
ofdmStream = ofdm_mod_onoff(qamStream, N, cp_size,freq_bins);

%% channel
ofdmStream = fftfilt(h, ofdmStream);                % alsof het signaal over het kanaal wordt gestuurd (signaal*TF)
rxOfdmStream = awgn(ofdmStream, SNR, 'measured');   % witte ruis erop

%% OFDM demodulation
rxQamStream = ofdm_demod_onoff(rxOfdmStream, N, cp_size, Hn_vector,freq_bins);

%% QAM demodulation
rxBitStream = qam_demod(rxQamStream,K);

%% Compute BER
rxBitStream = rxBitStream(1:length(bitStream));
berTransmission = ber(bitStream,rxBitStream);


%% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);

%% Plot images
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;




%% ON-OF bitloading

% -> errors bij elkaar in de buurt: lijnen over het beeld


