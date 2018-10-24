addpath(genpath('helper_functions'));
clearvars; hold on; close all;

%% White noise experiment
fs = 16000;
dftsize = 256;

sig = wgn(fs*2,1,0); %White noise

[simin,nbsecs,fs] = initparams(sig,fs);

sim('recplay');

rec=simout.signals.values;


window = ones(dftsize, 1);          % Rectangular window
noverlap = 16;                      % N samples to overlap between segments
                                    % overlap geeft nauwkeurigheid

figure('Name', 'spectograms');      % figure 1
subplot( 2, 1, 1 ); 
spectrogram(sig, window, noverlap, dftsize, fs, 'yaxis'); 
[sig_fft, df_sig, t_sig, psd_sig] = spectrogram(sig, window, noverlap, dftsize, fs);
title( 'original signal spectogram' );

subplot( 2, 1, 2 );
spectrogram(rec, window, noverlap, dftsize, fs, 'yaxis');
[rec_fft, df_rec, t_rec, psd_rec] = spectrogram(rec, window, noverlap, dftsize, fs);
title( 'recorded signal spectogram' );

psd_sig_avg = mean(psd_sig, 2)';
psd_rec_avg = mean(psd_rec, 2)';

figure('Name', 'PSD');              % figure 2
subplot(2,1, 1);
plot(df_sig, 10*log10(psd_sig_avg));  title( 'original signal PSD' ); xlabel( 'frequency(Hz)' ); ylabel( 'magnitude (dB)' );
subplot(2,1,2);
plot(df_rec, 10*log10(psd_rec_avg)); title( 'recorded signal PSD' ); xlabel( 'frequency(Hz)' ); ylabel( 'magnitude (dB)' );

%% IR2
h_order = 160;
sig_time = 2;
start = find(rec > 0.3*max(rec), 1, 'first')-20;    %geeft de eerst index die hieraan voldoet
                                                    %start is de index waar het signaal (ruis) groter is dan de ruis die er
                                                    %altijd is --> zeker dat het signaal daar start

y = rec(start+h_order :start+sig_time*fs+h_order);  %geselecteerde output
u = simin(2*fs+1+h_order:(2+sig_time)*fs+h_order+1);%input: neem de input vanaf waar je het
                                                    %signaal opstart (2*fs+1) en tot het einde van het
                                                    %signaal((2+sig_time)*fs)+ h_order want het sijpelt zoveel na
r = flipud(simin(2*fs+1 : 2*fs+h_order+1,1)).';     %moet geflipt en neem gewoon h_order owv definitie toeplitz elementen
U = toeplitz(u.', r);

%%%%%least squares 
h = mldivide(U, y);                                 %NOTE: least-square solutions in same plane -> x = A^(-1)*b 

%%%%%plot
t = (0: 1 : (length(h)-1)) /fs;                     % is een rijvector
fft_h = fft(h)';                                    % transpose want we willen een rijvector om te kunnen plotten
length_h = round((length(fft_h)-1)/2 );             % Take half the length
fft_h_trunc = abs(fft_h(1 : length_h));             % Use only first half of fft samples (=positive frequencies).
df_h = (0 : fs /(2*length_h) : (fs-1)/2);           % Frequency vector for x-axis in plots

figure('Name', 'IR2');
subplot(2, 1, 1);
    plot(t, h);
    title('impulse response');
    xlabel('time (sec)');
    ylabel('magnitude');
subplot(2, 1, 2);
    plot(df_h, 20*log10(fft_h_trunc));
    title('frequency response');
    xlabel('frequency(Hz)');
    ylabel('magnitude (dB)');

%% IR1
sig_pulse = zeros(1600,1); sig_pulse(1) = 1;        % impulse signal 
[simin, nbsecs, fs] = initparams(sig_pulse, fs);

sim('recplay')
rec_pulse=simout.signals.values;


impulse = find(rec_pulse==max(rec_pulse));          % vind de index (plaats van het sample) waar de waarde van de rec max wordt, dit is het begin van de impulsrespons
IR = rec_pulse(impulse-10:impulse+150);             % ietsje voor die index beginnen en ver genoeg erna eindigen om de hele respons te houden

 t = (0: 1 : (length(IR)-1)) /fs;                   % is een rijvector
 
% fft
fft_IR = fft(IR)';                                  % transpose want we willen een rijvector om te kunnen plotten
length_IR = round((length(fft_IR)-1)/2 );           % Take half the length
fft_IR_trunc = abs(fft_IR(1 : length_IR));          % Use only first half of fft samples (=positive frequencies).
df_IR = (0 : fs /(2*length_IR) : (fs-1)/2);         % Frequency vector for x-axis in plots



figure('Name', 'IR1');
subplot(2, 1, 1);
    plot(t, IR);
    title('impulse response');
    xlabel('time (sec)');
    ylabel('magnitude');
subplot(2, 1, 2);
    plot(df_IR, 20*log10(fft_IR_trunc));
    title('frequency response');
    xlabel('frequency(Hz)');
    ylabel('magnitude (dB)');

