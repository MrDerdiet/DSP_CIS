addpath(genpath('helper_functions'));
clearvars; hold on; close all; 

%% 2.1.1 Create and plot impulse responce

fs = 16000;
sig = zeros(1600,1); sig(1) = 1;                    % impulse signal 
[simin, nbsecs, fs] = initparams(sig, fs);

sim('recplay')
rec=simout.signals.values;


impulse = find(rec==max(rec));                      % vind de index (plaats van het samle) waar de waarde van de rec max wordt, dit is het begin van de impulsrespons
IR = rec(impulse-10:impulse+150);                   % ietsje voor die index beginnen en ver genoeg erna eindigen om de hele respons te houden

 t = (0: 1 : (length(IR)-1)) /fs;                   % is een rijvector
 
% fft
fft_IR = fft(IR)';                                  % transpose want we willen een rijvector om te kunnen plotten
length_IR = round((length(fft_IR)-1)/2 );           % Take half the length
fft_IR_trunc = abs(fft_IR(1 : length_IR));          % Use only first half of fft samples (=positive frequencies).
df_IR = (0 : fs /(2*length_IR) : (fs-1)/2);         % Frequency vector for x-axis in plots



figure;
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
    
%% 2.1.2
% 160 samples in onze figuur --> tot in het midden ongeveer is dus 80
% samples lang

%% 2.1.3
% more reflections door echo in cathedral
% signaal wordt luider uitgezonden
% signaal wordt beter opgenomen (better recording)
% ---> echte signaal wordt afgezwakt door de afstand DUS detectie is
% moeilijker +++ IR gaat langer zijn door de reflecties

%% 2.1.4
% buiten de omgeving zijn ook de speaker en de ontvanger belangrijk 
% DUS hardware
% OOK stereo/mono

%% 2.1.5
% witte ruis door kanaal geeft de karakteristieken van het kanaal op zich
% --> ontvangen witte ruis convolueren met ideaal, verwacht signaal geeft
% hetzelfde als de kanaalkarakteristieken convolueren met het ideaal,
% verwacht signaal --> DIT is hetzelfde als het signaal door het kanaal sturen


