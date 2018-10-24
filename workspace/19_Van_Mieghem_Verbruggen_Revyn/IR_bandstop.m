addpath(genpath('helper_functions'));
clearvars;
close all;

%%
%%%%%parameters
fs=16000;
h_order = 160;                                      %(bepaald bij IR1.m)  
sig_time = 2;

sig =  bandstop(wgn(fs*3, 1, 0),[700 3000],fs);     %deze doet vanzelf 60dB dus 40dB verzwakking is in orde
                                                    %en gewoon bandstop van
                                                    %white Gaussische ruis
[simin, nbsecs, fs] = initparams(sig, fs);
sim('recplay'); 
rec = simout.signals.values;



start = find(rec > 0.5*max(rec), 1, 'first')-20;    %geeft de eerst index die hieraan voldoet
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

%%%%%plot
t = (0: 1 : (length(h)-1)) /fs;                     % is een rijvector
fft_h = fft(h)';                                    % transpose want we willen een rijvector om te kunnen plotten
length_h = round((length(fft_h)-1)/2 );             % Take half the length
fft_h_trunc = abs(fft_h(1 : length_h));             % Use only first half of fft samples (=positive frequencies).
df_h = (0 : fs /(2*length_h) : (fs-1)/2);           % Frequency vector for x-axis in plots

figure;
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

%% CONCLUSION
% the freq between 700 en 3000 Hz give a slightly positivie freq response,
% while at the other freq, there is a strong negative response
% --> there's no noise transmitted between 700-3000 Hz, so no reflections,
% so in the freq response it seems that 700-3000Hz are really good to send
% data through
% ps: in IR2 the freq response is negative over all freq (because of
% transmitted white noise)

% it does change over the experiments, because enviromental noise changes
% all the time (the only thing that is recorded at those freq is the
% enviromental noise)

    