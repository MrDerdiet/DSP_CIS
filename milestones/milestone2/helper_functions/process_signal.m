function [freq_spectrum, df] = process_signal(signal, fs, dftsize, time, plot_signal)
% PROCESS_SIGNAL 
%   Calculate FFT of given signal
%   Plot signal in both time and frequency domain if asked
%% INPUT
%   'signal'        
%       -inut data vector 
%   'fs'            
%       -input data samplerate
%   'dftsize'       
%       -fft window size, if >signal_length, the signal is padded with 0
%   'time'          
%       -start time (in s) for fft window
%   'plot'          
%       0 -> no plot
%       1 -> time plot: x-axis = time
%       2 -> time plot: x-axis = sample
%% OUTPUT
%   'freq_spectrum' 
%       -calculated spectrum magnitude vector 
%       -negative frequencies omitted
%   'df'
%       -calculated spectrum frequency vector  

%% PARAMETER CHECK
if nargin < 2
    error('Need at least "signal" and "fs" (samplerate) to process');
end
if nargin < 5
    if nargin < 4
        if (nargin < 3)
            dftsize = length(signal)-1;
        end
        if (dftsize == 0)
            dftsize = length(signal)-1;
        end
        time = 0;
    end
    plot_signal = 0;
end


%% CODE
% Take signal segment at starting time 'time', with length 'dftsize'. 
% Add zeroes if the segment is shorter then 'dftsize'
if (time*fs + dftsize < length(signal)) 
    sig_segment = signal(1+time*fs: time*fs+dftsize );
else
    pad = (time*fs) + dftsize - length(signal);
    sig_segment1 = signal(1+time*fs: length(signal));
    sig_segment = [sig_segment1; zeros(pad, 1)];
end

% fft
fft_sig = fft(sig_segment');
Lsig = round((length(fft_sig)-1)/2 );       % Take half the length
fft_sig_trunc = abs(fft_sig(1 : Lsig));     % Use only first half of fft samples (=positive frequencies).
df_sig = (0 : fs /(2*Lsig) : (fs-1)/2);     % Frequency vector for x-axis in plots

% plot
if (plot_signal)
    figure;
    if (plot_signal == 1)
        t = (0: 1 : (length(signal)-1)) /fs;
        % Time domain plot
        subplot(2, 1, 1); plot(t, signal); title('signal time domain plot / filter impulse response'); xlabel('time (sec)'); ylabel('magnitude');
        % Frequency domain plot
        subplot(2, 1, 2); plot(df_sig, 20*log10(fft_sig_trunc)); title('signal spectrum / filter frequency response'); xlabel('frequency(Hz)'); ylabel('magnitude (dB)');
    else
        t = (0: 1 : (length(signal)-1));
        % Time domain plot
%         subplot(2, 1, 1); plot(t, signal); title('signal time domain plot / filter impulse response'); xlabel('time (sample)'); ylabel('magnitude');
        subplot(2, 1, 1); stem(t, signal, '.'); title('signal time domain plot / filter impulse response'); xlabel('time (sample)'); ylabel('magnitude');
        % Frequency domain plot
        subplot(2, 1, 2); plot(df_sig, 20*log10(fft_sig_trunc)); title('signal spectrum / filter frequency response'); xlabel('frequency(Hz)'); ylabel('magnitude (dB)');
    end
end

% return values
freq_spectrum = fft_sig_trunc;
% freq_spectrum = fft_sig;
df = df_sig;
end
