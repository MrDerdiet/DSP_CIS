
clearvars; hold on; close all;
addpath(genpath('help_functions'), genpath('data'), genpath('mod_func'));

% %% DEMO_1
% % [7.1.1] Construct a QAM symbol sequence of length 1000 = Xk at a single tone k
% % Choose a complex value for the Hk 
% % Generate the received symbol sequence Yk based on (1).
% 
% n = 1000;
% K = 4;
% qam_symbols = qammod((0:2^K-1)',2^K,'bin','UnitAveragePower', true); % generate the different possible qam symboles
% seq = randi([0 1], 1, K*n);
% seq_qam = qam_mod(seq, K);
% 
% H = (0.666 + 0.666j);
% 
% % [7.1.2] Implement the adaptive filter as shown in Fig1 using LMS or NLMS 
% % omit the  noise.   
% % Note  that  the  adaptive  filter in our case has only one (complex) tap.  
% % Choose an initial value for Wk = (1/Hk*)+ delta
% % delta corresponds to a small deviation such that the decision 
% % device still makes correct decisions 
% % (the conjugation of Hk is explained later).
% delta = (1+1j);
% mu = [0.05; 0.5; 5]; %5
% alpha = 10;
% 
% W = zeros(n+1, length(mu));
% W(1,:) = 1/conj(H)+delta;
% W_error = zeros(n, length(mu));
% Xtilde = zeros(n, length(mu));
% Xcirconflex = zeros(n, length(mu));
% legend_names = cell.empty(length(mu), 0);
% 
% Y = H*seq_qam;
% for i=1:length(mu)
%     
%     for k=1:n
%         Xtilde(k,i) = conj(W(k,i))*Y(k);
%         [~, index] = min(abs(qam_symbols-Xtilde(k,i)));
%         Xcirconflex(k,i) = qam_symbols(index);
%         W(k+1,i) =  W(k,i) + mu(i)/(alpha+conj(Y(k))*Y(k))*Y(k)*conj(Xcirconflex(k,i)-Xtilde(k,i));
%         W_error(k,i) = W(k,i) -1/conj(H);
%     end
%     legend_names{i} = strcat('Mu = ',num2str(mu(i)));
% end
% error = Xcirconflex - Xtilde;
% 
% % [7.1.4] Generate a figure that plots the error signal (y-axis)
% % i.e.  the difference between the adaptive filter coefficient
% % Wk* and the inverse channel coefficient 1/Hk,
% % over the iterations (x-axis), for different stepsizes.
% figure('Name', 'W_error signal (constant Hk)');
% subplot(3,1,1);
% plot(real(W_error)); title('Real')
% subplot(3,1,2)
% plot(imag(W_error)); title('imag');
% subplot(3,1,3)
% plot(abs(W_error)); title('abs');
% legend(legend_names);
% 
% figure('Name', 'QAM error signal (constant Hk)');
% subplot(3,1,1);
% plot(real(error)); title('Real')
% subplot(3,1,2)
% plot(imag(error)); title('imag');
% subplot(3,1,3)
% plot(abs(error)); title('abs');
% legend(legend_names);
% 
% % [7.1.5]
% % What happens when the channel changes (slow and fast)?  
% % What does this mean for the OFDM transmission?
% % Only if the adaptive filter converges correctly, you can proceed to
% % the next exercise!
% H = zeros(n+1, 1);
% H(1) = (0.666 + 0.666j);
% for i=1:length(mu)
%     for k=1:n
%         Y(k) = H(k)*seq_qam(k);
%         Xtilde(k,i) = conj(W(k,i))*Y(k);
%         [~, index] = min(abs(qam_symbols-Xtilde(k,i)));
%         Xcirconflex(k,i) = qam_symbols(index);
%         W(k+1,i) =  W(k,i) + mu(i)/(alpha+conj(Y(k))*Y(k))*Y(k)*conj(Xcirconflex(k,i)-Xtilde(k,i));
%         W_error(k,i) = W(k,i) -1/conj(H(k));
% %         H(k+1) = H(k) + (randi([-1000; 1000]) + randi([-1000; 1000])*j)/100000;
%         H(k+1) = H(k) + (randn(1) + randn(1)*1j)/500; 
%     end
%     legend_names{i} = strcat('Mu = ',num2str(mu(i)));
% end
% error = Xcirconflex - Xtilde;
% 
% figure('Name', 'W_error signal (variable Hk)');
% subplot(3,1,1);
% plot(real(W_error)); title('Real')
% subplot(3,1,2)
% plot(imag(W_error)); title('imag');
% subplot(3,1,3)
% plot(abs(W_error)); title('abs');
% legend(legend_names);
% 
% figure('Name', 'QAM Error (variable Hk)');
% subplot(3,1,1);
% plot(real(error)); title('Real')
% subplot(3,1,2)
% plot(imag(error)); title('imag');
% subplot(3,1,3)
% plot(abs(error)); title('abs');
% legend(legend_names);
% 
% figure('Name', 'Hk');
% subplot(3,1,1);
% plot(real(H)); title('Real')
% subplot(3,1,2)
% plot(imag(H)); title('imag');
% subplot(3,1,3)
% plot(abs(H)); title('abs');
% 
% %% DEMO_2
% 

clearvars; 

% Convert BMP image to bitstream
[bitStream_in, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');

for bit_loading = 0: 1
    %% PAUSE
    fprintf('Press any key to continue (in CW)');
    %pause;

    close all;
    clearvars -except bit_loading bitStream_in imageData colorMap imageSize bitsPerPixel; 
        
    %% Parameters
   
    tic;

    fs = 16000;
    N = 512; 
    K = 4; 
    Lt = 3;
    BW = 50;
    threshold = BW/100;
    alpha = 1;
    mu = 1;
    cp_size = 176;
    L = 160; %channel length
    alpha_test = [0.01; 0.1; 1; 10];
    mu_test = [0.01; 0.1; 1; 5];
    
    [pulse, ~] = sinusoid(1000, 1, 1, fs);
    %% channel est for freq_bins
    if bit_loading
        seq = randi([0 1], (N/2-1)*K, 1); % random bits
        trainblock = qam_mod(seq, K); 
        trainblocks = repmat(trainblock,20,1);
        Tx = ofdm_mod_tb(trainblocks, N, cp_size);

        
        [simin,nbsecs,fs,pulse ] = initparams(Tx,fs,L, pulse);

        sim('recplay'); 
        rec = simout.signals.values;

        [Rx_2] = alignIO(rec,pulse,L);
        Rx_2 = Rx_2(1:length(Tx));

        [seq_qam , H] = ofdm_demod_tb(Rx_2, N, cp_size, trainblock);

        receive_deqam = qam_demod(seq_qam, K); 
        trainblocks_deqam = qam_demod(trainblocks, K); 
        

        freq_bins = ofdm_freq_bins(H, N, threshold);
    else
        freq_bins = ones((N/2-1),1);
    end


    %% Generate trainblock, known by both transmitter and receiver
    n_trainblock    = (N/2-1)*K;   % contains N/2-1 random QAM symbols
    trainblock_bit  = randi([0 1], 1, n_trainblock);
    trainblock_qam  = qam_mod(trainblock_bit, K);

    %% Modulate
    qamStream_in = qam_mod(bitStream_in, K);
    [qamStream_in_pad] = pad_input(qamStream_in, length(freq_bins));
    Tx = ofdm_mod_adaptive(qamStream_in_pad, N, cp_size, freq_bins, trainblock_qam, Lt);

    %% Play and record audio
    [simin, nbsecs, fs, pulse] = initparams(Tx,fs,L, pulse);
    sim('recplay')
    out=simout.signals.values;

    %% Align in- and output
    Rx = alignIO(out, pulse,L);
    rest = mod(length(Rx),(N + cp_size));
    Rx = Rx(1: end-rest);

    %% Test for different mu/alpha
    fprintf('\nbegin test \n');
    BER_best = 1;
    i_best =0; j_best = 0;
    for j= 1:length(mu_test)
        for i= 1:length(alpha_test)
            [qamStream_out_pad, ~] = ofdm_demod_adaptive(Rx, N, cp_size, freq_bins, trainblock_qam, Lt, K, alpha_test(i), mu_test(j));
            qamStream_out = qamStream_out_pad(1: length(qamStream_in));
            bitStream_out = qam_demod(qamStream_out, K);

            [ber_IO, ~] = ber(bitStream_in, bitStream_out);
            fprintf(['K=', num2str(K), '\tmu=', num2str(mu_test(j)), '\talpha=', num2str(alpha_test(i)), '\n']);
            fprintf(['\tBER= ', num2str(ber_IO), '\n']);
            if(ber_IO < BER_best)
                fprintf(['test', num2str(i)]);
                i_best = i;
                j_best = j;
                BER_best = ber_IO;
            end
        end
    end
    alpha = alpha_test(i_best);
    mu = mu_test(j_best);
    fprintf('end test \n\n');

    %% Demodulate
    % alpha = 0.01;
    % mu = 1;
    [qamStream_out_pad, W_out] = ofdm_demod_adaptive(Rx, N, cp_size, freq_bins, trainblock_qam, Lt, K, alpha, mu);
    qamStream_out = qamStream_out_pad(1: length(qamStream_in));
    bitStream_out = qam_demod(qamStream_out, K);

    %% BER
    [ber_IO, error_vector] = ber(bitStream_in, bitStream_out);

    fprintf(['mu=', num2str(mu), '\t\talpha=', num2str(alpha), '\n']);
    toc; 
    fprintf(['BER= ', num2str(ber_IO), '\n\n']);

    %% Output image
    % imageRx = bitstreamtoimage(bitStream_out, imageSize, bitsPerPixel);
    % colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;

    %% Output filter/channel response
W = [zeros(1, size(W_out, 2)); W_out; zeros(1,size(W_out, 2)); flipud(conj(W_out))];
w = ifft(W);

H_out = 1./conj(W_out);
H_out_plot = [zeros(1, size(H_out, 2)); H_out];


H = [zeros(1, size(H_out, 2)); H_out; zeros(1, size(H_out, 2)); flipud(conj(H_out))];
h = ifft(W);


%% Visualize
time = 0.00001;


figure('Name','visualization of the demodulation','units','normalized','outerposition',[0 0 1 1])
subplot(2,2,2);colormap(colorMap); image(imageData); axis image; title('The transmitted image'); 
for i = 1:size(h, 2) 
    subplot(2,2,1)
    plot(h(:,i))
    title('Channel in time domain')
%     axis([0 h_order -1 1])
    axis([0 L -0.5 0.5])
    
    subplot(2,2,3)
    plot( (1:N/2) *fs/(N/2-1)/2, 20*log10(abs(H_out_plot(:,i))));
    title('Channel in frequency domain')
    axis([-inf, inf, -40 40]);
    ylabel('Magnitude')
    xlabel('Frequency')
    
    if i*(sum(freq_bins)*K) < sum(bitStream_out)
        imageRx = bitstreamtoimage(bitStream_out(1:i*sum(freq_bins)*K), imageSize, bitsPerPixel);
    else
        imageRx = bitstreamtoimage(bitStream_out, imageSize, bitsPerPixel);
    end
    subplot(2,2,4); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;
%     pause(time)
end    

%% [7.2.4] Slowly increase and decrease the volume during the transmission. 
% Can your (N)LMS filter track the changes in the channel? 
% If not, can you explain why?
%   -> Only in ideal circumstances
%       -> Tracking to slow: cannot track change
%       -> Tracking to fast: filter unstable
%       -> PD controller?

%% [7.2.5] Slowly move the microphone around during the transmission. 
% Can your (N)LMS filter track the changes in the channel? 
% If not, can you explain why?
%   -> see above
end




