% OFDM 802.11 PHY 
% End-to-end Simulation
% ECE257A Assignment 1
% Author: Ayush Mukund Jamdar (A69032160)

% TO CHECK: 
% Run this code directly in MATLAB. 
% A summary will be printed in the command window. 
% Some plots might be displayed. Too see all plots, run those particular 
% sections.

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% TRANSMISSION
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% I. Packet Construction and OFDM Modulation
% Step (a) - QPSK Modulation
Nbits = 4160; % should be even
bits = randi([0, 1], Nbits, 1);  % a random vector of size Nx1 

% convert bits to QPSK symbols
qpsk_syms_tx = zeros(Nbits/2, 1);
for i = 1:2:Nbits
    twoBits = string(bits(i)) + string(bits(i+1));
    qi = ceil(i/2);  % qpsk index
    switch twoBits
        case '00'
            qpsk_syms_tx(qi) = 1 + 0j;
        case '01'
            qpsk_syms_tx(qi) = 0 + 1j;
        case '11'
            qpsk_syms_tx(qi) = 0 - 1j;
        case '10'
            qpsk_syms_tx(qi) = -1 + 0j;
    end
end

% Group QPSK Symbols into OFDM Symbols
% Map each set of 48 data symbols into a block of 64 samples
% Represent as a matrix of Num_of_ofdm_syms x 64
N_ofdm_syms = ceil(size(qpsk_syms_tx, 1) / 48);  % groups of 48
ofdm_syms_tx = zeros(N_ofdm_syms, 64);  % nulls are assigned 0 by default

% pilots at indices = 8, 22, 44, 58
pilot_indices = [8 22 44 58];
ofdm_syms_tx(:, pilot_indices) = 1+0j;

% group qpsk into ofdm
ofdm_sym_index = 1;
data_indices = cat(2, 2:7, 9:21, 23:27, 39:43, 45:57, 59:64);

for i = 1:48:size(qpsk_syms_tx, 1)
    if i+47 > size(qpsk_syms_tx, 1) 
        % not always will there be groups of 48 exactly
        % so we append zeros to the last incomplete group
        last_qpskSyms = qpsk_syms_tx(i:end);
        ofdm_syms_tx(ofdm_sym_index, data_indices) = cat(1, last_qpskSyms, ...
            zeros(48-size(last_qpskSyms, 1), 1));
    else
        ofdm_syms_tx(ofdm_sym_index, data_indices) = qpsk_syms_tx(i:i+47);
    end
    ofdm_sym_index = ofdm_sym_index + 1;
end

%% 
% Step (b) FFT and Cyclic Prefix
xt_tx_nocp = ifft(ofdm_syms_tx, 64, 2); % fft of each row
cyclic_prefix_xt = xt_tx_nocp(:, end-15:end);
xt_tx = cat(2, cyclic_prefix_xt, xt_tx_nocp); % a matrix of N_ofdm_symsx80
% this is only the data - now add STF LTF 

%%
% Step (c) STF and LTF
stf_sym = zeros(1, 64);
% the stf/ltf sequence given in 802.11 docs is fftshifted
% hence, move the two parts around
stf_sym(1, 39:64) = sqrt(13/6) * [0 0 1+1j 0 0 0 -1-1j 0 0 0 1+1j ...
    0 0 0 -1-1j 0 0 0 -1-1j 0 0 0 1+1j 0 0 0];
stf_sym(1, 1:27) = sqrt(13/6) * [0 0 0 0 -1-1j 0 0 0 -1-1j 0 ...
    0 0 1+1j 0 0 0 1+1j 0 0 0 1+1j 0 0 0 1+1j 0 0];

ltf_sym = zeros(1, 64);
ltf_sym(1, 39:64) = [1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 ...
    -1 1 -1 1 1 1 1];
ltf_sym(1, 1:27) = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 ...
    1 -1 1 1 1 1];

% the above symbols are in the frequency domain
% convert them to time domain before appending to the packet
% add cyclic prefix 64+16 = 80 samples per OFDM symbol
stf_t = ifft(stf_sym, 64); stf_t = cat(2, stf_t(end-15:end), stf_t);
ltf_t = ifft(ltf_sym, 64); ltf_t = cat(2, ltf_t(end-15:end), ltf_t);

% we need two such symbols
stf_t = cat(2, stf_t, stf_t); ltf_t = cat(2, ltf_t, ltf_t);

% plot STF to visualise
plot(cat(2, real(stf_t), real(ltf_t)), LineWidth=1.5);
xline(160, 'red', LineWidth=2);
grid on;
xlabel('Sample Index');
ylabel('Magnitude');
title('STF ------ LTF'); 

packet_tx = cat(2, stf_t, ltf_t, reshape(xt_tx.', 1, [])); 
% flatten the ofdm matrix and attach it to training fields

%%
% Step (d)
% Power Spectrum Density of OFDM Symbols
% 64 pt fft; 20MHz sampling rate
[psd, f] = pwelch(packet_tx, [], [], 64, 20e6);
psd = fftshift(psd); 
f = cat(1, -flip(f(1:32)), f(1:32));
plot(f/1e6, 10*log10(psd));
xlabel('Frequency (MHz)');
ylabel('PSD Magnitude (dB)');
title('Power Spectral Density of OFDM Symbols');
grid on;

%%
fprintf('Transmission data:\n');
fprintf('# Bits sent: %d\n', Nbits);
fprintf('# QPSK symbols sent: %d\n', length(qpsk_syms_tx));
fprintf('# OFDM symbols in packet: %d\n', size(ofdm_syms_tx, 1));
fprintf('Packet sent successfully!\n');

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% CHANNEL MODELLING
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% II. Packet Transmission and Channel Distortion
idle_period = zeros(1, 100);
packet_idletime = cat(2, idle_period, packet_tx);

attenuation = 1e-5; 
phase_distortion = exp(-pi*3j/4);
% channel distortion
packet_rx = packet_idletime * attenuation * phase_distortion; 

% Carrier frequency offset
% Note: this is actually fc.ts = phase drift
freq_offset = 0.00017;
k = 0:length(packet_tx)-1;
cfo_array = cat(2, ones(1, 100), exp(-1j*2*pi*freq_offset*k));
yt_rx_prenoise = packet_rx .* cfo_array; 

% add zero mean gaussian noise
variance = 1e-14;
yt_rx = yt_rx_prenoise + randn(1, length(yt_rx_prenoise))*sqrt(variance);

% plot the packet's initial samples
plot(abs(yt_rx(1:400)));
grid on;
xlabel('Sample Index');
ylabel('Magnitude');
title('Received Signal: Idle-STF-LTF');


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% RECEPTION
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% III. Packet Detection
% Self Correlation
r_m = zeros(1, length(yt_rx)-16);
e_m = zeros(1, length(yt_rx)-16);
for i = 17:length(yt_rx)-16
    r_m(1, i-16) = dot(yt_rx(i:i+15), yt_rx(i-16:i-16+15));
    e_m(1, i-16) = dot(yt_rx(i:i+15), yt_rx(i:i+15));
end

% Plot
plot(abs(r_m)); hold on;
plot(e_m); grid on;
legend('R(m)', 'E(m)');
hold off;

%% IV. Packet Synchronization
% Cross Correlation
stf_16 = stf_t(1, 1:16);
r_m_cross = zeros(1, length(yt_rx)-16);
detected_stf_ind = zeros(1, 10);  % to store STF sample indices
stf_count = 0; 

for i = 1:length(yt_rx)-16
    r_m_cross(1, i) = dot(yt_rx(i:i+15), stf_16);
    if (abs(r_m_cross(1, i)) > 1.5e-6) && stf_count < 10
        stf_count = stf_count + 1;
        detected_stf_ind(1, stf_count) = i;
    end
end

% we now know STF start indices (#10) in the received yt

% Plot
plot(abs(r_m_cross(1:300))); grid on;
title('Cross Correlation')
hold off;

%% V. Estimation
% (a) CFO
% Calculate frequency offset using LTF
% We have 2.5 repeated LTF symbols
ltf_start_ind = detected_stf_ind(end) + 16;

% phase drift = fc.Ts
% we find an average using the first LTF symbol
phase_drift = 0; % init
for k=0:63
    phase_drift = phase_drift + angle(yt_rx(ltf_start_ind+k)/yt_rx(ltf_start_ind+k+64))/(64*2*pi);
end
phase_drift = phase_drift / 64; % averaging
% ts = find out
% cfo = phase_drift / ts;
fprintf('\nReception data: \n');
fprintf('Detected STF starting time sample: %d\n', detected_stf_ind(1));
disp('Actual Phase Drift: ' + string(freq_offset));
disp('Estimated Phase Drift: ' + string(phase_drift));

%%
% (b) Channel Estimation
% first LTF symbol was used to phase drift
% we'll use the second LTF symbol to estimate channel
% this includes frequency offset correction
% LTF is known to the receiver
cfo_correctn_ltf = 240+16:319;
h_channel = yt_rx(ltf_start_ind+96:ltf_start_ind+159) .* ...
    (exp(2j*pi*phase_drift*cfo_correctn_ltf)) ./ ltf_t(97:160);

subplot(2, 1, 1);
plot(abs(h_channel)); grid on;
xlabel('OFDM Carrier');
ylabel('Magnitude');
title('Magnitude Distortion');
subplot(2, 1, 2);
plot(angle(h_channel)); grid on;
xlabel('OFDM Carrier');
ylabel('Phase');
title('Phase Distortion');
sgtitle('Channel Distortion to Subcarriers');


%% VI. Decoding
ofdm_start_ind = ltf_start_ind + 160;
qpsk_syms_rx = [];
ofdm_sym_num_rx = 0;

% group 64 samples to run fft
for i = ofdm_start_ind+16:80:length(yt_rx)
    % for each ofdm symbol 
    % equalization
    cfo_correctn_rx = 320+16+80*(ofdm_sym_num_rx):320+80*(ofdm_sym_num_rx+1)-1;
    ofdm_sym_num_rx = ofdm_sym_num_rx + 1;
    xt_rx_ofdm = yt_rx(i:i+63) .* exp(2j*pi*phase_drift*cfo_correctn_rx) ./ h_channel;
    % removed cyclic prefix; corrected frequency offset

    % fft: one ofdm symbol
    ofdm_sym_rx = fft(xt_rx_ofdm);
   
    % decode: map to QPSK
    % this contains nulls and pilots which will get assigned to nearest 
    % qpsk symbol - there will be no nulls
    decoded_sym = []; % 64
    for j=1:length(ofdm_sym_rx)
        % nearest neighbour decoding logic
        realpart = real(ofdm_sym_rx(j));
        imagpart = imag(ofdm_sym_rx(j));
        if (imagpart<realpart) && (imagpart>=-realpart)
            decoded_sym = [decoded_sym 1+0j];
        elseif (imagpart>=realpart) && (imagpart>-realpart)
            decoded_sym = [decoded_sym 0+1j];
        elseif (imagpart>realpart) && (imagpart<=-realpart)
            decoded_sym = [decoded_sym -1+0j];
        elseif (imagpart<=realpart) && (imagpart<-realpart)
            decoded_sym = [decoded_sym 0-1j];
        else
            decoded_sym = [decoded_sym 0+0j];
        end
    end

    % this ofdm symbol contains nulls and pilots
    % extract the data points
    qpsk_syms_rx = [qpsk_syms_rx  decoded_sym(data_indices)];
end 

% we had appended zeros when the last set of 48 data points for one OFDM 
% symbol could not be completed; remove those data from received symbols
qpsk_syms_rx = qpsk_syms_rx(1:length(qpsk_syms_tx));
fprintf('# QPSK symbols received: %d\n', length(qpsk_syms_rx));
fprintf('Packet received successfully!\n');


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% VI. Evaluation
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% find error rate
error_rate = sum(qpsk_syms_rx ~= qpsk_syms_tx.') / length(qpsk_syms_rx);
fprintf('Error Rate: %.2f%%\n', error_rate * 100);
