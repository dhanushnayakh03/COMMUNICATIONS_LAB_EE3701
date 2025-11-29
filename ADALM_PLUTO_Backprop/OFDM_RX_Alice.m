% close all; clear; clc; j = 1i;
% Global_Parameters_PLS;
% 
% %% Dual Source Receiver for Simultaneous Bob and Eve Detection
% % This receiver can simultaneously synchronize with both transmitters
% % using their orthogonal preambles
% 
% %% Load or Generate Received Signal
% % In practice, this would come from Hardware_RX or simulation
% % For testing, we'll create a synthetic received signal
% 
% % Load both TX signals
% load('TX_signal_Bob.mat');
% load('TX_signal_Eve.mat');
% 
% % Create combined received signal with different channels and timing offsets
% SNR_dB = 20;
% delay_Bob = 50; % samples
% delay_Eve = 75; % samples
% 
% % Channel coefficients
% h_Bob = [0.8, 0.3*exp(1i*pi/4), 0.1*exp(1i*pi/3)];
% h_Eve = [0.7, 0.4*exp(1i*pi/6), 0.15*exp(1i*pi/2)];
% 
% % Apply channels
% RX_Bob = conv(TX_signal_Bob, h_Bob);
% RX_Eve = conv(TX_signal_Eve, h_Eve);
% 
% % Add delays and combine
% max_len = max(length(RX_Bob) + delay_Bob, length(RX_Eve) + delay_Eve);
% RX_combined = zeros(1, max_len);
% RX_combined(delay_Bob + (1:length(RX_Bob))) = RX_combined(delay_Bob + (1:length(RX_Bob))) + RX_Bob;
% RX_combined(delay_Eve + (1:length(RX_Eve))) = RX_combined(delay_Eve + (1:length(RX_Eve))) + RX_Eve;
% 
% % Add noise
% signal_power = mean(abs(RX_combined).^2);
% noise_power = signal_power / (10^(SNR_dB/10));
% noise = sqrt(noise_power/2) * (randn(size(RX_combined)) + 1i*randn(size(RX_combined)));
% RX_signal = RX_combined + noise;
% 
% %% Matched Filtering with Orthogonal Preambles
% 
% % Bob's matched filter (EVEN subcarriers)
% Bob_short_preamble_time = ifft(ifftshift(Parameters_struct.Short_preamble_Bob_Frequency));
% Bob_short_preamble = repmat(Bob_short_preamble_time(1:16), 1, 10);
% Bob_matched_filter = conj(fliplr(Bob_short_preamble));
% 
% % Eve's matched filter (ODD subcarriers)
% Eve_short_preamble_time = ifft(ifftshift(Parameters_struct.Short_preamble_Eve_Frequency));
% Eve_short_preamble = repmat(Eve_short_preamble_time(1:16), 1, 10);
% Eve_matched_filter = conj(fliplr(Eve_short_preamble));
% 
% % Cross-correlations
% corr_Bob = abs(conv(RX_signal, Bob_matched_filter));
% corr_Eve = abs(conv(RX_signal, Eve_matched_filter));
% 
% % Normalize correlations
% corr_Bob = corr_Bob / max(corr_Bob);
% corr_Eve = corr_Eve / max(corr_Eve);
% 
% %% Detect Timing Offsets
% threshold = 0.7;
% 
% % Find peaks for Bob
% [pks_Bob, locs_Bob] = findpeaks(corr_Bob, 'MinPeakHeight', threshold, 'MinPeakDistance', 100);
% if ~isempty(locs_Bob)
%     timing_offset_Bob = locs_Bob(1) - length(Bob_matched_filter);
%     detected_Bob = true;
% else
%     timing_offset_Bob = 0;
%     detected_Bob = false;
% end
% 
% % Find peaks for Eve
% [pks_Eve, locs_Eve] = findpeaks(corr_Eve, 'MinPeakHeight', threshold, 'MinPeakDistance', 100);
% if ~isempty(locs_Eve)
%     timing_offset_Eve = locs_Eve(1) - length(Eve_matched_filter);
%     detected_Eve = true;
% else
%     timing_offset_Eve = 0;
%     detected_Eve = false;
% end
% 
% %% Fine Timing and Channel Estimation using Long Preambles
% 
% if detected_Bob
%     % Extract Bob's long preamble
%     Bob_long_start = timing_offset_Bob + length(Bob_short_preamble);
%     Bob_long_preamble_time = ifft(ifftshift(Parameters_struct.Long_preamble_Bob_Frequency));
% 
%     if Bob_long_start + 160 <= length(RX_signal)
%         RX_Bob_long = RX_signal(Bob_long_start + (1:160));
% 
%         % Channel estimation for Bob (using two repetitions)
%         RX_Bob_long1 = RX_Bob_long(33:96);
%         RX_Bob_long2 = RX_Bob_long(97:160);
% 
%         H_Bob_time = (RX_Bob_long1 + RX_Bob_long2) / (2 * Bob_long_preamble_time);
%         H_Bob_freq = fft(fftshift(H_Bob_time));
%     else
%         H_Bob_freq = ones(1, Parameters_struct.N_FFT);
%         disp('Warning: Bob long preamble extraction failed');
%     end
% else
%     H_Bob_freq = ones(1, Parameters_struct.N_FFT);
% end
% 
% if detected_Eve
%     % Extract Eve's long preamble
%     Eve_long_start = timing_offset_Eve + length(Eve_short_preamble);
%     Eve_long_preamble_time = ifft(ifftshift(Parameters_struct.Long_preamble_Eve_Frequency));
% 
%     if Eve_long_start + 160 <= length(RX_signal)
%         RX_Eve_long = RX_signal(Eve_long_start + (1:160));
% 
%         % Channel estimation for Eve (using two repetitions)
%         RX_Eve_long1 = RX_Eve_long(33:96);
%         RX_Eve_long2 = RX_Eve_long(97:160);
% 
%         H_Eve_time = (RX_Eve_long1 + RX_Eve_long2) / (2 * Eve_long_preamble_time);
%         H_Eve_freq = fft(fftshift(H_Eve_time));
%     else
%         H_Eve_freq = ones(1, Parameters_struct.N_FFT);
%         disp('Warning: Eve long preamble extraction failed');
%     end
% else
%     H_Eve_freq = ones(1, Parameters_struct.N_FFT);
% end
% 
% %% Payload Extraction and Decoding
% 
% N_FFT = Parameters_struct.N_FFT;
% Bob_pilot_idx = Parameters_struct.Bob_pilot_idx;
% Eve_pilot_idx = Parameters_struct.Eve_pilot_idx;
% data_idx = Parameters_struct.data_idx;
% 
% % Bob's Payload Extraction
% if detected_Bob
%     Bob_payload_start = timing_offset_Bob + 320; % After short + long preambles
% 
%     % Symbol 1
%     if Bob_payload_start + 80 <= length(RX_signal)
%         Bob_sym1_rx = RX_signal(Bob_payload_start + (1:80));
%         Bob_sym1_no_cp = Bob_sym1_rx(17:80); % Remove CP
%         Bob_sym1_freq = fftshift(fft(Bob_sym1_no_cp));
% 
%         % Equalization using pilot-based channel estimate
%         H_Bob_pilots = H_Bob_freq(Bob_pilot_idx);
%         Bob_sym1_eq = Bob_sym1_freq ./ H_Bob_freq;
% 
%         % Extract data
%         Bob_data1_eq = Bob_sym1_eq(data_idx);
%         Bob_data1_decoded = pskdemod(Bob_data1_eq, 4, pi/4);
%     else
%         Bob_data1_decoded = [];
%     end
% 
%     % Symbol 2
%     if Bob_payload_start + 160 <= length(RX_signal)
%         Bob_sym2_rx = RX_signal(Bob_payload_start + 80 + (1:80));
%         Bob_sym2_no_cp = Bob_sym2_rx(17:80);
%         Bob_sym2_freq = fftshift(fft(Bob_sym2_no_cp));
%         Bob_sym2_eq = Bob_sym2_freq ./ H_Bob_freq;
%         Bob_data2_eq = Bob_sym2_eq(data_idx);
%         Bob_data2_decoded = pskdemod(Bob_data2_eq, 4, pi/4);
%     else
%         Bob_data2_decoded = [];
%     end
% else
%     Bob_data1_decoded = [];
%     Bob_data2_decoded = [];
% end
% 
% % Eve's Payload Extraction
% if detected_Eve
%     Eve_payload_start = timing_offset_Eve + 320;
% 
%     % Symbol 1
%     if Eve_payload_start + 80 <= length(RX_signal)
%         Eve_sym1_rx = RX_signal(Eve_payload_start + (1:80));
%         Eve_sym1_no_cp = Eve_sym1_rx(17:80);
%         Eve_sym1_freq = fftshift(fft(Eve_sym1_no_cp));
%         Eve_sym1_eq = Eve_sym1_freq ./ H_Eve_freq;
%         Eve_data1_eq = Eve_sym1_eq(data_idx);
%         Eve_data1_decoded = pskdemod(Eve_data1_eq, 4, pi/4);
%     else
%         Eve_data1_decoded = [];
%     end
% 
%     % Symbol 2
%     if Eve_payload_start + 160 <= length(RX_signal)
%         Eve_sym2_rx = RX_signal(Eve_payload_start + 80 + (1:80));
%         Eve_sym2_no_cp = Eve_sym2_rx(17:80);
%         Eve_sym2_freq = fftshift(fft(Eve_sym2_no_cp));
%         Eve_sym2_eq = Eve_sym2_freq ./ H_Eve_freq;
%         Eve_data2_eq = Eve_sym2_eq(data_idx);
%         Eve_data2_decoded = pskdemod(Eve_data2_eq, 4, pi/4);
%     else
%         Eve_data2_decoded = [];
%     end
% else
%     Eve_data1_decoded = [];
%     Eve_data2_decoded = [];
% end
% 
% %% Performance Analysis
% 
% % Calculate BER for Bob
% if ~isempty(Bob_data1_decoded)
%     Bob_errors1 = sum(Bob_data1_decoded ~= Parameters_struct.data_Payload_Bob_1.data_Payload_1);
%     Bob_BER1 = Bob_errors1 / length(Bob_data1_decoded);
% else
%     Bob_BER1 = NaN;
% end
% 
% if ~isempty(Bob_data2_decoded)
%     Bob_errors2 = sum(Bob_data2_decoded ~= Parameters_struct.data_Payload_Bob_2.data_Payload_2);
%     Bob_BER2 = Bob_errors2 / length(Bob_data2_decoded);
% else
%     Bob_BER2 = NaN;
% end
% 
% % Calculate BER for Eve
% if ~isempty(Eve_data1_decoded)
%     Eve_errors1 = sum(Eve_data1_decoded ~= Parameters_struct.data_Payload_Eve_1.data_Payload_1);
%     Eve_BER1 = Eve_errors1 / length(Eve_data1_decoded);
% else
%     Eve_BER1 = NaN;
% end
% 
% if ~isempty(Eve_data2_decoded)
%     Eve_errors2 = sum(Eve_data2_decoded ~= Parameters_struct.data_Payload_Eve_2.data_Payload_2);
%     Eve_BER2 = Eve_errors2 / length(Eve_data2_decoded);
% else
%     Eve_BER2 = NaN;
% end
% 
% %% Visualization
% 
% figure('Name', 'Dual Source Detection and Synchronization', 'Position', [50 50 1600 900]);
% 
% % Correlation plots
% subplot(3,3,1);
% plot(corr_Bob, 'b', 'LineWidth', 1.5);
% hold on;
% plot(corr_Eve, 'r', 'LineWidth', 1.5);
% if detected_Bob
%     xline(locs_Bob(1), 'b--', 'LineWidth', 2);
% end
% if detected_Eve
%     xline(locs_Eve(1), 'r--', 'LineWidth', 2);
% end
% yline(threshold, 'k--');
% title('Orthogonal Preamble Correlations');
% xlabel('Sample'); ylabel('Normalized Correlation');
% legend('Bob (EVEN)', 'Eve (ODD)', 'Location', 'best');
% grid on;
% 
% % Channel estimates
% subplot(3,3,2);
% plot(abs(H_Bob_freq), 'b', 'LineWidth', 1.5);
% hold on;
% plot(abs(H_Eve_freq), 'r', 'LineWidth', 1.5);
% title('Estimated Channel Frequency Response');
% xlabel('Subcarrier'); ylabel('Magnitude');
% legend('Bob Channel', 'Eve Channel');
% grid on;
% 
% subplot(3,3,3);
% plot(angle(H_Bob_freq), 'b', 'LineWidth', 1.5);
% hold on;
% plot(angle(H_Eve_freq), 'r', 'LineWidth', 1.5);
% title('Channel Phase Response');
% xlabel('Subcarrier'); ylabel('Phase (rad)');
% legend('Bob Channel', 'Eve Channel');
% grid on;
% 
% % Constellation diagrams
% if ~isempty(Bob_data1_eq)
%     subplot(3,3,4);
%     plot(Bob_data1_eq, 'b.', 'MarkerSize', 10);
%     title('Bob Symbol 1 - Equalized');
%     xlabel('Real'); ylabel('Imaginary');
%     axis equal; grid on;
% end
% 
% if ~isempty(Eve_data1_eq)
%     subplot(3,3,5);
%     plot(Eve_data1_eq, 'r.', 'MarkerSize', 10);
%     title('Eve Symbol 1 - Equalized');
%     xlabel('Real'); ylabel('Imaginary');
%     axis equal; grid on;
% end
% 
% % Preamble spectra
% subplot(3,3,6);
% stem(abs(fftshift(Parameters_struct.Short_preamble_Bob_Frequency)), 'b');
% hold on;
% stem(abs(fftshift(Parameters_struct.Short_preamble_Eve_Frequency)), 'r');
% title('Orthogonal Short Preambles - Frequency Domain');
% xlabel('Subcarrier'); ylabel('Magnitude');
% legend('Bob (EVEN)', 'Eve (ODD)');
% grid on;
% 
% % Timing detection details
% subplot(3,3,7);
% text(0.1, 0.9, 'Detection Results:', 'FontSize', 12, 'FontWeight', 'bold');
% text(0.1, 0.75, ['Bob Detected: ', num2str(detected_Bob)], 'FontSize', 10, 'Color', 'blue');
% text(0.1, 0.65, ['Bob Timing Offset: ', num2str(timing_offset_Bob)], 'FontSize', 10, 'Color', 'blue');
% text(0.1, 0.55, ['Bob Peak: ', num2str(max(corr_Bob), '%.3f')], 'FontSize', 10, 'Color', 'blue');
% text(0.1, 0.40, ['Eve Detected: ', num2str(detected_Eve)], 'FontSize', 10, 'Color', 'red');
% text(0.1, 0.30, ['Eve Timing Offset: ', num2str(timing_offset_Eve)], 'FontSize', 10, 'Color', 'red');
% text(0.1, 0.20, ['Eve Peak: ', num2str(max(corr_Eve), '%.3f')], 'FontSize', 10, 'Color', 'red');
% axis off;
% 
% % BER Results
% subplot(3,3,8);
% text(0.1, 0.9, 'BER Performance:', 'FontSize', 12, 'FontWeight', 'bold');
% text(0.1, 0.75, 'Bob:', 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'blue');
% text(0.1, 0.65, ['  Symbol 1 BER: ', num2str(Bob_BER1, '%.4f')], 'FontSize', 10, 'Color', 'blue');
% text(0.1, 0.55, ['  Symbol 2 BER: ', num2str(Bob_BER2, '%.4f')], 'FontSize', 10, 'Color', 'blue');
% text(0.1, 0.40, 'Eve:', 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'red');
% text(0.1, 0.30, ['  Symbol 1 BER: ', num2str(Eve_BER1, '%.4f')], 'FontSize', 10, 'Color', 'red');
% text(0.1, 0.20, ['  Symbol 2 BER: ', num2str(Eve_BER2, '%.4f')], 'FontSize', 10, 'Color', 'red');
% axis off;
% 
% % System info
% subplot(3,3,9);
% text(0.1, 0.9, 'System Configuration:', 'FontSize', 12, 'FontWeight', 'bold');
% text(0.1, 0.75, ['SNR: ', num2str(SNR_dB), ' dB'], 'FontSize', 10);
% text(0.1, 0.65, ['FFT Size: ', num2str(N_FFT)], 'FontSize', 10);
% text(0.1, 0.55, ['Bob Pilots: ', num2str(Bob_pilot_idx)], 'FontSize', 9);
% text(0.1, 0.45, ['Eve Pilots: ', num2str(Eve_pilot_idx)], 'FontSize', 9);
% text(0.1, 0.30, 'Orthogonality:', 'FontSize', 11, 'FontWeight', 'bold');
% text(0.1, 0.20, '• Short: Frequency Division', 'FontSize', 9);
% text(0.1, 0.10, '• Long: Code Division', 'FontSize', 9);
% axis off;
% 
% %% Display Results
% disp('========================================');
% disp('DUAL SOURCE DETECTION RESULTS');
% disp('========================================');
% disp(['SNR: ', num2str(SNR_dB), ' dB']);
% disp(' ');
% disp('Bob (Legitimate TX):');
% disp(['  Detected: ', num2str(detected_Bob)]);
% if detected_Bob
%     disp(['  Timing Offset: ', num2str(timing_offset_Bob), ' samples']);
%     disp(['  Correlation Peak: ', num2str(max(corr_Bob), '%.4f')]);
%     disp(['  Symbol 1 BER: ', num2str(Bob_BER1, '%.4e')]);
%     disp(['  Symbol 2 BER: ', num2str(Bob_BER2, '%.4e')]);
% end
% disp(' ');
% disp('Eve (Eavesdropper/Helper TX):');
% disp(['  Detected: ', num2str(detected_Eve)]);
% if detected_Eve
%     disp(['  Timing Offset: ', num2str(timing_offset_Eve), ' samples']);
%     disp(['  Correlation Peak: ', num2str(max(corr_Eve), '%.4f')]);
%     disp(['  Symbol 1 BER: ', num2str(Eve_BER1, '%.4e')]);
%     disp(['  Symbol 2 BER: ', num2str(Eve_BER2, '%.4e')]);
% end
% disp('========================================');




function [M_n, Threshold_graph, H_Bob, H_Eve, channel_correlation, ...
          RX_Bob_Payload_1, RX_Bob_Payload_2, RX_Eve_Payload_1, RX_Eve_Payload_2, ...
          BER_Bob, BER_Eve, H_direct] = OFDM_RX_Alice(RX, Parameters_struct)
%% OFDM_RX_Alice - Receiver for Alice in Physical Layer Security System
% 
% This function processes the superimposed signals from Bob and Eve,
% separates them using orthogonal pilot patterns, and estimates both channels.
%
% Channel Tap Logging Feature:
%   - Logs Bob and Eve channel taps (time domain) when BER ≈ 0
%   - Only logs first 250 frames
%   - Separate logs: 'Bob_channel_taps_log.mat' and 'Eve_channel_taps_log.mat'
%
% INPUTS:
%   RX                - Received signal [1 x N_samples]
%   Parameters_struct - System parameters structure
%
% OUTPUTS:
%   M_n                  - Packet detection metric
%   Threshold_graph      - Detection threshold visualization
%   H_Bob                - Bob's channel estimate (time domain, padded to 64)
%   H_Eve                - Eve's channel estimate (time domain, padded to 64)
%   channel_correlation  - Correlation coefficient between channels
%   RX_Bob_Payload_1/2   - Bob's equalized constellation points (44 each)
%   RX_Eve_Payload_1/2   - Eve's equalized constellation points (44 each)
%   BER_Bob              - Bit Error Rate for Bob's channel
%   BER_Eve              - Bit Error Rate for Eve's channel
%   H_direct             - Direct channel estimate at pilot locations

%% ========== INITIALIZATION ==========
Debug_mode = 'off';
if strcmp(Debug_mode, 'on')
    clearvars -except Debug_mode;
    close all; clc;
    Global_Parameters_PLS;
    load('RX_test_signal.mat'); % Load test RX signal
end

j = 1i;
N_FFT = Parameters_struct.N_FFT;

% Determine payload size from Parameters_struct or default to 44
if isfield(Parameters_struct, 'N_data')
    N_data = Parameters_struct.N_data;
else
    N_data = 44; % Default payload size
end

%% ========== CHANNEL TAP LOGGING SETUP ==========
persistent frame_counter Bob_tap_log Eve_tap_log;

% Initialize on first call
if isempty(frame_counter)
    frame_counter = 0 ; 
Bob_tap_log = struct([]);   % Always reset
Eve_tap_log = struct([]);

    Bob_tap_log = struct('frame_num', {}, 'taps', {}, 'BER', {}, 'timestamp', {});
    Eve_tap_log = struct('frame_num', {}, 'taps', {}, 'BER', {}, 'timestamp', {});
end

% Increment frame counter
frame_counter = frame_counter + 1;

% BER threshold for "close to zero"
BER_THRESHOLD = 1e-4; % Adjust as needed (0.001 = 0.1%)
MAX_FRAMES_TO_LOG = 50;

%% ========== PULSE SHAPING FILTER ==========
% Root Raised Cosine matched filter
rolloff = 0.5;
L_RRC = 6;
OVR = 2;
RRC = rcosdesign(rolloff, L_RRC, OVR, 'sqrt'); % [1x13]
RX_signal = conv(RX, RRC); % Matched filtering [1x3012]

%% ========== PACKET DETECTION (SCHMIDL-COX ALGORITHM) ==========
% Uses repetitive structure of short preamble for synchronization
D = 16; % Delay (half of short preamble slot)
L = 32; % Window length for averaging

C_n = zeros(1,length(RX)-D+1-L);
P_n = zeros(1,length(RX)-D+1-L);
C_k = zeros(1,L);
P_k = zeros(1,L);

for n=1:length(RX)-D+1-L
    for k=1:L
        C_k(k) = RX(n+k-1)*complex(RX(n+k-1+D));
        P_k(k) = abs(RX(n+k-1+D))^2;
    end
    C_n(n) = sum(C_k);
    P_n(n) = sum(P_k);
end
M_n = (abs(C_n).^2)./(P_n.^2);

%% ========== PACKET START DETECTION ==========
Threshold = 0.7;
loc = find(M_n > Threshold);

if isempty(loc)
    warning('No packet detected! Check signal strength or threshold.');
    % Return default values
    M_n = zeros(1, 100);
    Threshold_graph = Threshold * ones(1, 100);
    H_Bob = zeros(1, N_FFT);
    H_Eve = zeros(1, N_FFT);
    channel_correlation = 0;
    RX_Bob_Payload_1 = zeros(1, N_data);
    RX_Bob_Payload_2 = zeros(1, N_data);
    RX_Eve_Payload_1 = zeros(1, N_data);
    RX_Eve_Payload_2 = zeros(1, N_data);
    BER_Bob = 0.5;
    BER_Eve = 0.5;
    H_direct = zeros(N_FFT,1);
    return;
end

% Find packet front (large gap in detection metric)
temp_1 = [loc, 0];
temp_2 = [0, loc];
temp_3 = temp_1 - temp_2;
Packet_Front = find(temp_3 > 300);

if isempty(Packet_Front)
    Packet_Front_idx = loc(1);
else
    Packet_Front_idx = loc(Packet_Front);
end

% Verify packet is sustained above threshold
Length_over_Threshold = 230;
idx = Packet_Front_idx(1) + L_RRC + 1; % Default

for x = 1:length(Packet_Front_idx)
    if Packet_Front_idx(x) + Length_over_Threshold <= length(M_n)
        if M_n(Packet_Front_idx(x) + Length_over_Threshold) > Threshold
            idx = Packet_Front_idx(x) + L_RRC + 1;
            break;
        end
    end
end

% Visualization marker
Threshold_graph = Threshold * ones(1, length(M_n));
if idx - L_RRC - 1 <= length(Threshold_graph)
    Threshold_graph(idx - L_RRC - 1) = 1.15;
end

%% ========== FRAME EXTRACTION AND DOWNSAMPLING ==========
% Frame structure: Short(160) + Long(160) + Payload1(80) + Payload2(80) = 480 samples
Frame_length = 480;

% Extract frame and downsample by OVR
if idx + OVR*Frame_length - 1 > length(RX_signal)
    warning('Packet position exceeds signal length. Adjusting...');
    idx = max(1, length(RX_signal) - OVR*Frame_length);
end

Frame_DWN_sampling = RX_signal(idx:OVR:min(idx + OVR*Frame_length - 1, length(RX_signal)));

% Ensure frame length is correct
if length(Frame_DWN_sampling) < Frame_length
    warning('Frame too short after downsampling. Padding with zeros.');
    Frame_DWN_sampling = [Frame_DWN_sampling, zeros(1, Frame_length - length(Frame_DWN_sampling))];
end
Frame_DWN_sampling = Frame_DWN_sampling(1:Frame_length); % Trim to exact length

%% ========== COARSE CFO ESTIMATION ==========
% Uses correlation between repeated short preamble slots
Short_preamble_slot_length = 16;

slot5 = Frame_DWN_sampling(Short_preamble_slot_length*5 + 1 : Short_preamble_slot_length*6);
slot6 = Frame_DWN_sampling(Short_preamble_slot_length*6 + 1 : Short_preamble_slot_length*7);
z_coarse = slot5 * slot6';

f_Coarse_est = (-1/(2*pi*Short_preamble_slot_length*Parameters_struct.Ts)) * angle(z_coarse);
Frame_After_Coarse = Frame_DWN_sampling .* exp(-j*2*pi*f_Coarse_est*Parameters_struct.Ts*(0:Frame_length-1));

%% ========== FINE CFO ESTIMATION ==========
% Uses correlation between two long preamble symbols
long_preamble_1 = Frame_After_Coarse(Short_preamble_slot_length*12 + 1 : Short_preamble_slot_length*16);
long_preamble_2 = Frame_After_Coarse(Short_preamble_slot_length*16 + 1 : Short_preamble_slot_length*20);
z_fine = long_preamble_1 * long_preamble_2';

f_Fine_est = (-1/(2*pi*64*Parameters_struct.Ts)) * angle(z_fine);
Frame_After_Fine = Frame_After_Coarse .* exp(-j*2*pi*f_Fine_est*Parameters_struct.Ts*(0:Frame_length-1));

%% ========== INITIAL CHANNEL ESTIMATION (PREAMBLE-BASED) ==========
% This provides a coarse estimate; refined estimates come from pilots
Long_preamble_1_FFT = fftshift(fft(long_preamble_1));
Long_preamble_2_FFT = fftshift(fft(long_preamble_2));
H_est_preamble = 0.5 * (Long_preamble_1_FFT + Long_preamble_2_FFT) .* ...
                 conj(Parameters_struct.Long_preamble_slot_Frequency);

%% ========================================================================
%% ========== PAYLOAD 1 PROCESSING ==========
%% ========================================================================
RX_Payload_1_time = Frame_After_Fine(321:400);      % Extract payload 1 [1x80]
RX_Payload_1_no_CP = RX_Payload_1_time(17:80);     % Remove cyclic prefix [1x64]
RX_Payload_1_Frequency = fftshift(fft(RX_Payload_1_no_CP)); % FFTshifted [1xN_FFT]

%% --- Setup & Bob pilots (common across payloads) ---
Bob_pilot_idx     = Parameters_struct.Bob_pilot_idx(:).';       % ensure row vector
Bob_pilot_symbols = Parameters_struct.Bob_pilot_symbols(:).';   % row vector
Lch_Bob = Parameters_struct.Lch_Bob;
Np_Bob = length(Bob_pilot_idx);

% Sanity checks
assert(numel(Bob_pilot_idx) == numel(Bob_pilot_symbols), 'Bob pilot idx and symbols length mismatch');
assert(all(Bob_pilot_idx >= 1 & Bob_pilot_idx <= N_FFT), 'Bob_pilot_idx out of range');

% Extract pilot observations (column vectors)
Yp_Bob_1 = RX_Payload_1_Frequency(Bob_pilot_idx).';  % [Np x 1]
Xp_Bob   = Bob_pilot_symbols.';                     % [Np x 1]

% Convert pilot indices to centered frequency bins (for fftshifted representation)
% MATLAB indexing p=1..N corresponds to k_centered = (p-1) - N/2  (for even N)
k_centered_Bob = (Bob_pilot_idx - 1).' - (N_FFT/2);  % [Np x 1] column

% Build Vandermonde matrix (vectorized)
F_Bob = exp(-1j * 2*pi * (k_centered_Bob) * (0:(Lch_Bob-1)) / N_FFT);  % [Np x Lch_Bob]

% LS estimate (time-domain taps)
h_Bob_1 = F_Bob \ (Yp_Bob_1 ./ Xp_Bob);  % [Lch_Bob x 1]

% Reconstruct full frequency response with zero-padding and fftshift to match RX ordering
h_padded = [h_Bob_1; zeros(N_FFT - Lch_Bob, 1)];     % [N_FFT x 1]
H_Bob_recon_1 = fftshift(fft(h_padded)).';           % [1 x N_FFT]

%% --- Eve (Payload 1) ---
Eve_pilot_idx     = Parameters_struct.Eve_pilot_idx(:).';
Eve_pilot_symbols = Parameters_struct.Eve_pilot_symbols(:).';
Lch_Eve = Parameters_struct.Lch_Eve;
Np_Eve = length(Eve_pilot_idx);

assert(numel(Eve_pilot_idx) == numel(Eve_pilot_symbols), 'Eve pilot idx and symbols length mismatch');
assert(all(Eve_pilot_idx >= 1 & Eve_pilot_idx <= N_FFT), 'Eve_pilot_idx out of range');

Yp_Eve_1 = RX_Payload_1_Frequency(Eve_pilot_idx).';  % [Np x 1]
Xp_Eve   = Eve_pilot_symbols.';                     % [Np x 1]

k_centered_Eve = (Eve_pilot_idx - 1).' - (N_FFT/2);
F_Eve = exp(-1j * 2*pi * (k_centered_Eve) * (0:(Lch_Eve-1)) / N_FFT);

h_Eve_1 = F_Eve \ (Yp_Eve_1 ./ Xp_Eve);   % [Lch_Eve x 1]

h_padded_e = [h_Eve_1; zeros(N_FFT - Lch_Eve, 1)];
H_Eve_recon_1 = fftshift(fft(h_padded_e)).';        % [1 x N_FFT]

%% --- Zero-Forcing Equalization for Payload 1 ---
H_Bob_recon_1(abs(H_Bob_recon_1) < 1e-12) = 1e-12;
H_Eve_recon_1(abs(H_Eve_recon_1) < 1e-12) = 1e-12;

RX_Bob_Payload_1_Frequency_Eq = RX_Payload_1_Frequency ./ H_Bob_recon_1; % [1xN_FFT]
RX_Eve_Payload_1_Frequency_Eq = RX_Payload_1_Frequency ./ H_Eve_recon_1; % [1xN_FFT]

%% ========================================================================
%% ========== PAYLOAD 2 PROCESSING ==========
%% ========================================================================
RX_Payload_2_time = Frame_After_Fine(401:480);      % Extract payload 2 [1x80]
RX_Payload_2_no_CP = RX_Payload_2_time(17:80);     % Remove cyclic prefix [1x64]
RX_Payload_2_Frequency = fftshift(fft(RX_Payload_2_no_CP)); % FFTshifted [1xN_FFT]

%% --- Bob (Payload 2) ---
Yp_Bob_2 = RX_Payload_2_Frequency(Bob_pilot_idx).';  % [Np x 1]
h_Bob_2  = F_Bob \ (Yp_Bob_2 ./ Xp_Bob);             % [Lch_Bob x 1]

h_padded = [h_Bob_2; zeros(N_FFT - Lch_Bob, 1)];
H_Bob_recon_2 = fftshift(fft(h_padded)).';          % [1 x N_FFT]

%% --- Eve (Payload 2) ---
Yp_Eve_2 = RX_Payload_2_Frequency(Eve_pilot_idx).';  % [Np x 1]
h_Eve_2  = F_Eve \ (Yp_Eve_2 ./ Xp_Eve);             % [Lch_Eve x 1]

h_padded_e = [h_Eve_2; zeros(N_FFT - Lch_Eve, 1)];
H_Eve_recon_2 = fftshift(fft(h_padded_e)).';        % [1 x N_FFT]

%% --- Zero-Forcing Equalization for Payload 2 ---
H_Bob_recon_2(abs(H_Bob_recon_2) < 1e-12) = 1e-12;
H_Eve_recon_2(abs(H_Eve_recon_2) < 1e-12) = 1e-12;

RX_Bob_Payload_2_Frequency_Eq = RX_Payload_2_Frequency ./ H_Bob_recon_2; % [1xN_FFT]
RX_Eve_Payload_2_Frequency_Eq = RX_Payload_2_Frequency ./ H_Eve_recon_2; % [1xN_FFT]

%% ========== DIRECT CHANNEL ESTIMATE ==========
H_direct = zeros(N_FFT,1); 
H_direct(Bob_pilot_idx) = Yp_Bob_1 ./ Xp_Bob;

%% ========================================================================
%% ========== CHANNEL ANALYSIS AND CORRELATION ==========
%% ========================================================================

% Average channel estimates across both payload symbols
h_Bob_avg = 0.5 * (h_Bob_1 + h_Bob_2); % [Lch_Bob x 1]
h_Eve_avg = 0.5 * (h_Eve_1 + h_Eve_2); % [Lch_Eve x 1]

% Store padded versions for output
H_Bob = [h_Bob_avg.', zeros(1, N_FFT - Lch_Bob)]; % [1x64]
H_Eve = [h_Eve_avg.', zeros(1, N_FFT - Lch_Eve)]; % [1x64]

% Calculate normalized correlation coefficient
correlation_numerator = abs(h_Bob_avg' * h_Eve_avg);
correlation_denominator = norm(h_Bob_avg) * norm(h_Eve_avg);

if correlation_denominator < eps
    channel_correlation = 0;
else
    channel_correlation = correlation_numerator / correlation_denominator;
end

%% ========================================================================
%% ========== DATA DEMODULATION AND SYMBOL MAPPING ==========
%% ========================================================================

% Get data subcarrier indices (excluding pilots, DC, and virtual carriers)
if isfield(Parameters_struct, 'data_idx')
    data_idx = Parameters_struct.data_idx; % Use predefined indices
else
    % Auto-generate data indices if not provided (44 data subcarriers)
    % Exclude: DC(33), Virtual(1-6, 60-64), Bob pilots, Eve pilots
    all_idx = 1:64;
    exclude_idx = [1:6, 33, 60:64, Bob_pilot_idx, Eve_pilot_idx];
    data_idx = setdiff(all_idx, exclude_idx);
    data_idx = data_idx(1:min(N_data, length(data_idx))); % Limit to N_data
end

%% --- Bob's Data Demodulation ---
% Payload 1
RX_Bob_data_1_eq = RX_Bob_Payload_1_Frequency_Eq(data_idx); % [1 x N_data]
RX_Bob_data_1_demod = pskdemod(RX_Bob_data_1_eq, 4, pi/4); % QPSK demodulation

% Payload 2
RX_Bob_data_2_eq = RX_Bob_Payload_2_Frequency_Eq(data_idx); % [1 x N_data]
RX_Bob_data_2_demod = pskdemod(RX_Bob_data_2_eq, 4, pi/4);

%% --- Eve's Data Demodulation ---
% Payload 1
RX_Eve_data_1_eq = RX_Eve_Payload_1_Frequency_Eq(data_idx); % [1 x N_data]
RX_Eve_data_1_demod = pskdemod(RX_Eve_data_1_eq, 4, pi/4);

% Payload 2
RX_Eve_data_2_eq = RX_Eve_Payload_2_Frequency_Eq(data_idx); % [1 x N_data]
RX_Eve_data_2_demod = pskdemod(RX_Eve_data_2_eq, 4, pi/4);

%% Store equalized constellation points for visualization
RX_Bob_Payload_1 = RX_Bob_data_1_eq;
RX_Bob_Payload_2 = RX_Bob_data_2_eq;
RX_Eve_Payload_1 = RX_Eve_data_1_eq;
RX_Eve_Payload_2 = RX_Eve_data_2_eq;

%% ========================================================================
%% ========== BIT ERROR RATE CALCULATION ==========
%% ========================================================================

% Bob's BER - Handle different possible data storage formats
if isfield(Parameters_struct, 'data_Payload_Bob_1') && isfield(Parameters_struct, 'data_Payload_Bob_2')
    % Check if data is in a nested structure
    if isstruct(Parameters_struct.data_Payload_Bob_1)
        data_Bob_1_ref = Parameters_struct.data_Payload_Bob_1.data_Payload_1;
        data_Bob_2_ref = Parameters_struct.data_Payload_Bob_2.data_Payload_2;
    else
        data_Bob_1_ref = Parameters_struct.data_Payload_Bob_1;
        data_Bob_2_ref = Parameters_struct.data_Payload_Bob_2;
    end

    % Trim to match received data length
    len_match = min([length(data_Bob_1_ref), length(RX_Bob_data_1_demod)]);
    errors_Bob_1 = abs(sign(data_Bob_1_ref(1:len_match) - RX_Bob_data_1_demod(1:len_match)));

    len_match = min([length(data_Bob_2_ref), length(RX_Bob_data_2_demod)]);
    errors_Bob_2 = abs(sign(data_Bob_2_ref(1:len_match) - RX_Bob_data_2_demod(1:len_match)));

    Error_bits_Bob = sum(errors_Bob_1) + sum(errors_Bob_2);
    Total_bits_Bob = length(errors_Bob_1) + length(errors_Bob_2);
    BER_Bob = Error_bits_Bob / Total_bits_Bob;
else
    BER_Bob = NaN; % No reference data available
end

% Eve's BER
if isfield(Parameters_struct, 'data_Payload_Eve_1') && isfield(Parameters_struct, 'data_Payload_Eve_2')
    % Check if data is in a nested structure
    if isstruct(Parameters_struct.data_Payload_Eve_1)
        data_Eve_1_ref = Parameters_struct.data_Payload_Eve_1.data_Payload_1;
        data_Eve_2_ref = Parameters_struct.data_Payload_Eve_2.data_Payload_2;
    else
        data_Eve_1_ref = Parameters_struct.data_Payload_Eve_1;
        data_Eve_2_ref = Parameters_struct.data_Payload_Eve_2;
    end

    % Trim to match received data length
    len_match = min([length(data_Eve_1_ref), length(RX_Eve_data_1_demod)]);
    errors_Eve_1 = abs(sign(data_Eve_1_ref(1:len_match) - RX_Eve_data_1_demod(1:len_match)));

    len_match = min([length(data_Eve_2_ref), length(RX_Eve_data_2_demod)]);
    errors_Eve_2 = abs(sign(data_Eve_2_ref(1:len_match) - RX_Eve_data_2_demod(1:len_match)));

    Error_bits_Eve = sum(errors_Eve_1) + sum(errors_Eve_2);
    Total_bits_Eve = length(errors_Eve_1) + length(errors_Eve_2);
    BER_Eve = Error_bits_Eve / Total_bits_Eve;
else
    BER_Eve = NaN;
end

%% ========================================================================
%% ========== CHANNEL TAP LOGGING (FOR FIRST 250 FRAMES WITH LOW BER) ==========
%% ========================================================================


MAX_BOB_LOGS = 50;
MAX_EVE_LOGS = 50;
if length(Bob_tap_log) <= MAX_BOB_LOGS || length(Eve_tap_log) <= MAX_EVE_LOGS
%% ---- Bob logging (independent) ----
if length(Bob_tap_log) < MAX_BOB_LOGS && ~isnan(BER_Bob) && BER_Bob <= BER_THRESHOLD
    log_entry_Bob.frame_num = frame_counter;
    log_entry_Bob.taps = h_Bob_avg(:).';
    log_entry_Bob.BER = BER_Bob;
    log_entry_Bob.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');

    Bob_tap_log(end+1) = log_entry_Bob;
    fprintf('[Frame %d] Bob: BER = %.6f  --> Logged (%d/%d)\n', ...
             frame_counter, BER_Bob, length(Bob_tap_log), MAX_BOB_LOGS);
end

%% ---- Eve logging (independent) ----
if length(Eve_tap_log) < MAX_EVE_LOGS && ~isnan(BER_Eve) && BER_Eve <= BER_THRESHOLD
    log_entry_Eve.frame_num = frame_counter;
    log_entry_Eve.taps = h_Eve_avg(:).';
    log_entry_Eve.BER = BER_Eve;
    log_entry_Eve.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');

    Eve_tap_log(end+1) = log_entry_Eve;
    fprintf('[Frame %d] Eve: BER = %.6f  --> Logged (%d/%d)\n', ...
             frame_counter, BER_Eve, length(Eve_tap_log), MAX_EVE_LOGS);
end

%% ---- Periodic & Safe Auto-Save ----
if mod(frame_counter,50)==0 || length(Bob_tap_log)==MAX_BOB_LOGS || length(Eve_tap_log)==MAX_EVE_LOGS
    if ~isempty(Bob_tap_log), save('Bob_channel_taps_log.mat', 'Bob_tap_log'); end
    if ~isempty(Eve_tap_log), save('Eve_channel_taps_log.mat', 'Eve_tap_log'); end
end
end


% %% ---- Final Exit (only when both FULL) ----
% if length(Bob_tap_log) >= MAX_BOB_LOGS && length(Eve_tap_log) >= MAX_EVE_LOGS
%     fprintf('\n=== DONE: Bob & Eve logs FULLY FILLED. Stopping simulation. ===\n');
%     save('Bob_channel_taps_log.mat', 'Bob_tap_log');
%     save('Eve_channel_taps_log.mat', 'Eve_tap_log');
% 
% end


%% ========================================================================
%% ========== DEBUG VISUALIZATION ==========
%% ========================================================================

if strcmp(Debug_mode, 'on')
    figure('Name', 'Alice RX - Detailed Analysis', 'Position', [50 50 1600 900]);

    % ===== Row 1: Raw Signal =====
    subplot(3, 5, 1);
    plot(RX, '.', 'MarkerSize', 4);
    title('RX Raw Constellation');
    xlabel('Real'); ylabel('Imaginary');
    axis([-1.5 1.5 -1.5 1.5]); axis square; grid on;

    subplot(3, 5, 2);
    plot(real(RX), 'b', 'LineWidth', 0.5);
    title('I Component'); xlabel('Sample'); ylabel('Amplitude');
    axis([1 length(RX) -1.5 1.5]); grid on;

    subplot(3, 5, 3);
    plot(imag(RX), 'r', 'LineWidth', 0.5);
    title('Q Component'); xlabel('Sample'); ylabel('Amplitude');
    axis([1 length(RX) -1.5 1.5]); grid on;

    subplot(3, 5, 4);
    [Spectrum_waveform, Welch_Spectrum_frequency] = pwelch(RX, [], [], [], ...
        1/Parameters_struct.Ts, 'centered', 'power');
    plot(Welch_Spectrum_frequency/1e6, pow2db(Spectrum_waveform), 'LineWidth', 1.5);
    title('Power Spectral Density');
    xlabel('Frequency (MHz)'); ylabel('Power (dB)');
    grid on;

    subplot(3, 5, 5);
    text(0.1, 0.8, 'System Info', 'FontSize', 12, 'FontWeight', 'bold');
    text(0.1, 0.6, sprintf('FFT Size: %d', N_FFT), 'FontSize', 9);
    text(0.1, 0.5, sprintf('Data SC: %d', N_data), 'FontSize', 9);
    text(0.1, 0.4, sprintf('f_{coarse}: %.2f Hz', f_Coarse_est), 'FontSize', 9);
    text(0.1, 0.3, sprintf('f_{fine}: %.2f Hz', f_Fine_est), 'FontSize', 9);
    text(0.1, 0.2, sprintf('Packet idx: %d', idx), 'FontSize', 9);
    axis off;

    % ===== Row 2: Detection & Channels =====
    subplot(3, 5, 6);
    plot(1:length(M_n), M_n, 'b', 'LineWidth', 1.5);
    hold on;
    plot(1:length(M_n), Threshold_graph, 'r--', 'LineWidth', 2);
    title('Packet Detection Metric');
    xlabel('Sample'); ylabel('M_n');
    legend('M_n', 'Threshold', 'Location', 'best');
    axis([1 length(M_n) 0 1.2]); grid on;

    subplot(3, 5, 7);
    stem(0:Lch_Bob-1, abs(h_Bob_avg), 'b', 'LineWidth', 2, 'MarkerSize', 8);
    title('Bob Channel Taps |h_{Bob}|');
    xlabel('Tap Index'); ylabel('Magnitude');
    grid on;

    subplot(3, 5, 8);
    stem(0:Lch_Eve-1, abs(h_Eve_avg), 'm', 'LineWidth', 2, 'MarkerSize', 8);
    title('Eve Channel Taps |h_{Eve}|');
    xlabel('Tap Index'); ylabel('Magnitude');
    grid on;

    subplot(3, 5, 9);
    bar([norm(h_Bob_avg), norm(h_Eve_avg)]);
    set(gca, 'XTickLabel', {'Bob', 'Eve'});
    title('Channel Norms');
    ylabel('||h||'); grid on;

    subplot(3, 5, 10);
    text(0.1, 0.8, 'Channel Statistics', 'FontSize', 12, 'FontWeight', 'bold');
    text(0.1, 0.6, sprintf('Correlation: %.4f', channel_correlation), 'FontSize', 10);
    text(0.1, 0.5, sprintf('||h_{Bob}||: %.4f', norm(h_Bob_avg)), 'FontSize', 10);
    text(0.1, 0.4, sprintf('||h_{Eve}||: %.4f', norm(h_Eve_avg)), 'FontSize', 10);
    if channel_correlation < 0.3
        text(0.1, 0.2, 'Security: HIGH', 'FontSize', 10, 'Color', 'g', 'FontWeight', 'bold');
    elseif channel_correlation < 0.6
        text(0.1, 0.2, 'Security: MEDIUM', 'FontSize', 10, 'Color', [1 0.5 0], 'FontWeight', 'bold');
    else
        text(0.1, 0.2, 'Security: LOW', 'FontSize', 10, 'Color', 'r', 'FontWeight', 'bold');
    end
    axis off;

    % ===== Row 3: Equalization Results =====
    subplot(3, 5, 11);
    plot([RX_Bob_Payload_1, RX_Bob_Payload_2], '*b', 'MarkerSize', 8);
    title(sprintf('Bob Equalized (%d symbols)\nBER = %.4f', length([RX_Bob_Payload_1, RX_Bob_Payload_2]), BER_Bob));
    xlabel('Real'); ylabel('Imaginary');
    axis([-1.5 1.5 -1.5 1.5]); axis square; grid on;

    subplot(3, 5, 12);
    plot([RX_Eve_Payload_1, RX_Eve_Payload_2], '*m', 'MarkerSize', 8);
    title(sprintf('Eve Equalized (%d symbols)\nBER = %.4f', length([RX_Eve_Payload_1, RX_Eve_Payload_2]), BER_Eve));
    xlabel('Real'); ylabel('Imaginary');
    axis([-1.5 1.5 -1.5 1.5]); axis square; grid on;

    subplot(3, 5, 13);
    H_Bob_freq = fft([h_Bob_avg; zeros(N_FFT - Lch_Bob, 1)]);
    plot(1:N_FFT, abs(H_Bob_freq), 'b', 'LineWidth', 1.5);
    hold on;
    stem(Bob_pilot_idx, abs(H_Bob_freq(Bob_pilot_idx)), 'r', 'LineWidth', 2);
    title('Bob |H(f)|');
    xlabel('Subcarrier'); ylabel('Magnitude');
    legend('Estimated', 'Pilots', 'Location', 'best');
    grid on;

    subplot(3, 5, 14);
    H_Eve_freq = fft([h_Eve_avg; zeros(N_FFT - Lch_Eve, 1)]);
    plot(1:N_FFT, abs(H_Eve_freq), 'm', 'LineWidth', 1.5);
    hold on;
    stem(Eve_pilot_idx, abs(H_Eve_freq(Eve_pilot_idx)), 'r', 'LineWidth', 2);
    title('Eve |H(f)|');
    xlabel('Subcarrier'); ylabel('Magnitude');
    legend('Estimated', 'Pilots', 'Location', 'best');
    grid on;

    subplot(3, 5, 15);
    % Phase comparison
    phase_Bob = angle(H_Bob_freq);
    phase_Eve = angle(H_Eve_freq);
    plot(1:N_FFT, phase_Bob, 'b', 'LineWidth', 1.5);
    hold on;
    plot(1:N_FFT, phase_Eve, 'm', 'LineWidth', 1.5);
    title('Channel Phase Response');
    xlabel('Subcarrier'); ylabel('Phase (rad)');
    legend('Bob', 'Eve', 'Location', 'best');
    grid on;

    set(gcf, 'Units', 'centimeters', 'position', [1 2 55 28]);
end

%% ========== FUNCTION END ==========
end