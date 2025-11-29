
function [M_n, Threshold_graph, H_Alice, ...
          RX_Alice_Payload_1, RX_Alice_Payload_2,...
          BER_Alice, H_direct] = OFDM_RX_Bob(RX, Parameters_struct)
%% OFDM_RX_Bob - Receiver for Bob in Physical Layer Security System
% 
% This function processes the superimposed signals from Alice and Eve,
% separates them using orthogonal pilot patterns, and estimates both channels.
%
% Channel Tap Logging Feature:
%   - Logs Alice and Eve channel taps (time domain) when BER ≈ 0
%   - Only logs first 250 frames
%   - Separate logs: 'Alice_channel_taps_log.mat' and 'Eve_channel_taps_log.mat'
%
% INPUTS:
%   RX                - Received signal [1 x N_samples]
%   Parameters_struct - System parameters structure
%
% OUTPUTS:
%   M_n                  - Packet detection metric
%   Threshold_graph      - Detection threshold visualization
%   H_Alice                - Alice's channel estimate (time domain, padded to 64)
%   H_Eve                - Eve's channel estimate (time domain, padded to 64)
%   channel_correlation  - Correlation coefficient between channels
%   RX_Alice_Payload_1/2   - Alice's equalized constellation points (44 each)
%   RX_Eve_Payload_1/2   - Eve's equalized constellation points (44 each)
%   BER_Alice              - Bit Error Rate for Alice's channel
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
    H_Alice = zeros(1, N_FFT);
    
   
    RX_Alice_Payload_1 = zeros(1, N_data);
    RX_Alice_Payload_2 = zeros(1, N_data);
    
    BER_Alice = 0.5;
    
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

%% --- Setup & Alice pilots (common across payloads) ---
Alice_pilot_idx     = Parameters_struct.Alice_pilot_idx(:).';       % ensure row vector
Alice_pilot_symbols = Parameters_struct.Alice_pilot_symbols(:).';   % row vector
Lch_Alice = Parameters_struct.Lch_Alice;
Np_Alice = length(Alice_pilot_idx);

% Sanity checks
assert(numel(Alice_pilot_idx) == numel(Alice_pilot_symbols), 'Alice pilot idx and symbols length mismatch');
assert(all(Alice_pilot_idx >= 1 & Alice_pilot_idx <= N_FFT), 'Alice_pilot_idx out of range');

% Extract pilot observations (column vectors)
Yp_Alice_1 = RX_Payload_1_Frequency(Alice_pilot_idx).';  % [Np x 1]
Xp_Alice   = Alice_pilot_symbols.';                     % [Np x 1]

% Convert pilot indices to centered frequency bins (for fftshifted representation)
% MATLAB indexing p=1..N corresponds to k_centered = (p-1) - N/2  (for even N)
k_centered_Alice = (Alice_pilot_idx - 1).' - (N_FFT/2);  % [Np x 1] column

% Build Vandermonde matrix (vectorized)
F_Alice = exp(-1j * 2*pi * (k_centered_Alice) * (0:(Lch_Alice-1)) / N_FFT);  % [Np x Lch_Alice]

% LS estimate (time-domain taps)
h_Alice_1 = F_Alice \ (Yp_Alice_1 ./ Xp_Alice);  % [Lch_Alice x 1]

% Reconstruct full frequency response with zero-padding and fftshift to match RX ordering
h_padded = [h_Alice_1; zeros(N_FFT - Lch_Alice, 1)];     % [N_FFT x 1]
H_Alice_recon_1 = fftshift(fft(h_padded)).';           % [1 x N_FFT]


%% --- Zero-Forcing Equalization for Payload 1 ---
H_Alice_recon_1(abs(H_Alice_recon_1) < 1e-12) = 1e-12;


RX_Alice_Payload_1_Frequency_Eq = RX_Payload_1_Frequency ./ H_Alice_recon_1; % [1xN_FFT]


%% ========================================================================
%% ========== PAYLOAD 2 PROCESSING ==========
%% ========================================================================
RX_Payload_2_time = Frame_After_Fine(401:480);      % Extract payload 2 [1x80]
RX_Payload_2_no_CP = RX_Payload_2_time(17:80);     % Remove cyclic prefix [1x64]
RX_Payload_2_Frequency = fftshift(fft(RX_Payload_2_no_CP)); % FFTshifted [1xN_FFT]

%% --- Alice (Payload 2) ---
Yp_Alice_2 = RX_Payload_2_Frequency(Alice_pilot_idx).';  % [Np x 1]
h_Alice_2  = F_Alice \ (Yp_Alice_2 ./ Xp_Alice);             % [Lch_Alice x 1]

h_padded = [h_Alice_2; zeros(N_FFT - Lch_Alice, 1)];
H_Alice_recon_2 = fftshift(fft(h_padded)).';          % [1 x N_FFT]


%% --- Zero-Forcing Equalization for Payload 2 ---
H_Alice_recon_2(abs(H_Alice_recon_2) < 1e-12) = 1e-12;

RX_Alice_Payload_2_Frequency_Eq = RX_Payload_2_Frequency ./ H_Alice_recon_2; % [1xN_FFT]


%% ========== DIRECT CHANNEL ESTIMATE ==========
H_direct = zeros(N_FFT,1); 
H_direct(Alice_pilot_idx) = Yp_Alice_1 ./ Xp_Alice;

%% ========================================================================
%% ========== CHANNEL ANALYSIS AND CORRELATION ==========
%% ========================================================================

% Average channel estimates across both payload symbols
h_Alice_avg = 0.5 * (h_Alice_1 + h_Alice_2); % [Lch_Alice x 1]

% Store padded versions for output
H_Alice = [h_Alice_avg.', zeros(1, N_FFT - Lch_Alice)]; % [1x64]



%% ========================================================================
%% ========== DATA DEMODULATION AND SYMBOL MAPPING ==========
%% ========================================================================

% Get data subcarrier indices (excluding pilots, DC, and virtual carriers)
if isfield(Parameters_struct, 'data_idx')
    data_idx = Parameters_struct.data_idx; % Use predefined indices
else
    % Auto-generate data indices if not provided (44 data subcarriers)
    % Exclude: DC(33), Virtual(1-6, 60-64), Alice pilots, Eve pilots
    all_idx = 1:64;
    exclude_idx = [1:6, 33, 60:64, Alice_pilot_idx, Eve_pilot_idx];
    data_idx = setdiff(all_idx, exclude_idx);
    data_idx = data_idx(1:min(N_data, length(data_idx))); % Limit to N_data
end

%% --- Alice's Data Demodulation ---
% Payload 1
RX_Alice_data_1_eq = RX_Alice_Payload_1_Frequency_Eq(data_idx); % [1 x N_data]
RX_Alice_data_1_demod = pskdemod(RX_Alice_data_1_eq, 4,pi/4); % QPSK demodulation

% Payload 2
RX_Alice_data_2_eq = RX_Alice_Payload_2_Frequency_Eq(data_idx); % [1 x N_data]
RX_Alice_data_2_demod = pskdemod(RX_Alice_data_2_eq, 4,pi/4);



%% Store equalized constellation points for visualization
RX_Alice_Payload_1 = RX_Alice_data_1_eq;
RX_Alice_Payload_2 = RX_Alice_data_2_eq;


%% ========================================================================
%% ========== BIT ERROR RATE CALCULATION ==========
%% ========================================================================

% Alice's BER - Handle different possible data storage formats
if isfield(Parameters_struct, 'data_Payload_Alice_1') && isfield(Parameters_struct, 'data_Payload_Alice_2')
    % Check if data is in a nested structure
    if isstruct(Parameters_struct.data_Payload_Alice_1)
        data_Alice_1_ref = Parameters_struct.data_Payload_Alice_1.data_Payload_1;
        data_Alice_2_ref = Parameters_struct.data_Payload_Alice_2.data_Payload_2;
    else
        data_Alice_1_ref = Parameters_struct.data_Payload_Alice_1;
        data_Alice_2_ref = Parameters_struct.data_Payload_Alice_2;
    end

    % Trim to match received data length
    len_match = min([length(data_Alice_1_ref), length(RX_Alice_data_1_demod)]);
    errors_Alice_1 = abs(sign(data_Alice_1_ref(1:len_match) - RX_Alice_data_1_demod(1:len_match)));

    len_match = min([length(data_Alice_2_ref), length(RX_Alice_data_2_demod)]);
    errors_Alice_2 = abs(sign(data_Alice_2_ref(1:len_match) - RX_Alice_data_2_demod(1:len_match)));

    Error_bits_Alice = sum(errors_Alice_1) + sum(errors_Alice_2);
    Total_bits_Alice = length(errors_Alice_1) + length(errors_Alice_2);
    BER_Alice = Error_bits_Alice / Total_bits_Alice;
else
    BER_Alice = NaN; % No reference data available
end



%% ========================================================================
%% ========== CHANNEL TAP LOGGING (FOR FIRST 250 FRAMES WITH LOW BER) ==========
%% ========================================================================


% MAX_Alice_LOGS = 50;
% MAX_EVE_LOGS = 50;
% if length(Alice_tap_log) <= MAX_Alice_LOGS || length(Eve_tap_log) <= MAX_EVE_LOGS
% %% ---- Alice logging (independent) ----
% if length(Alice_tap_log) < MAX_Alice_LOGS && ~isnan(BER_Alice) && BER_Alice <= BER_THRESHOLD
%     log_entry_Alice.frame_num = frame_counter;
%     log_entry_Alice.taps = h_Alice_avg(:).';
%     log_entry_Alice.BER = BER_Alice;
%     log_entry_Alice.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');
% 
%     Alice_tap_log(end+1) = log_entry_Alice;
%     fprintf('[Frame %d] Alice: BER = %.6f  --> Logged (%d/%d)\n', ...
%              frame_counter, BER_Alice, length(Alice_tap_log), MAX_Alice_LOGS);
% end
% 
% %% ---- Eve logging (independent) ----
% if length(Eve_tap_log) < MAX_EVE_LOGS && ~isnan(BER_Eve) && BER_Eve <= BER_THRESHOLD
%     log_entry_Eve.frame_num = frame_counter;
%     log_entry_Eve.taps = h_Eve_avg(:).';
%     log_entry_Eve.BER = BER_Eve;
%     log_entry_Eve.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');
% 
%     Eve_tap_log(end+1) = log_entry_Eve;
%     fprintf('[Frame %d] Eve: BER = %.6f  --> Logged (%d/%d)\n', ...
%              frame_counter, BER_Eve, length(Eve_tap_log), MAX_EVE_LOGS);
% end
% 
% %% ---- Periodic & Safe Auto-Save ----
% if mod(frame_counter,50)==0 || length(Alice_tap_log)==MAX_Alice_LOGS || length(Eve_tap_log)==MAX_EVE_LOGS
%     if ~isempty(Alice_tap_log), save('Alice_channel_taps_log.mat', 'Alice_tap_log'); end
%     if ~isempty(Eve_tap_log), save('Eve_channel_taps_log.mat', 'Eve_tap_log'); end
% end
% end


% %% ---- Final Exit (only when both FULL) ----
% if length(Alice_tap_log) >= MAX_Alice_LOGS && length(Eve_tap_log) >= MAX_EVE_LOGS
%     fprintf('\n=== DONE: Alice & Eve logs FULLY FILLED. Stopping simulation. ===\n');
%     save('Alice_channel_taps_log.mat', 'Alice_tap_log');
%     save('Eve_channel_taps_log.mat', 'Eve_tap_log');
% 
% end


%% ========================================================================
%% ========== DEBUG VISUALIZATION ==========
%% ========================================================================
if strcmp(Debug_mode, 'on')
    figure('Name', 'Alice RX - Detailed Analysis', 'Position', [50 50 1400 800]);

    % ===== Row 1: Raw Signal =====
    subplot(3, 4, 1);
    plot(RX, '.', 'MarkerSize', 4);
    title('RX Raw Constellation');
    xlabel('Real'); ylabel('Imaginary');
    axis([-1.5 1.5 -1.5 1.5]); axis square; grid on;

    subplot(3, 4, 2);
    plot(real(RX), 'b', 'LineWidth', 0.5);
    title('I Component'); xlabel('Sample'); ylabel('Amplitude');
    axis([1 length(RX) -1.5 1.5]); grid on;

    subplot(3, 4, 3);
    plot(imag(RX), 'r', 'LineWidth', 0.5);
    title('Q Component'); xlabel('Sample'); ylabel('Amplitude');
    axis([1 length(RX) -1.5 1.5]); grid on;

    subplot(3, 4, 4);
    [Spectrum_waveform, Welch_Spectrum_frequency] = pwelch(RX, [], [], [], ...
        1/Parameters_struct.Ts, 'centered', 'power');
    plot(Welch_Spectrum_frequency/1e6, pow2db(Spectrum_waveform), 'LineWidth', 1.5);
    title('Power Spectral Density');
    xlabel('Frequency (MHz)'); ylabel('Power (dB)');
    grid on;

    % ===== Row 2: Detection & Channel =====
    subplot(3, 4, 5);
    plot(1:length(M_n), M_n, 'b', 'LineWidth', 1.5);
    hold on;
    plot(1:length(M_n), Threshold_graph, 'r--', 'LineWidth', 2);
    title('Packet Detection Metric');
    xlabel('Sample'); ylabel('M_n');
    legend('M_n', 'Threshold', 'Location', 'best');
    axis([1 length(M_n) 0 1.2]); grid on;

    subplot(3, 4, 6);
    stem(0:Lch_Alice-1, abs(h_Alice_avg), 'b', 'LineWidth', 2, 'MarkerSize', 8);
    title('Alice Channel Taps |h_{Alice}|');
    xlabel('Tap Index'); ylabel('Magnitude');
    grid on;

    subplot(3, 4, 7);
    bar(norm(h_Alice_avg), 'b');
    set(gca, 'XTickLabel', {'Alice'});
    title('Channel Norm');
    ylabel('||h||'); grid on;

    subplot(3, 4, 8);
    text(0.1, 0.8, 'Alice Channel Stats', 'FontSize', 12, 'FontWeight', 'bold');
    text(0.1, 0.6, sprintf('||h_{Alice}||: %.4f', norm(h_Alice_avg)), 'FontSize', 10);
    text(0.1, 0.5, sprintf('Packet idx: %d', idx), 'FontSize', 10);
    text(0.1, 0.4, sprintf('f_{coarse}: %.2f Hz', f_Coarse_est), 'FontSize', 10);
    text(0.1, 0.3, sprintf('f_{fine}: %.2f Hz', f_Fine_est), 'FontSize', 10);
    axis off;

    % ===== Row 3: Equalization Results =====
    subplot(3, 4, 9);
    plot([RX_Alice_Payload_1, RX_Alice_Payload_2], '*b', 'MarkerSize', 8);
    title(sprintf('Alice Equalized (%d symbols)\nBER = %.4f', ...
        length([RX_Alice_Payload_1, RX_Alice_Payload_2]), BER_Alice));
    xlabel('Real'); ylabel('Imaginary');
    axis([-1.5 1.5 -1.5 1.5]); axis square; grid on;

    subplot(3, 4, 10);
    H_Alice_freq = fft([h_Alice_avg; zeros(N_FFT - Lch_Alice, 1)]);
    plot(1:N_FFT, abs(H_Alice_freq), 'b', 'LineWidth', 1.5);
    hold on;
    stem(Alice_pilot_idx, abs(H_Alice_freq(Alice_pilot_idx)), 'r', 'LineWidth', 2);
    title('Alice |H(f)|');
    xlabel('Subcarrier'); ylabel('Magnitude');
    legend('Estimated', 'Pilots', 'Location', 'best');
    grid on;

    subplot(3, 4, 11);
    plot(1:N_FFT, angle(H_Alice_freq), 'b', 'LineWidth', 1.5);
    title('Alice Phase Response');
    xlabel('Subcarrier'); ylabel('Phase (rad)');
    grid on;

    % General Info
    subplot(3, 4, 12);
    text(0.1, 0.8, 'System Info', 'FontSize', 12, 'FontWeight', 'bold');
    text(0.1, 0.6, sprintf('FFT Size: %d', N_FFT), 'FontSize', 9);
    text(0.1, 0.5, sprintf('Data SC: %d', N_data), 'FontSize', 9);
    axis off;

    set(gcf, 'Units', 'centimeters', 'position', [1 2 50 25]);
end

%% ========== FUNCTION END ==========
end