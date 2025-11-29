close all; clear; clc; j = 1i;
Global_Parameters_PLS;

%% OFDM TX for Bob (Legitimate Transmitter)
N_FFT = Parameters_struct.N_FFT;

%% Short Preamble
Short_preamble_slot_Time = ifft(ifftshift(Parameters_struct.Short_preamble_slot_Frequency));
Short_preamble = repmat(Short_preamble_slot_Time(1:16), 1, 10); % [1x160]

%% Long Preamble
Long_preamble_slot_Time = ifft(ifftshift(Parameters_struct.Long_preamble_slot_Frequency));
Long_preamble = [Long_preamble_slot_Time(33:64), Long_preamble_slot_Time, Long_preamble_slot_Time]; % [1x160]

%% Payload Symbol 1 with Bob's Orthogonal Pilots
M = 4; % QPSK
data1_mod = pskmod(Parameters_struct.data_Payload_Bob_1.data_Payload_1, M, pi/4); % [1x34]

% Initialize frequency domain symbol
Payload1_Frequency = zeros(1, N_FFT);

% Insert virtual subcarriers (zeros at edges and DC)
% DC at position 33

% Insert Bob's pilots at his designated positions
Bob_pilot_idx = Parameters_struct.Bob_pilot_idx;
Bob_pilot_symbols = Parameters_struct.Bob_pilot_symbols;
Payload1_Frequency(Bob_pilot_idx) = Bob_pilot_symbols;

% Insert data on data subcarriers
data_idx = Parameters_struct.data_idx;
Payload1_Frequency(data_idx) = data1_mod;

% Transform to time domain
data1_slot_Time = ifft(ifftshift(Payload1_Frequency)); % [1x64]
data1_TX_payload = [data1_slot_Time(49:64), data1_slot_Time]; % Add CP [1x80]

%% Payload Symbol 2 with Bob's Orthogonal Pilots
data2_mod = pskmod(Parameters_struct.data_Payload_Bob_2.data_Payload_2, M, pi/4); % [1x34]

Payload2_Frequency = zeros(1, N_FFT);
Payload2_Frequency(Bob_pilot_idx) = Bob_pilot_symbols;
Payload2_Frequency(data_idx) = data2_mod;

data2_slot_Time = ifft(ifftshift(Payload2_Frequency));
data2_TX_payload = [data2_slot_Time(49:64), data2_slot_Time]; % [1x80]

%% Frame Combination
Frame = [Short_preamble, Long_preamble, data1_TX_payload, data2_TX_payload]; % [1x480]

%% Oversampling
OVR = 2;
Frame_OVR = zeros(1, length(Frame) * OVR);
for i = 1:length(Frame)
    Frame_OVR(2*i-1:2*i) = Frame(i);
end

%% Root Raised Cosine Filter
rolloff = 0.5;
L = 6;
RRC = rcosdesign(rolloff, L, OVR, 'sqrt');

%% TX Signal for Bob
TX_signal_Bob = conv(Frame_OVR, RRC); % [1x972]

%% Visualization
figure('Name', 'Bob TX Signal', 'Position', [100 100 1200 600]);

subplot(2,3,1);
plot(real(TX_signal_Bob));
title('Bob TX - Real Part');
xlabel('Sample'); ylabel('Amplitude');
grid on;

subplot(2,3,2);
plot(imag(TX_signal_Bob));
title('Bob TX - Imaginary Part');
xlabel('Sample'); ylabel('Amplitude');
grid on;

subplot(2,3,3);
plot(TX_signal_Bob, '.');
title('Bob TX - Constellation');
xlabel('Real'); ylabel('Imaginary');
axis equal; grid on;

subplot(2,3,4);
plot(abs(fftshift(fft(Payload1_Frequency))));
title('Bob Payload 1 - Frequency Domain');
xlabel('Subcarrier'); ylabel('Magnitude');
hold on;
stem(Bob_pilot_idx, abs(Bob_pilot_symbols), 'r', 'LineWidth', 2);
legend('Data + Pilots', 'Bob Pilots');
grid on;

subplot(2,3,5);
[pxx, f] = pwelch(TX_signal_Bob, [], [], [], 1/Parameters_struct.Ts, 'centered');
plot(f/1e6, pow2db(pxx));
title('Bob TX - Power Spectral Density');
xlabel('Frequency (MHz)'); ylabel('Power (dB)');
grid on;

subplot(2,3,6);
text(0.1, 0.7, 'Bob TX Configuration:', 'FontSize', 12, 'FontWeight', 'bold');
text(0.1, 0.5, ['Pilot Indices: ', num2str(Bob_pilot_idx)], 'FontSize', 10);
text(0.1, 0.4, ['Pilot Symbols: All +1'], 'FontSize', 10);
text(0.1, 0.3, ['Data Subcarriers: ', num2str(length(data_idx))], 'FontSize', 10);
text(0.1, 0.2, ['Frame Length: ', num2str(length(TX_signal_Bob)), ' samples'], 'FontSize', 10);
axis off;

%% Save TX Signal
save('TX_signal_Bob.mat', 'TX_signal_Bob');
disp('Bob TX Signal Generated and Saved');
disp(['Signal Length: ', num2str(length(TX_signal_Bob)), ' samples']);