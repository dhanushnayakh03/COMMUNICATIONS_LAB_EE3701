%% Parameters for Physical Layer Security System
% Alice (RX), Tx1 (TX1), Tx2 (TX2)

Parameters_struct.CenterFrequency = 850e6; % 915 MHz
Parameters_struct.Bandwidth = 20e6; % 20 MHz
Parameters_struct.Ts = 64/Parameters_struct.Bandwidth; % 50 ns

%% OFDM Parameters
Parameters_struct.N_FFT = 64;
Parameters_struct.N_CP = 16;
Parameters_struct.N_data = 44; % Data subcarriers per symbol

%% Pilot Configuration for Tx1 and Tx2
% Tx1 uses pilot subcarriers: [8, 16, 24, 32, 40, 48, 56]
% Tx2 uses pilot subcarriers: [10, 18, 26, 34, 42, 50, 58]
% These are interleaved to avoid collision

Parameters_struct.Tx1_pilot_idx = [16,24,42,50]; % 7 pilots
Parameters_struct.Tx2_pilot_idx = [18,26,40,48]; % 7 pilots
Parameters_struct.Tx1_pilot_symbols = ones(1, 4); % BPSK [1,1,1,1,1,1,1]
Parameters_struct.Tx2_pilot_symbols = -ones(1, 4); % BPSK [-1,-1,-1,-1,-1,-1,-1]

%% Data subcarrier indices (avoiding pilots and DC)
% Total subcarriers: 1-64 (with fftshift: -32 to +31)
% After fftshift mapping: DC at 33, Virtual: [1:6, 60:64]
% Tx1 pilots at: [8, 16, 24, 32, 40, 48, 56]
% Tx2 pilots at: [10, 18, 26, 34, 42, 50, 58]
% Data carriers: remaining active subcarriers

% Data subcarriers (excluding DC=33, virtual carriers, and all pilots)
Parameters_struct.data_idx = [7:15,17,19:23, 25, 27:32, 34:39 , 41 , 43:47 , 49,51:59];
% This gives 34 data subcarriers

%% Load or Generate Payload Data
% For Tx1's data
if exist('data_Payload_1.mat', 'file')
    data_Payload_Tx1_1 = load('data_Payload_1'); 
else
    data_Payload_Tx1_1 = randi([0 3], 1, 34); % QPSK symbols
end

if exist('data_Payload_2.mat', 'file')
    data_Payload_Tx1_2 = load('data_Payload_2');
else
    data_Payload_Tx1_2 = randi([0 3], 1, 34);
end

% For Tx2's data
if exist('data_Payload_1.mat', 'file')
     data_Payload_Tx2_1 = load('data_Payload_1');
else
    data_Payload_Tx2_1 = randi([0 3], 1, 34);
end

if exist('data_Payload_2.mat', 'file')
    data_Payload_Tx2_2 = load('data_Payload_2');
else
    data_Payload_Tx2_2 = randi([0 3], 1, 34);
end

Parameters_struct.data_Payload_Tx1_1 = data_Payload_Tx1_1;
Parameters_struct.data_Payload_Tx1_2 = data_Payload_Tx1_2;
Parameters_struct.data_Payload_Tx2_1 = data_Payload_Tx2_1;
Parameters_struct.data_Payload_Tx2_2 = data_Payload_Tx2_2;



%% Short Preamble (same for both)
S_k = sqrt(13/6)*[0,0,1+1i,0,0,0,-1-1i,0,0,0,1+1i,0,0,0,-1-1i,0,0,0,-1-1i,0,0,0,1+1i,0,0,0,0,0,0,0,-1-1i,0,0,0,-1-1i,0,0,0,1+1i,0,0,0,1+1i,0,0,0,1+1i,0,0,0,1+1i,0,0];
virtual_subcarrier = zeros(1, Parameters_struct.N_FFT - length(S_k));
Short_preamble_slot_Frequency = [virtual_subcarrier(1:6), S_k, virtual_subcarrier(7:11)];
Parameters_struct.Short_preamble_slot_Frequency = Short_preamble_slot_Frequency;

%% Long Preamble (same for both - for timing sync)
L_k = [1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,0,1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1];
Long_preamble_slot_Frequency = [virtual_subcarrier(1:6), L_k, virtual_subcarrier(7:11)];
Parameters_struct.Long_preamble_slot_Frequency = Long_preamble_slot_Frequency;

%% Channel Estimation Parameters
Parameters_struct.Lch_Tx1 = 4; % Max channel taps for Tx1 (from 7 pilots)
Parameters_struct.Lch_Tx2 = 4; % Max channel taps for Tx2 (from 7 pilots)

disp('Physical Layer Security Parameters Loaded');
disp(['Tx1 Pilot Indices: ', num2str(Parameters_struct.Tx1_pilot_idx)]);
disp(['Tx2 Pilot Indices: ', num2str(Parameters_struct.Tx2_pilot_idx)]);
disp(['Number of Data Subcarriers: ', num2str(length(Parameters_struct.data_idx))]);