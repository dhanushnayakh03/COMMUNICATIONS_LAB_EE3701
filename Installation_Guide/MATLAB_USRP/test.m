%% USRP Signal Reception
clear; clc; close all;

% Config
centerFreq = 2.45e9;
sampleRate = 10e6;
gain = 30;

% Setup Receiver
rx = comm.SDRuReceiver('Platform', 'B210', ...
    'CenterFrequency', centerFreq, ...
    'Gain', gain, ...
    'SampleRate', sampleRate, ...
    'SamplesPerFrame', 10000, ...
    'OutputDataType', 'double');

% Capture & Plot
fprintf('Receiving...\n');
data = rx();
release(rx);

% Plot PSD
[psd, f] = pwelch(data, [], [], [], sampleRate, 'centered');
plot(f/1e6, 10*log10(psd));
title('Power Spectral Density');
xlabel('Frequency (MHz)'); ylabel('Power (dB)');
grid on;