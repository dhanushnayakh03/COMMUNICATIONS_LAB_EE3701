clear; close all; clc; j = 1i;
Global_Parameters_PLS;

global state;
state = 1;  % Running flag

%% GUI
f = figure('Name', 'Alice RX', 'NumberTitle', 'off', ...
    'Units', 'centimeters', 'Position', [1 2 55 28]);
uicontrol('Style', 'pushbutton', 'String', 'Stop!', ...
    'Position', [1475 15 100 60], 'Callback', @stop_rx);

%% Hardware Parameters for Alice (USRP B210)
rx_object_Alice = comm.SDRuReceiver(...
    'Platform', 'B210', ...
    'SerialNum', '3459F2C', ...
    'CenterFrequency',800e6, ...
    'Gain', 60, ...
    'MasterClockRate', 20e6, ...
    'DecimationFactor', 500, ...
    'ChannelMapping', 1, ...
    'SamplesPerFrame', 4000, ...
    'OutputDataType', 'double');


Ready_Time = 0;
scale = 1024;
Run_time_number = 1;


disp('Starting continuous reception... Press Stop! or Ctrl+C to stop.');

%% ---- Main Loop ----
while state
    try
        [data_rx_raw, len] = rx_object_Alice();

        if len > 0
            Run_time_number = Run_time_number + 1;

            if Run_time_number > Ready_Time
                % RX Raw Processing
                data_rx_scaled = data_rx_raw ./ scale;
                RX = data_rx_scaled.';

                % ---- Plot only every Nth frame ----
                if mod(Run_time_number, 10) == 0
                    sample_rate = Parameters_struct.Bandwidth;

                    subplot(3, 4, 1);
                    plot(RX, '.');
                    title('RX Raw - Constellation');
                    axis([-1.5 1.5 -1.5 1.5]); axis square; grid on;

                    subplot(3, 4, 2);
                    plot(real(RX)); title('RX - I Component');
                    axis([1 length(RX) -1e-2 1e-2]); grid on;

                    subplot(3, 4, 3);
                    plot(imag(RX)); title('RX - Q Component');
                    axis([1 length(RX) -1e-2 1e-2]); grid on;

                    % Power Spectral Density
                    [Spectrum_waveform, Welch_Spectrum_frequency] = pwelch(RX, [], [], [], sample_rate, 'centered', 'power');
                    subplot(3, 4, 4);
                    plot(Welch_Spectrum_frequency/1e6, pow2db(Spectrum_waveform));
                    title('Welch PSD');
                    xlabel('Frequency (MHz)'); ylabel('Power (dB)');
                    axis([-sample_rate/2e6, sample_rate/2e6, -100, -10]); grid on;
                end

                drawnow limitrate;  % Non-blocking update

                % ---- Processing Section ----
                [M_n, Threshold_graph, H_Bob, H_Eve, channel_correlation, ...
                    RX_Bob_Payload_1, RX_Bob_Payload_2, RX_Eve_Payload_1, RX_Eve_Payload_2, ...
                    BER_Bob, BER_Eve] = OFDM_RX_Alice(RX, Parameters_struct);

                % ---- Update other plots occasionally ----
                if mod(Run_time_number, 10) == 0
                    subplot(3, 4, 5);
                    plot(1:length(M_n), M_n, 1:length(M_n), Threshold_graph);
                    title('Packet Detection'); axis([1, length(M_n), 0, 1.2]);
                    legend('M_n', 'Threshold'); grid on;

                    subplot(3, 4, 6);
                    H_Bob_taps = H_Bob(1:Parameters_struct.Lch_Bob);
                    stem(0:Parameters_struct.Lch_Bob-1, abs(H_Bob_taps), 'b', 'LineWidth', 2);
                    title('Bob Channel Taps'); grid on;

                    subplot(3, 4, 7);
                    H_Eve_taps = H_Eve(1:Parameters_struct.Lch_Eve);
                    stem(0:Parameters_struct.Lch_Eve-1, abs(H_Eve_taps), 'm', 'LineWidth', 2);
                    title('Eve Channel Taps'); grid on;

                    subplot(3, 4, 8); cla;
                    text(0.1, 0.8, 'Channel Statistics:', 'FontSize', 12, 'FontWeight', 'bold');
                    text(0.1, 0.6, sprintf('Correlation: %.4f', channel_correlation), 'FontSize', 10);
                    text(0.1, 0.5, sprintf('||h_{Bob}||: %.4f', norm(H_Bob_taps)), 'FontSize', 10);
                    text(0.1, 0.4, sprintf('||h_{Eve}||: %.4f', norm(H_Eve_taps)), 'FontSize', 10);
                    text(0.1, 0.2, sprintf('Frame: %d', Run_time_number), 'FontSize', 10);
                    axis off;

                    subplot(3, 4, 9);
                    plot([RX_Bob_Payload_1, RX_Bob_Payload_2], '*b', 'MarkerSize', 8);
                    title(sprintf('Bob After Equalization\nBER = %.4f', BER_Bob));
                    axis([-1.5 1.5 -1.5 1.5]); axis square; grid on;

                    subplot(3, 4, 10);
                    plot([RX_Eve_Payload_1, RX_Eve_Payload_2], '*m', 'MarkerSize', 8);
                    title(sprintf('Eve After Equalization\nBER = %.4f', BER_Eve));
                    axis([-1.5 1.5 -1.5 1.5]); axis square; grid on;

                    subplot(3, 4, 11);
                    H_Bob_freq = fft(H_Bob);
                    plot(1:Parameters_struct.N_FFT, abs(H_Bob_freq), 'b');
                    title('Bob Channel - Freq Response'); grid on;

                    subplot(3, 4, 12);
                    H_Eve_freq = fft(H_Eve);
                    plot(1:Parameters_struct.N_FFT, abs(H_Eve_freq), 'm');
                    title('Eve Channel - Freq Response'); grid on;

                    drawnow limitrate;
                end
            end
        end

    catch ME
        fprintf(2, 'Error: %s\n', ME.message);
    end
end

release(rx_object_Alice);
close all;
disp('Alice Reception Complete');

%% --- Callback Function ---
function stop_rx(~, ~)
    global state;
    state = 0;
end
