% ===================== Receiver (Bob) =====================
clc ; close all ; 

rx_object_Alice = comm.SDRuReceiver(...
    'Platform', 'B210', ...
    'SerialNum', '3459F2C', ...
    'CenterFrequency', 850e6, ...
    'Gain', 40, ...
    'MasterClockRate', 20e6, ...
    'DecimationFactor', 500, ...
    'ChannelMapping', 2, ...
    'SamplesPerFrame', 4000, ...
    'OutputDataType', 'double');

fs = 20e6 / 500;   % Sampling rate after decimation
disp(['RX Sampling Rate: ', num2str(fs), ' Hz']);

release(rx_object_Alice);

figure;
while true
    [rxSig, len] = rx_object_Alice();
    if len > 0
        subplot(2,1,1);
        plot(real(rxSig));
        title('Received Signal (Real Part)');
        xlabel('Samples'); ylabel('Amplitude');
        grid on;

        subplot(2,1,2);
        pwelch(rxSig, [], [], [], fs, 'centered');
        title('Power Spectral Density');
        drawnow;
    end
end

release(rx_object_Alice);
