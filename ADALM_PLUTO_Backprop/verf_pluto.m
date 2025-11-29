clc ; clear all ; 


rxPluto = sdrrx('Pluto', ...
    'CenterFrequency', 850e6, ...
    'GainSource', 'Manual', ...
    'Gain', 40, ...
    'BasebandSampleRate', 66.7e3, ...   % ≈ 66.7 kHz
    'SamplesPerFrame', 4000, ...
    'OutputDataType', 'double');

fs = 66.7e3;
disp(['RX Sampling Rate: ', num2str(fs), ' Hz']);

release(rxPluto);

figure;
while true
    [rxSig, len] = rxPluto();
    if len > 0
        subplot(2,1,1);
        plot(real(rxSig(1:100)));
        title('Received Signal (Real Part)');
        xlabel('Samples'); ylabel('Amplitude');
        grid on;

        subplot(2,1,2);
        pwelch(rxSig, [], [], [], fs, 'centered', 'power');
        title('Power Spectral Density');
        drawnow;
    end
end

release(rxPluto);
