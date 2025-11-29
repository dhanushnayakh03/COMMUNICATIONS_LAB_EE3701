% Hardware_Tx_Dual.m
% Transmit two streams (perfectly synchronized) on B210 CH1 and CH2
% Uses transmitRepeat mode. Requires TX_signal_Dual.mat produced by OFDM_Tx_Dual.m

clear; close all; clc;
Global_Parameters_PLS;

%% Load dual signal
S = load('TX_signal_Dual.mat');
TX_signal_Dual = S.TX_signal_Dual;           % [N x 2] single
% Ensure length >= 4096 for transmitRepeat
if size(TX_signal_Dual,1) < 4096
    reps = ceil(4096 / size(TX_signal_Dual,1));
    TX_signal_Dual = repmat(TX_signal_Dual, reps, 1);
end

%% SDRu Transmitter (single device, two channels for tight sync)
fc   = Parameters_struct.CenterFrequency;
tx   = comm.SDRuTransmitter( ...
        'Platform',            'B210', ...
        'SerialNum',           '34414D2', ...   % your device
        'CenterFrequency',     fc, ...
        'Gain',                40, ...          % scalar applies to both channels
        'MasterClockRate',     30e6, ...
        'InterpolationFactor', 10, ...
        'ChannelMapping',      [1 2]);          % CH1 & CH2 simultaneously

%% UI
figure('Name','Dual TX (CH1 & CH2)','NumberTitle','off');
uicontrol('Style','text','Position',[40 140 220 35], 'String','Transmitting on CH1 & CH2', ...
          'FontSize',16,'HorizontalAlignment','left','BackgroundColor',[0.85 0.95 1]);
btn = uicontrol('String','Stop!','Position',[100 60 100 50], 'Callback','set(gcf,''UserData'',1);');
set(gcf,'Units','centimeters','position',[3 3 8 6]);

state = 0;
disp('Starting dual-channel transmitRepeat...');

try
    % Start continuous, synchronous playback on both channels
    while state == 0
        step(tx, TX_signal_Dual);   % [N x 2] single -> CH1, CH2
        drawnow;
        state = ~isempty(get(gcf,'UserData'));
    end
catch ME
    warning("TX error: %s", ME);
end

release(tx);
disp('Transmission stopped.');