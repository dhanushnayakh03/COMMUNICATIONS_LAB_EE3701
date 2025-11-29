clear; close all; clc; j = 1i;
Global_Parameters_PLS;

%% Hardware Parameters for Alice
Mode = 'transmitRepeat'; % Select Mode
tx_object_Alice = sdrtx('Pluto', 'Gain', -10, 'CenterFrequency', Parameters_struct.CenterFrequency);

%% Button Setting
figure('Name', 'Alice TX', 'NumberTitle', 'off');
TransmittingDisplay = uicontrol('Style', 'text', 'Position', [55, 150, 155, 35], ...
    'String', 'Alice Transmitting', 'FontSize', 20, 'HorizontalAlignment', 'left', ...
    'BackgroundColor', [0.8 0.9 1.0]);
button = uicontrol; % Generate GUI button
set(button, 'String', 'Stop!', 'Position', [80 50 100 60]); % Add "Stop!" text
set(gcf, 'Units', 'centimeters', 'position', [3 3 7 6]); % Set the position of GUI

%% Load Alice's TX Signal
load('TX_signal_Alice'); % Load TX_signal_Alice

% transmitRepeat Mode requires data length >= 4096
TX_Hardware_Alice = repmat(TX_signal_Alice.', 5, 1); % [4860x1]

state = 1;

%% Main Transmission Loop
switch Mode
    case 'step'
        while (state == 1)
            step(tx_object_Alice, TX_Hardware_Alice);
            % ----- Button Behavior -----%
            set(button, 'Callback', 'setstate0_TX_Alice'); % Set the reaction of pushing button
            drawnow;
        end
        release(tx_object_Alice);
        
    case 'transmitRepeat'
        transmitRepeat(tx_object_Alice, TX_Hardware_Alice);
        % ----- Button Behavior -----%
        set(button, 'Callback', 'setstate0_TX_Alice'); % Set the reaction of pushing button
end

disp('Alice Transmission Started');
disp('Press Stop button to terminate transmission');