clear; close all; clc; j = 1i;
Global_Parameters_PLS;

%% Hardware Parameters for Eve
Mode = 'transmitRepeat'; % Select Mode
tx_object_Eve = sdrtx('Pluto', 'Gain', -10, 'CenterFrequency', Parameters_struct.CenterFrequency);

%% Button Setting
figure('Name', 'Eve TX', 'NumberTitle', 'off');
TransmittingDisplay = uicontrol('Style', 'text', 'Position', [55, 150, 155, 35], ...
    'String', 'Eve Transmitting', 'FontSize', 20, 'HorizontalAlignment', 'left', ...
    'BackgroundColor', [1.0 0.9 0.8]);
button = uicontrol; % Generate GUI button
set(button, 'String', 'Stop!', 'Position', [80 50 100 60]); % Add "Stop!" text
set(gcf, 'Units', 'centimeters', 'position', [11 3 7 6]); % Set the position of GUI

%% Load Eve's TX Signal
load('TX_signal_Eve'); % Load TX_signal_Eve

% transmitRepeat Mode requires data length >= 4096
TX_Hardware_Eve = repmat(TX_signal_Eve.', 5, 1); % [4860x1]

state = 1;

%% Main Transmission Loop
switch Mode
    case 'step'
        while (state == 1)
            step(tx_object_Eve, TX_Hardware_Eve);
            % ----- Button Behavior -----%
            set(button, 'Callback', 'setstate0_TX_Eve'); % Set the reaction of pushing button
            drawnow;
        end
        release(tx_object_Eve);
        
    case 'transmitRepeat'
        transmitRepeat(tx_object_Eve, TX_Hardware_Eve);
        % ----- Button Behavior -----%
        set(button, 'Callback', 'setstate0_TX_Eve'); % Set the reaction of pushing button
end

disp('Eve Transmission Started');
disp('Press Stop button to terminate transmission');