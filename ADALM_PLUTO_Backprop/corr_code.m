clear; clc;

% Load previously saved logs
load('Bob_channel_taps_log.mat', 'Bob_tap_log');
load('Eve_channel_taps_log.mat', 'Eve_tap_log');

% Extract only taps (Assuming each has field .taps)
Bob_taps = vertcat(Bob_tap_log.taps);   % (N x L) matrix
Eve_taps = vertcat(Eve_tap_log.taps);   % (M x L)

fprintf('Loaded %d Bob logs, %d Eve logs\n', size(Bob_taps,1), size(Eve_taps,1));

%% --- Correlation Analysis ---
% We compute average correlation across frames
corr_vals = zeros(min(size(Bob_taps,1),size(Eve_taps,1)),1);

for i = 1:length(corr_vals)
    hb = Bob_taps(i,:).';     % Bob channel vector
    he = Eve_taps(i,:).';     % Eve channel vector
    
    % Normalize channels (remove magnitude bias)
    hb_n = hb / norm(hb);
    he_n = he / norm(he);
    
    % Compute correlation (inner product magnitude)
    corr_vals(i) = (hb_n' * he_n);
end

%% --- Statistics ---
fprintf('\n=== Correlation Results ===\n');
fprintf('Mean Correlation |h_B · h_E| = %.4f\n', abs(mean(corr_vals)));
fprintf('Max Correlation                = %.4f\n', max(corr_vals));
fprintf('Min Correlation                = %.4f\n', min(corr_vals));

%% --- Plot for Visualization ---
figure; 
plot(abs(corr_vals),'LineWidth',1.5);
title('Bob-Eve Channel Correlation per Frame');
xlabel('Frame Index'); ylabel('|Correlation|');
grid on;
