
%% Full OFDM TX/RX with pilot-based LS channel estimation + LMMSE smoothing
% Paste into a single MATLAB .m file and run.

clear; close all; clc;

%% Parameters / TX design
alpha = 0.125;    % RRC rolloff
L     = 30;       % truncated impulse length (-L/2 : L/2)
Nc    = 2;       % samples per symbol (upsampling)
istrt = floor(L/2);
n     = -istrt:1/Nc:istrt;
pt    = truncRRC(n, alpha, 0);   % root RRC pulse (tx & rx)

% OFDM / carriers
N   = 64;
Na  = floor(N/2*(1-alpha));
Nu  = 2*Na + 1;        % number of active carriers (including DC)
Nsc = N - Nu;         % virtual carriers (unused)
Np  = 30;             % cyclic prefix length

% Modulation
M      = 4;
mb     = log2(M);
psklev = exp(1j*2*pi*(0:M-1)/M + 1j*pi*0.25);

% Data size
nSym   = 1e3;       % data symbols (excluding training)
symtx_data = randi([0,M-1], nSym, 1);

% Training (frequency-domain) - one OFDM symbol of pilots
training_equalizer = randi([0,M-1], Nu, 1);   % pilot symbols on active carriers
[tr_eq_gray, ~] = togray(training_equalizer, mb);
known_ci = psklev(tr_eq_gray+1);              % complex pilot values

% PN training/preamble used for timing (as before)
pn_length = N/2;
pn_sequence = -1 + 2.*randi([0, 1], 1, pn_length);
training_symbol = repmat(pn_sequence, 1, 2);
training_symbol_with_prefix = [training_symbol(end - Np + 1:end), training_symbol];


%% Build TX symbol stream
symtx = [training_equalizer; symtx_data];     % prepend training OFDM symbol
[symtxg, ~] = togray(symtx, mb);
ci          = psklev(symtxg+1);               % mapped complex symbols

% pad if required to reshape
r = rem(length(ci), Nu);
if r~=0
    fprintf('Discarding last %d symbols to make full OFDM blocks\n', r);
    ci = ci(1:end-r);
end

% Serial-to-parallel for active carriers
ci_par = reshape(ci, Nu, []);
Nbl = size(ci_par,2);   % number of OFDM blocks (including training)

% Insert virtual carriers (zeros) into N-length vector per block
ci_N = [ci_par(1:Na+1,:); zeros(Nsc, Nbl); ci_par(Na+2:end,:)];

% IDFT (note: using scaled IFFT to match your earlier code)
a_N = N * ifft(ci_N);

% Add cyclic prefix
a_NT = [a_N(end-Np+1:end,:); a_N];

% Serialize & prepend PN training symbol (time-domain training symbol for coarse sync)
a_NT = a_NT(:);
a_NT = [training_symbol_with_prefix.'; a_NT];

% TX pulse shaping (no synthetic channel applied)
a_up = upsample(a_NT, Nc)/Nc;
txSig = filter(pt, 1, a_up);   % use root-RRC pulse for TX shaping


%% Radio / SDR TX-RX or loopback capture
% YOU CAN USE PLUTO (as before). If you do not have Pluto attached,
% replace the SDRunet block with a simulated AWGN channel for offline testing.
usePluto = true;

if usePluto
    SampleRate = 1e6;
    SamplesPerRXFrame = 2^14;
    FramesToCollect = 5;
    FrequencyOffset = 0;

    rx = sdrrx('Pluto','SamplesPerFrame',SamplesPerRXFrame,...
        'BasebandSampleRate',SampleRate,'OutputDataType','double');
    tx = sdrtx('Pluto','Gain',0,'BasebandSampleRate',SampleRate);
    tx.CenterFrequency = tx.CenterFrequency + FrequencyOffset;
    tx.transmitRepeat(txSig);

    saved = zeros(FramesToCollect*SamplesPerRXFrame,1);
    ofe = 0;
    for g=1:SamplesPerRXFrame:FramesToCollect*SamplesPerRXFrame
        [data1,len,of] = rx();
        saved(g:g+SamplesPerRXFrame-1) = data1;
        if of
            ofe = ofe + 1;
        end
    end
    fprintf('Overflow events: %d of %d\n', ofe, FramesToCollect);
    rxSig = saved;
else
    % Loopback / simulated channel (if no Pluto) - simple AWGN channel
    SNRdB = 30;
    txSig_power = mean(abs(txSig).^2);
    noise_power = txSig_power / (10^(SNRdB/10));
    rxSig = txSig + sqrt(noise_power/2) * (randn(size(txSig)) + 1j*randn(size(txSig)));
    fprintf('Simulated loopback with AWGN SNR = %d dB\n', SNRdB);
end


%% RX: matched filtering, downsample, remove filter delay
rxdem = filter(pt, 1, rxSig);          % match filter with same root-RRC
rx_ds = downsample(rxdem, Nc);
rx_sym = rx_ds(L+1:end);               % remove filter transients as before


%% Timing synchronization (Method A) - unchanged logic
Ltrain = N/2;
finestra = 4000;
Pgrande = zeros(1, finestra);
Rf = zeros(1, finestra);
Mf = zeros(1, finestra);
M1 = zeros(1, finestra);

for d = 1: finestra
    r_asterisco = conj(rx_sym(d:d+Ltrain-1));
    r_base = rx_sym(d+Ltrain:d+2*Ltrain-1);
    Pgrande(d)= sum(r_asterisco .* r_base);
    r_N = rx_sym(d:d+N-1);
    Rf(d) = 0.5*sum(r_N.*conj(r_N));
    Mf(d) = (abs(Pgrande(d))^2) / (Rf(d)^2);
end
risultato = zeros(1, finestra);
for d= 1:finestra
    for k = 0 : Np
        if d-k >= 1
            risultato(d) = risultato(d) + Mf(d-k);
        end
    end
end
M1 = (1/(Np+1)) .* risultato;
figure, plot(M1), title("M1, Frame synchronization metric");

[picco,idx]=max(M1);
[piccok,idxk]=maxk(M1,50);
secondomax = idx + 400; % default fallback
for indicemassimo = 1:50
    if abs(idxk(indicemassimo)-idx) > 300
        secondomax = idxk(indicemassimo);
        break;
    end
end
if idx > secondomax
    change = idx;
    idx = secondomax;
    secondomax = change;
end
primomax = idx;
rx_sym = rx_sym(primomax+N : secondomax-Np-1);


%% Frame chop, remove CP, FFT
r = rem(length(rx_sym), N+Np);
a_NTr = rx_sym(1:end-r);
a_NTr = reshape(a_NTr, N+Np, []);
a_Nr  = a_NTr(Np+1:end,:);

ci_Nr = 1/N * (fft(a_Nr));    % N x Nblocks received OFDM spectra

% Extract active carriers only (same mapping as TX)
ci_parr = [ci_Nr(1:Na+1,:); ci_Nr(end-Na+1:end,:)];   % Nu x Nblocks


%% Pilot-based LS estimation (use first received OFDM block as training)
% Received pilot vector (active carriers) is the first column
Ypilot = ci_parr(:,1);               % Nu x 1 received pilots
H_ls = Ypilot ./ known_ci(:);        % LS estimate (Nu x 1)

% Print LS estimates (complex values)
fprintf('--- LS channel estimates (first %d active carriers) ---\n', min(20,length(H_ls)));
disp(H_ls(1:min(20,end)));


%% Estimate noise variance (simple heuristic)
% Here we estimate noise variance from the variance of the pilot residuals
% after removing mean; this is a heuristic.
resid = Ypilot - mean(Ypilot);  
sigma2_est = max(1e-8, mean(abs(resid).^2) * 0.1);


%% LMMSE smoothing across frequency (Wiener smoother)
% Frequency correlation model: Rhh(i,j) = rho^|i-j|
rho = 0.95;    % tune depending on channel frequency correlation (0..1)
idx_v = (1:Nu).';
Rhh = rho .^ (abs(idx_v - idx_v.'));      % Nu x Nu prior covariance
Ryy = Rhh + sigma2_est * eye(Nu);         % covariance + noise floor

W = Rhh / Ryy;                            % Wiener matrix
H_lmmse = W * H_ls;                       % Nu x 1 smoothed LMMSE estimate

% Print LMMSE-corrected estimates (complex values)
fprintf('--- LMMSE channel estimates (first %d active carriers) ---\n', min(20,length(H_lmmse)));
disp(H_lmmse(1:min(20,end)));


%% Extra Plots: Channel Estimates

% 1) LS Channel Estimate (coefficients + heatmap)
figure('Name','LS Channel Estimate','NumberTitle','off','Position',[200 200 1000 400]);

subplot(1,2,1);
plot(abs(H_ls),'-o','DisplayName','|H_{LS}|'); hold on;
plot(angle(H_ls),'-x','DisplayName','∠H_{LS}'); grid on;
xlabel('Active subcarrier index'); ylabel('Magnitude / Phase');
legend; title('LS channel (magnitude & phase)');

subplot(1,2,2);
imagesc(1,1:Nu,abs(H_ls(:).')); colormap jet; colorbar;
xlabel('Single pilot OFDM symbol'); ylabel('Subcarrier index');
title('LS channel heatmap (magnitude)');
set(gca,'YDir','normal');

% 2) LMMSE Channel Estimate (coefficients + heatmap)
figure('Name','LMMSE Channel Estimate','NumberTitle','off','Position',[200 200 1000 400]);

subplot(1,2,1);
plot(abs(H_lmmse),'-o','DisplayName','|H_{LMMSE}|'); hold on;
plot(angle(H_lmmse),'-x','DisplayName','∠H_{LMMSE}'); grid on;
xlabel('Active subcarrier index'); ylabel('Magnitude / Phase');
legend; title('LMMSE channel (magnitude & phase)');

subplot(1,2,2);
imagesc(1,1:Nu,abs(H_lmmse(:).')); colormap jet; colorbar;
xlabel('Single pilot OFDM symbol'); ylabel('Subcarrier index');
title('LMMSE channel heatmap (magnitude)');
set(gca,'YDir','normal');

%% Use LMMSE-corrected channel for equalization
sym_da_eq_par = ci_parr(:,2:end);    % exclude the pilot block, data blocks only
NdataBlocks = size(sym_da_eq_par, 2);

equalized = zeros(size(sym_da_eq_par));
for blk = 1:NdataBlocks
    equalized(:,blk) = sym_da_eq_par(:,blk) ./ H_lmmse;   % element-wise division
end

ci_eqserie = equalized(:);      % serial equalized constellation points
cis = ci_parr(:);

figure, scatter(real(ci_eqserie), imag(ci_eqserie), 10, '.'), title("Equalized Symbols (LMMSE-based)");
figure, scatter(real(cis), imag(cis), 10, '.'), title("All received carriers (incl. pilot)");


%% Detection & BER (compare to original transmitted data)
% detected constellation mapping
dist_vec     = abs(psklev(:) - ci_eqserie.').^2;
[~, sym_idx] = min(dist_vec);
det_symg     = (sym_idx - 1).';
[det_sym, ~] = fromgray(det_symg, mb);

% Compare to transmitted data (symtx_data)
l = min(length(det_sym), length(symtx_data));
correct_detection = sum(det_sym(1:l) == symtx_data(1:l));
error_rate = 1 - correct_detection / l;
fprintf('BER (symbol error rate) using LMMSE-estimated channel: %.6f\n', error_rate);


%% Optional: Display a small table of complex channel coefficients
fprintf('\nLS (first 10) vs LMMSE (first 10):\n');
for k = 1:min(10, Nu)
    fprintf('k=%2d  LS: % .4e + % .4e j   |  LMMSE: % .4e + % .4e j\n', ...
        k, real(H_ls(k)), imag(H_ls(k)), real(H_lmmse(k)), imag(H_lmmse(k)));
end


%% PAPR on transmitted signal (diagnostic)
papr_Test = calculatePAPR(txSig);
fprintf('TX PAPR: %.2f dB\n', papr_Test);


%% Supporting functions (local)
function papr = calculatePAPR(signal)
    peakPower = max(abs(signal))^2;
    averagePower = mean(abs(signal))^2;
    papr = 10 * log10(peakPower / averagePower);
end

function pRRC = truncRRC(x,a,tau)
    x = x - tau;
    pRRC = (1-a) * sinc(x * (1-a));
    pRRC = pRRC + a * sinc(a*x + 0.25).*cos(pi * (x + 0.25));
    pRRC = pRRC + a * sinc(a*x - 0.25).*cos(pi * (x - 0.25));
end

function [gd, gb] = togray(in, len)
    if ischar(in)
        in = bin2dec(in);
    end
    gd = bitxor(in, bitshift(in,-1));
    gb = dec2bin(gd, len);
end

function [de, bi] = fromgray(in, len)
    if ischar(in)
        in = bin2dec(in);
    end
    in = in(:);
    n_syms = length(in);
    b = zeros(n_syms,len);
    b(:,1) = bitget(in, repmat(len,[n_syms,1]));
    for sh = 2 : len
        in_loc = bitget(in, repmat(len-sh+1,[n_syms,1]));
        b(:,sh) = bitxor(b(:,sh-1), in_loc);
    end
    bi = num2str(b);
    de = bin2dec(bi);
end
