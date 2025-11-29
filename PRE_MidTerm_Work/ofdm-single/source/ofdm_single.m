%% Clearing variables
clear;close all;clc

%% OFDM MODULATOR
% Root of raised cosine filter parameters
alpha = 0.125;    % Rolloff factor
L     = 30;       % Truncated impulse length (-L/2:L/2)
Nc    = 10;       % Number of samples per symbol

% TX filter design
istrt = floor(L/2); 
n     = -istrt:1/Nc:istrt;
pt    = truncRRC(n, alpha, 0);

% RX filter design: p(t) = p(-t)
pmt = pt;

% Amplitudes of the paths
a = [0.8; 0.4; 0.4];
% Delays of the paths
tau = [0; -1; 2];

% Frequency-dispersive channel
ptau = sum(a.*truncRRC(n, alpha, tau), 1)/(sum(a));

% Computation of the number of useful subcarriers
N   = 64;                  
Na  = floor(N/2*(1-alpha));  
Nu  = 2*Na + 1;              
Nsc = N - Nu;                
Np  = 30;                    

% M-PSK constellation
M      = 4;                         
mb     = log2(M);                   
psklev = exp(1j*2*pi*(0:M-1)/M);    
nSym   = 1e3;                       

% Generation of data input and symbol mapping
symtx_data = randi([0,M-1], nSym, 1);

% Training sequence for equalizer
training_equalizer = randi([0,M-1], Nu, 1);
symtx = [training_equalizer; symtx_data];
[symtxg, ~] = togray(symtx, mb);      
ci          = psklev(symtxg+1);       

% Known training symbols for EQ
[tr_eq_gray, ~] = togray(training_equalizer,mb);
known_ci = psklev(tr_eq_gray+1);

% PN sequence for timing estimation
pn_length = N/2;
pn_sequence = -1 + 2.*randi([0, 1], 1, pn_length);
training_symbol = repmat(pn_sequence, 1, 2);
training_symbol_with_prefix = [training_symbol(end - Np + 1:end), training_symbol];

% Serial to parallel
r = rem(nSym, Nu);
if r~=0
    txt = sprintf('discarding last %d symbols', r);
    disp(txt);
    ci = ci(1:end-r);
end
ci_par = reshape(ci, Nu, []);
Nbl = size(ci_par,2);

% Virtual carrier insertion
ci_N = [ci_par(1:Na+1,:); zeros(Nsc, Nbl); ci_par(Na+2:end,:)];

% IDFT
a_N = N*ifft(ci_N);

% Cyclic prefix insertion
a_NT = [a_N(end-Np+1:end,:); a_N];

% P/S
a_NT = a_NT(:);

% Training symbol insertion
a_NT = [training_symbol_with_prefix.'; a_NT];

% TX filtering: upsample + filter
a_up = upsample(a_NT, Nc)/Nc;
txSig = filter(ptau, 1, a_up);

%% Radio parameters
SampleRate = 1e6;
SamplesPerRXFrame = 2^14; 
FramesToCollect = 5;
FrequencyOffset = 0;

%% Set up PlutoSDR
rx = sdrrx('Pluto','SamplesPerFrame',SamplesPerRXFrame,...
    'BasebandSampleRate',SampleRate,'OutputDataType','double');
tx = sdrtx('Pluto','Gain',0,...
    'BasebandSampleRate',SampleRate);
tx.CenterFrequency = tx.CenterFrequency + FrequencyOffset;
tx.transmitRepeat(txSig);

% Capture
saved = zeros(FramesToCollect*SamplesPerRXFrame,1);
ofe = 0;
for g=1:SamplesPerRXFrame:FramesToCollect*SamplesPerRXFrame
    [data1,len,of] = rx();
    saved(g:g+SamplesPerRXFrame-1) = data1;
    if of
        ofe = ofe + 1;
    end
end
fprintf('Overflow events: %d of %d\n',ofe,FramesToCollect)
rxSig = saved;

%% Receiver
rxdem = filter(pmt, 1, rxSig);
rx_ds = downsample(rxdem, Nc);
rx_sym = rx_ds(L+1:end);

%% Timing Sync Method A
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
figure, plot(M1), title("M1, Frame synchronization metric")

[picco,idx]=max(M1);
[piccok,idxk]=maxk(M1,50);
for indicemassimo = 1:50
    if abs(idxk(indicemassimo)-idx)>300
        secondomax = idxk(indicemassimo);
        break;
    end
end
if idx > secondomax
    change = idx;
    idx = secondomax;
    secondomax=change;
end
primomax = idx;
rx_sym = rx_sym(primomax+N:secondomax-Np-1);

%% S/P + remove CP
r = rem(length(rx_sym), N+Np);
a_NTr = rx_sym(1:end-r);
a_NTr = reshape(a_NTr, N+Np, []);
a_Nr  = a_NTr(Np+1:end,:);

%% FFT
ci_Nr = 1/N*(fft(a_Nr));

%% Remove unused carriers
ci_parr = [ci_Nr(1:Na+1,:); ci_Nr(end-Na+1:end,:)];

%% Equalization
Hstima = [ci_parr(:,1)].' ./ known_ci;
sym_da_eq_par = ci_parr(:,2:end);
equalized = sym_da_eq_par ./ Hstima.';
ci_eqserie = equalized(:);
cis = ci_parr(:);

figure, scatter(real(ci_eqserie), imag(ci_eqserie)), title("Equalized Symbols, (removed PN sequence)")
figure, scatter(real(cis), imag(cis)), title("NON Equalized Symbols")

%% Detection
dist_vec     = abs(psklev(:) - ci_eqserie.').^2;
[~, sym_idx] = min(dist_vec);
det_symg     = sym_idx - 1;
[det_sym, ~] = fromgray(det_symg, mb);

%% BER
l = min(size(det_sym,1) , 1000);
err = symtx_data(1:l) == det_sym(1:l);
correct_detection = sum(err);
bad_detect = l - correct_detection;
error_rate = bad_detect/l;
fprintf('BER: %.2f \n', error_rate);

%% PAPR
papr_Test = calculatePAPR(rxSig);
fprintf('PAPR: %.2f dB\n', papr_Test);

%% IQ Diagram
iq(1,:) = real(txSig);
iq(2,:) = imag(txSig);
figure, plot(iq(1,:),iq(2,:), "-"), title("IQ Diagram")

%% Power Spectrum
txSigNoF = filter(pt, 1, a_up);
[pxx, f] = pwelch(txSigNoF,[],[],[],SampleRate);
pxx_shifted = circshift(pxx, size(f,1)/2);
fplot = f- ones(size(f,1),1).*max(f)/2;
fplot = fplot/SampleRate;
figure, plot(fplot, 10*log10(pxx_shifted)), xlabel('Normalized Frequency f/Fs'), ylabel("PSD (dB)");
xticks(-0.5:0.1:0.5); grid on;
title('Power Spectrum');

%% Functions
function papr = calculatePAPR(signal)
    peakPower = max(abs(signal))^2;
    averagePower = mean(abs(signal))^2;
    papr = 10 * log10(peakPower / averagePower);
end

function pRRC=truncRRC(x,a,tau)
x = x - tau;
pRRC = (1-a)*sinc(x*(1-a));
pRRC = pRRC + a* sinc(a*x + 0.25).*cos(pi * (x + 0.25));
pRRC = pRRC + a* sinc(a*x - 0.25).*cos(pi * (x - 0.25));
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




