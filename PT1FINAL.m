clc;
close all;

% Parameters
numSubcarriers = 120;       % Total number of subcarriers
user1Subcarriers = 1:120;   % Subcarriers assigned to User-1

% Initialize subcarrier mapping vector
subcarrierMapping = zeros(numSubcarriers, 1);

% Assign subcarriers to User-1
subcarrierMapping(user1Subcarriers) = 1;

% Display the subcarrier mapping
disp(subcarrierMapping);

% Set parameters for Zadoff-Chu sequence generation
sequenceLength = 120;  % Length of the ZC sequence
rootIdx = 5;           % Root index of the ZC sequence

% Generate the Zadoff-Chu sequence for User-1
n = 0:sequenceLength-1;
zcSequence = exp(-1i * pi * rootIdx * n .* (n + 1) / sequenceLength);

% Display the Zadoff-Chu sequence for User-1
disp(zcSequence);

% Parameters for OFDM waveform generation
user = 1;                % User index
fftSize = 2048;          % FFT size
cpLength = 144;          % Cyclic prefix length

% Generate the OFDM waveform for User-1
ofdmSignal = ifft(zcSequence, fftSize);   % Perform IFFT on the Zadoff-Chu sequence
ofdmSignal = [ofdmSignal(end-cpLength+1:end) ofdmSignal];   % Add cyclic prefix

% Display the OFDM waveform for User-1
disp(ofdmSignal);
plot(abs(ofdmSignal));

% Generate complex Gaussian random variables for each tap
h1 = (randn(1,1)+1j*randn(1,1))/sqrt(2);
h2 = (randn(1,1)+1j*randn(1,1))/sqrt(2);
h3 = (randn(1,1)+1j*randn(1,1))/sqrt(2);
h4 = (randn(1,1)+1j*randn(1,1))/sqrt(2);

% Combine taps to form broadband channel coefficients
h = [0.3*h1 0.8*h2 0.2*h3 0.5*h4];

% Display the generated coefficients
disp(h);

% Convolution of fading and transmitted signal
conv_op = conv(ofdmSignal, h);
plot(abs(conv_op));
M = length(conv_op);

%% Receiver Operations %%
% CP removal
L = 4;
rx_signal = conv_op((L+3):length(conv_op));

% Computing FFT
decoded_ofdm = fft(rx_signal);

% Zero-forcing equalization
H = fft(h, fftSize);
H_inv = 1 ./ H;

% QAM demodulation
qam_demod = qamdemod(decoded_ofdm, 2);
