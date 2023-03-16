%% Octave Band Filtering: Lab P-14: 3 Lab Excercise

%% 3.1

close all
clear

%% (a) Generate the impusle response of a length-25 bandpass filter with 
% w_c = 0.2pi, and plot h[n] with a stem plot.

% Initialize values
L = 25;                 % Length of 25
n = 0:L-1;              % n from 0 to 24
w_c = 0.2*pi;           % Center frequency

% Define impulse response using equation (2)
h = 2/L*cos(w_c.*n);

% Plot
figure
clf
stem(n,h)
title('Impulse response of a FIR filter with L = 25 and w_c = 0.2\pi')

%% (b) Compute the frequency response of the length-25 BPF from the 
% previous part and plot its magnitude and phase versus w over the range
% -2pi <= w <= 2pi.

% Compute frequency response using fft with L points
hf = fft(h,L);

% Generate frequency vector with 25 points
w = linspace(-2*pi,2*pi,L);

% Plot magnitude and phase
figure
clf
subplot(2,1,1)
plot(w,abs(hf))
title('Magnitude')
xlabel('Frequency (radians)')
ylabel('Magnitude')
subplot(2,1,2)
plot(w,angle(hf))
title('Phase')
xlabel('Frequency')
ylabel('Phase (radians)')

%% (c)  Use the magnitude response plot to describe the passband of the 
% BPF. For example, measure the bandwidth at a convenient point such as 
% 50% of the peak.

% Plot a line on the magnitude plot that corresponds to half of the peak
% value.
subplot(2,1,1)
hold on
yline(max(abs(hf))/2);

% After measuring the bandwidth using the plot, it's found that the
% passband of the BPF has a bandwidth of  5.52 - 4.292 = 1.228 centered
% around the center frequency of 0.2pi = 0.628


