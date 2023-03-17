%% Octave Band Filtering: Lab P-14: 3 Lab Excercise

close all
clear

%% 3.1a)
% Generate the impulse response of a length-25 bandpass filter with 
% w_c = 0.2pi, and plot h[n] with a stem plot.

% Initialize values
L = 25;                 % Length of 25
n = 0:L-1;              % n from 0 to 24
w_c = 0.2*pi;           % Center frequency

% Define impulse response using equation (2)
h = 2/L*cos(w_c.*n);

% Plot
figure(1)
clf
stem(n,h)
title('Impulse response of a FIR filter with L = 25 and w_c = 0.2\pi')

%% 3.1b) 
% Compute the frequency response of the length-25 BPF from the 
% previous part and plot its magnitude and phase versus w over the range
% -2pi <= w <= 2pi.

% Compute frequency response using fft with L points
hf = fft(h,L);

% Generate frequency vector with 25 points
w = linspace(-2*pi,2*pi,L);

% Plot magnitude and phase
figure(2)
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

%% 3.1c)  
% Use the magnitude response plot to describe the passband of the 
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

%% 3.2a)  
% Plot the magnitude response of the 25-point BPF filter defined in 
% (2). For this plot, it is sufficient to use a horizontal frequency axis 
% that extends from -pi to pi.

% Plot 
figure(3)
plot()

%% 3.2b) 
% Use the stem function to place vertical markers at several points
% on the frequency response at frequencies {0, +-0.2pi, +-0.5pi}.

hold on
stem(pi*[-0.5,-0.2,0,0.2,0.5],0.9*ones(1,5),'r.')
hold on

%% 3.4a) 
% Generate x[n] in a vector called xx.

%% 3.4b) 
% Filter xx to obtain yy.

%% 3.4c) 
% Make stem plots of xx and yy.

%% 3.4d) 
% Use the frequency response to validate that the output signal has 
% the correct magnitude and phase in each of the three regions where the 
% input has different frequencies.

%% 3.4e) 
% Comment: observe the transient effect at the transitions (n D 200 
% and n D 400). This is due to the start-up of the FIR filter as it 
% encounters a new sinusoid. How long does the transient last (in samples)?

