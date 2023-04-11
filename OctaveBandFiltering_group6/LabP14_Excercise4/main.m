%% Octave Band Filtering: Lab P-14: 4 Lab Excercise

close all
clear

%% 4.1a) 

% Initial Values
wc = 0.4*pi;            % Center frequency
L = 40;                 % Length of filter
n = linspace(0,L-1,L);  % Vector n (defined as 0<=n<L)

% Calculate impulse response of bandpass filter
h = 2/L * cos(wc.*n);

% Calculate frequency response
hf = fft(h,L);

% Perform fftshift
hf = fftshift(hf);

% Create frequency vector for plot
w = linspace(-pi,pi,L);

% Plot
figure(1)
clf
subplot(2,1,1)
plot(w,abs(hf))                 % Magnitude
title('|H(e^{jw})|)')
subtitle('w_c = 0.4\pi and L = 40')
xlabel('Frequency (radians)')
ylabel('Magnitude')
subplot(2,1,2)
plot(w,angle(hf)*180/pi)        % Phase
title('Phase')
xlabel('Frequency (radians)')
ylabel('Phase (degrees)')

%% 4.1b) 

% Look at magnitude plot, and plot a line at 0.5 to determine bandwidth of
% passband.
subplot(2,1,1)
hold on
yline(0.5)
hold off

% Pass band interects at w = 1.45 and 1.28886  = 0.16114. Bandwidth is then
% 0.16114 radians.

%% 4.1c)

% Develop BPF for L = 20 with same wc.
L = 20;                 % Length of filter
n = linspace(0,L-1,L);  % Vector n (defined as 0<=n<L)

% Calculate impulse response of bandpass filter
h1 = 2/L * cos(wc.*n);

% Create frequency vector for L = 20
w1 = linspace(-pi,pi,L);

% Develop BPF for L = 80 with same wc.
L = 80;                 % Length of filter
n = linspace(0,L-1,L);  % Vector n (defined as 0<=n<L)

% Calculate impulse response of bandpass filter
h2 = 2/L * cos(wc.*n);

% Create frequency vector for L = 20
w2 = linspace(-pi,pi,L);

% Calculate frequency response of each BPF (h1, h2)
hf1 = fft(h1);
hf2 = fft(h2);

% Perform fftshift for each frequency response
hf1 = fftshift(hf1);
hf2 = fftshift(hf2);

% Plot BPF with L = 20
figure(2)
clf
subplot(2,1,1)
plot(w1,abs(hf1))               % Magnitude
title('|H(e^{jw})|')
subtitle('w_c = 0.4\pi and L = 20')
xlabel('Frequency (radians)')
ylabel('Magnitude')
hold on
yline(0.5)                      % Create measure line to measure bandwidth
hold off
subplot(2,1,2)
plot(w1,angle(hf1)*180/pi)      % Phase
title('Phase')
xlabel('Frequency (radians)')
ylabel('Phase (degrees)')

% Pass band interects at w = 1.65347, 1.32278. Bandwidth is then 
% 0.33069 radians.

% Plot BPF with L = 80
figure(3)
clf
subplot(2,1,1)
plot(w2,abs(hf2))               % Magnitude
title('|H(e^{jw})|')
subtitle('w_c = 0.4\pi and L = 80')
xlabel('Frequency (radians)')
ylabel('Magnitude')
hold on
yline(0.5)                      % Create measure line to measure bandwidth
hold off
subplot(2,1,2)
plot(w2,angle(hf2)*180/pi)      % Phase
title('Phase')
xlabel('Frequency (radians)')
ylabel('Phase (degrees)')

% Pass band interects at w = 1.352078, 1.27254. Bandwidth is then 1/80*2pi =
% 0.0785 radians/sec.

% It can then be surmised that the bandwidth of the passband is inversely
% proportional to the filter length. When L is doubled, the bandwidth is
% halved and vice versa.

%% 4.2a)

% Define initial values
wc = 0.25*pi;
L = 41;
n = linspace(0,L-1,L);  % Vector n (defined as 0<=n<L)

% Calculate Hamming bandpass filter
hh = (0.54-0.46.*cos(2.*pi.*n./(L-1))).*cos(wc.*(n-(L-1)./2));

% Calculate fft
hhf = fft(hh,L);

% Perform fftshift
hhf = fftshift(hhf);

% Create frequency vector
w = linspace(-pi,pi,L);

% Plot magnitude and phase
figure(4)
clf
subplot(2,1,1)
plot(w,abs(hhf))               % Magnitude
title('|H(e^{jw})|')
xlabel('Frequency (radians)')
ylabel('Magnitude')
subplot(2,1,2)
plot(w,angle(hhf))             % Phase
title('Phase')
xlabel('Frequency (radians)')
ylabel('Phase (degrees)')

% Create vector with the frequencies of interest
wi = [0; 0.1*pi; 0.25*pi; 0.4*pi; 0.5*pi; 0.75*pi;];

%ind = zeros(length(wi));
mag = zeros(length(wi), 1);
phase = zeros(length(wi), 1);

for i = 1:length(wi) 
    ind = find(w == wi(i));
    mag(i) = abs(hhf(ind));
    phase(i) = angle(hhf(ind));
end 

table(wi, mag, phase)
%% 4.2b)

% Look at magnitude plot, and plot a line at 0.5 to determine bandwidth of
% passband.
subplot(2,1,1)
hold on
yline(0.5 * max(abs(hhf)));
hold off

% Pass band interects at w = 0.665 and 0.955. Bandwidth is then
% 0.29 radians, and the passband ranges from .212pi to .304pi.

% L = 21 
% Define initial values
L = 21;
n = linspace(0,L-1,L);  % Vector n (defined as 0<=n<L)

% Calculate Hamming bandpass filter
hh = (0.54-0.46.*cos(2.*pi.*n./(L-1))).*cos(wc.*(n-(L-1)./2));

% Calculate fft
hhf_L21 = fft(hh,L);

% Perform fftshift
hhf_L21 = fftshift(hhf_L21);

% Create frequency vector
w = linspace(-pi,pi,L);

% Plot magnitude and phase
figure(5)
clf
subplot(2,1,1)
plot(w,abs(hhf_L21))               % Magnitude
title('|H(e^{jw})|')
subtitle('L = 21')
xlabel('Frequency (radians)')
subplot(2,1,2)
plot(w,angle(hhf_L21))             % Phase
title('Phase')
xlabel('Frequency (radians)')
ylabel('Phase (degrees)')

subplot(2,1,1)
hold on
yline(0.5 * max(abs(hhf_L21)));
hold off

% Pass band interects at w = 0.484 and 1.1484. Bandwidth is then
% 0.6644 radians, and the passband ranges from .154pi to .366pi.

% L = 81 
% Define initial values
L = 81;
n = linspace(0,L-1,L);  % Vector n (defined as 0<=n<L)

% Calculate Hamming bandpass filter
hh = (0.54-0.46.*cos(2.*pi.*n./(L-1))).*cos(wc.*(n-(L-1)./2));

% Calculate fft
hhf_L81 = fft(hh,L);

% Perform fftshift
hhf_L81 = fftshift(hhf_L81);

% Create frequency vector
w = linspace(-pi,pi,L);

% Plot magnitude and phase
figure(6)
clf
subplot(2,1,1)
plot(w,abs(hhf_L81))               % Magnitude
title('|H(e^{jw})|')
subtitle('L = 81')
xlabel('Frequency (radians)')
subplot(2,1,2)
plot(w,angle(hhf_L81))             % Phase
title('Phase')
xlabel('Frequency (radians)')
ylabel('Phase (degrees)')

subplot(2,1,1)
hold on
yline(0.5 * max(abs(hhf_L81)));
hold off

% Pass band interects at w = 0.72575 and 0.86965. Bandwidth is then
% 0.1439 radians, and the passband ranges from .231pi to .2768pi.

% When L is halved, the bandwidth increases, but when L is doubled, the 
% bandwidth decreases. 

%% 4.2c)
% x[n] = 2 + 2cos(0.1*pi*n+pi/3) + cos(0.25*pi*n - pi/3)
% Breaking down x[n], we have the following frequencies and phases in the 
% individual parts of the signal.

% Frequency (radians)
% w = [0 , 0.1*pi, 0.25*pi]
% phase = [0, pi/3, -pi/3]

% Looking at the magnitude and phase of the impulse response of the BPF at 
% those specific frequencies, we have:
% mag = [0.08, 0.077701, 10.733]
% phase = [3.1416, -2.9883, -2.7585]

% We can then use these coefficients and phases to determine what the 
% output y[n] would be where y[n] is the convolution of x[n] and hh[n] 
% through multiplication of the magnitudes and addition of the phases.
% y[n] = (2*0.08) + (2*0.077701)cos(0.1*pi*n+pi/3 - 2.9883) + 
%        (1*10.733)cos(0.25*pi*n - pi/3 - 2.7585)

% The above gives us:
% y[n] = 0.16 + (0.1554)cos(0.1*pi*n+pi/3 - 2.9883) + 
%        (10.7330)cos(0.25*pi*n - pi/3 - 2.7585)

% Amplitudes are three signal components are 0.16 for w = 0, 0.1554 for w =
% 0.1*pi, and 10.733 for w = 0.25*pi. The amplitude of 10.733 indicates
% that this frequency is within the passband while the others are within
% the stopband of the filter.

%% 4.2d)
% This filter has an impulse response that centers on +-0.25*pi and the
% length of the filter determines the bandwidth of the passband. When the
% length of the filter is 41, the passband ranges from .211pi to .304pi and
% -.211pi to -.304pi.
% This indicates that signal components that have frequencies within this
% range will be able to pass through the filter, while those that are
% outside of this frequency range will be reduced/rejected.
