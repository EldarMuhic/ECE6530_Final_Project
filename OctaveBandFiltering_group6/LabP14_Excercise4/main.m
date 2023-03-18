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
plot(w,abs(hf))       % Magnitude
title('|H(e^{jw})|)')
subtitle('w_c = 0.4\pi and L = 40')
xlabel('Frequency (radians)')
ylabel('Magnitude')
subplot(2,1,2)
plot(w,angle(hf)*180/pi)     % Phase
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
title('Frequency response of BPF with w_c = 0.4\pi with L = 20')
clf
subplot(2,1,1)
plot(w1,abs(hf1))       % Magnitude
title('|H(e^{jw})|')
subtitle('w_c = 0.4\pi and L = 20')
xlabel('Frequency (radians)')
hold on
yline(0.5)              % Create measure line to measure bandwidth
hold off
subplot(2,1,2)
plot(w1,angle(hf1)*180/pi)     % Phase
title('Phase')
xlabel('Frequency (radians)')
ylabel('Phase (degrees)')

% Pass band interects at w = 1.65347, 1.32278. Bandwidth is then 
% 0.33069 radians.

% Plot BPF with L = 80
figure(3)
title('Frequency response of BPF with w_c = 0.4\pi with L = 80')
clf
subplot(2,1,1)
plot(w2,abs(hf2))       % Magnitude
title('|H(e^{jw})|')
subtitle('w_c = 0.4\pi and L = 80')
xlabel('Frequency (radians)')
hold on
yline(0.5)              % Create measure line to measure bandwidth
hold off
subplot(2,1,2)
plot(w2,angle(hf2)*180/pi)     % Phase
title('Phase')
xlabel('Frequency (radians)')
ylabel('Phase (degrees)')

% Pass band interects at w = 1.352078, 1.27254. Bandwidth is then 1/80*2pi =
% 0.0785 radians/sec.

% It can then be surmised that the bandwidth of the passband is inversely
% proportional to the filter length. When L is doubled, the bandwidth is
% halved and vice versa.

%% 4.2a)

%% 4.2b)

%% 4.2c)

%% 4.2d)