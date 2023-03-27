%% Octave Band Filtering: Lab P-14: 5 Lab Excercise

%% 5.1a)

% Initialize A4
A4 = 440;

% Create vectors to save upper and lower frequencies of octaves
upper_freq = zeros(5,1);
lower_freq = zeros(5,1);

% Calculate C4 and other frequencies of C notes from octaves 2-6 (lower
% frequencies)
lower_freq(3) = 440 * 2^(-9/12);        % Middle C (C4)
lower_freq(4) = lower_freq(3) * 2;      % C5
lower_freq(5) = lower_freq(4) * 2;      % C6
lower_freq(2) = lower_freq(3) * 2^(-1); % C3
lower_freq(1) = lower_freq(2) * 2^(-1); % C2

% Calculate B notes of respective octaves to determine upper frequencies.
upper_freq(2) = lower_freq(3) * 2^(-1/12);  % B3
upper_freq(3) = upper_freq(2) * 2;          % B4
upper_freq(4) = upper_freq(3) * 2;          % B5
upper_freq(5) = upper_freq(4) * 2;          % B6
upper_freq(1) = upper_freq(2) * 2^-1;       % B2

disp('Frequencies in Hz')
table(lower_freq,upper_freq)

disp('Frequencies in radians')
table(lower_freq*2*pi,upper_freq*2*pi)

%% 5.1b)
% Compute the center frequencies of each octave (Octaves 2-6)
center_freq = zeros(5,1);

for i = 1:5
    center_freq(i) = (upper_freq(i)+lower_freq(i))/2;
end

%% 5.2a)
% Our strategy for determining beta will be performed by using MATLAB's max
% function to measure the peak value of the unscaled frequency response,
% and then have MATLAB compute beta to scale the peak to be one. The code
% will look similar to the following.

% beta = 1/max(abs(hf_unscaled));
% hf_scaled = beta*hf_unscaled;

%% 5.2b)
% Sampling frequency (Hz)
fs = 8000;

% Calculate w_c values from center_freq
w_c = center_freq./fs.*2;

%% Design of BPF for octave 2.
L = 197;                         % Bandwidth = 134.33Hz-55.45Hz
n = linspace(0,L-1,L);

% Calculate impulse response
h = (0.54-0.46.*cos(2.*pi.*n./(L-1))).*cos(w_c(1)*pi.* ...
    (n-(L-1)./2));

% Calculate frequency impulse response
hf = fft(h,L);

% Perform fftshift
hf_unscaled = fftshift(hf);

% Scale using calculated beta
beta = 1/max(abs(hf_unscaled));
hf_scaled2 = beta*hf_unscaled;

% Create frequency vector
w2 = linspace(-pi,pi,L);

% figure(1)
% clf
% plot(w2,abs(hf_scaled2))

%% Design of BPF for octave 3.
L = 125;                         % Bandwidth = 246.7538Hz-132.5442Hz
n = linspace(0,L-1,L);

% Calculate impulse response
h = (0.54-0.46.*cos(2.*pi.*n./(L-1))).*cos(w_c(2)*pi.* ...
    (n-(L-1)./2));

% Calculate frequency impulse response
hf = fft(h,L);

% Perform fftshift
hf_unscaled = fftshift(hf);

% Scale using calculated beta
beta = 1/max(abs(hf_unscaled));
hf_scaled3 = beta*hf_unscaled;

% Create frequency vector
w3 = linspace(-pi,pi,L);

% figure(1)
% clf
% plot(w3,abs(hf_scaled3))

%% Design of BPF for octave 4.
L = 55;                         % Bandwidth = 537.94Hz-225.62Hz
n = linspace(0,L-1,L);

% Calculate impulse response
h = (0.54-0.46.*cos(2.*pi.*n./(L-1))).*cos(w_c(3)*pi.* ...
    (n-(L-1)./2));

% Calculate frequency impulse response
hf = fft(h,L);

% Perform fftshift
hf_unscaled = fftshift(hf);

% Scale using calculated beta
beta = 1/max(abs(hf_unscaled));
hf_scaled4 = beta*hf_unscaled;

% Create frequency vector
w4 = linspace(-pi,pi,L);

% figure(1)
% clf
% plot(w4,abs(hf_scaled4))

%% Design of BPF for octave 5.
L = 35;                         % Bandwidth = 1018Hz-551.3127Hz
n = linspace(0,L-1,L);

% Calculate impulse response
h = (0.54-0.46.*cos(2.*pi.*n./(L-1))).*cos(w_c(4)*pi.* ...
    (n-(L-1)./2));

% Calculate frequency impulse response
hf = fft(h,L);

% Perform fftshift
hf_unscaled = fftshift(hf);

% Scale using calculated beta
beta = 1/max(abs(hf_unscaled));
hf_scaled5 = beta*hf_unscaled;

% Create frequency vector
w5 = linspace(-pi,pi,L);

% figure(1)
% clf
% plot(w5,abs(hf_scaled5))

%% Design of BPF for octave 6.
L = 15;                         % Bandwidth = 2158.1Hz-1031.3Hz
n = linspace(0,L-1,L);

% Calculate impulse response
h = (0.54-0.46.*cos(2.*pi.*n./(L-1))).*cos(w_c(5)*pi.* ...
    (n-(L-1)./2));

% Calculate frequency impulse response
hf = fft(h,L);

% Perform fftshift
hf_unscaled = fftshift(hf);

% Scale using calculated beta
beta = 1/max(abs(hf_unscaled));
hf_scaled6 = beta*hf_unscaled;

% Create frequency vector
w6 = linspace(-pi,pi,L);

% figure(1)
% clf
% plot(w6,abs(hf_scaled6))

%% 5.2c)
% Plot each of the BPFs on the same plot.
figure(1)
clf
plot(w2((end+1)/2:end),abs(hf_scaled2((end+1)/2:end)),'-x')
hold on
plot(w3((end+1)/2:end),abs(hf_scaled3((end+1)/2:end)),'-x')
plot(w4((end+1)/2:end),abs(hf_scaled4((end+1)/2:end)),'-x')
plot(w5((end+1)/2:end),abs(hf_scaled5((end+1)/2:end)),'-x')
plot(w6((end+1)/2:end),abs(hf_scaled6((end+1)/2:end)),'-x')
title('Frequency Response of BPFs')
xlabel('\omega (radians)')
ylabel('Normalized Magnitude')
xline(center_freq/fs*2*pi,'--')
yline(0.5,'--b')
legend('Octave 2','Octave 3','Octave 4','Octave 5','Octave 6','Center Frequencies')
hold off

% The plot generated here shows the frequency responses of each BPF, the
% center frequencies of each octave, and the passband level used to
% determine the width of each filter by modifying L.

%% 5.2d)
% In determining the width of each BPF by modifying L, there were some
% sacrifices that had to be made and some assumptions used. In the case of
% the sacrifices, there were moments when determining L that changing
% L by 1 would greatly modify the passband and cause it to either shrink or
% grow greatly on the upper and lower frequency side. In all of these
% cases, the L was chosen that made the passband more narrow. One of the
% assumptions that was used heavily was that the notes being passed into
% the system would exclusively be those within the octaves specified. This
% assumption was used for the upper and lower BPF octaves, as the passband
% could cover frequencies outside of the desired frequencies in the lower
% or upper direction for octaves 2 and 6, respectively.

%% 5.3a)
% Generate a signal xx with the following characteristics:
% t = 0   --> 0.25, cos(2pi(220)t)
% t = 0.3 --> 0.55, cos(2pi(880)t)
% t = 0.6 --> 0.85, cos(2pi(440)t) + cos(2pi(1760)t)

% Sampling frequency
fs = 8000;

% Generate each of the time periods. "n" indicates time where there is no
% signal.
t1 = 0:1/fs:0.25-1/fs;
t2n = 0.25:1/fs:0.3-1/fs;
t3 = 0.3:1/fs:0.55-1/fs;
t4n = 0.55:1/fs:0.6-1/fs;
t5 = 0.6:1/fs:0.85-1/fs;

% Generate each of the individual time period signals.
x1 = cos(2*pi*220*t1);
x2n = 0*t2n;
x3 = cos(2*pi*880*t3);
x4n = 0*t4n;
x5 = cos(2*pi*440*t5)+cos(2*pi*1760*t5);

% Concatenate and generate time vector for all periods
xx = [x1 x2n x3 x4n x5];
t = linspace(0,0.85-1/fs,.85*fs);


%% 5.3b)
% Filter the signal through each of the five filters.
yy2 = conv(xx,hf_scaled2);
yy3 = conv(xx,hf_scaled3);
yy4 = conv(xx,hf_scaled4);
yy5 = conv(xx,hf_scaled5);
yy6 = conv(xx,hf_scaled6);


%% 5.3c)
% Plot each of the five outputs in a single figure
figure(2)
clf
subplot(5,1,1)
plot(t,abs(yy2))
subplot(5,1,2)
plot(t,abs(yy3))
subplot(5,1,3)
plot(t,abs(yy4))
subplot(5,1,4)
plot(t,abs(yy5))
subplot(5,1,5)
plot(t,abs(yy6))

%% 5.3d)
% Use the frequency responses to validate that the correct magnitude and
% phase are outputted in each of the three regions.

% 

%% 5.3e)

%% 5.4a)

%% 5.4b)

%% 5.4c)

%% 5.4d)
