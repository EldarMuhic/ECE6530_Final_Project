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
lower_freq(2) = lower_freq(5) * 2^(-1); % C3
lower_freq(1) = lower_freq(2) * 2^(-1); % C2

% Calculate B notes of respective octaves to determine upper frequencies.
upper_freq(2) = lower_freq(3) * 2^(-1/12);  % B3
upper_freq(3) = upper_freq(2) * 2;          % B4
upper_freq(4) = upper_freq(3) * 2;          % B5
upper_freq(5) = upper_freq(4) * 2;          % B6
upper_freq(1) = upper_freq(5) * 2^-1;       % B2

disp('Frequencies in Hz')
table(lower_freq,upper_freq)

disp('Frequencies in radians')
table(lower_freq*2*pi,upper_freq*2*pi)

%% 5.1b)
% Compute the center frequencies of each octave (Octaves 2-6)
center_freq = zeros(5,1);

for i = 1:5
    center_freq(i) = upper_freq(i)-lower_freq(i);
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

% Design of BPF for octave 2.
L = 30;                         % THIS L VALUE IS NOT CORRECT YET
n = linspace(0,L-1,L);

% Calculate impulse response
h = (0.54-0.46.*cos(2.*pi.*n./(L-1))).*cos(center_freq(1)./fs.*pi.* ...
    (n-(L-1)./2));

% Calculate frequency impulse response
hf = fft(h,L);

% Perform fftshift
hf_unscaled = fftshift(hf);

% Scale using calculated beta
beta = 1/max(abs(hf_unscaled));
hf_scaled = beta*hf_unscaled;

% Create frequency vector
w = linspace(-pi,pi,L);

figure(1)
clf
plot(w,abs(hf_scaled))

% For L = 35, passband stretches from to 979.4013 Hz

%% 5.2c)

%% 5.2d)

%% 5.3a)

%% 5.3b)

%% 5.3c)

%% 5.3d)

%% 5.3e)

%% 5.4a)

%% 5.4b)

%% 5.4c)

%% 5.4d)
