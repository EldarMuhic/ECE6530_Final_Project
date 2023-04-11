%% Octave Band Filtering: Lab P-14: 5 Lab Excercise

clear
close all
clc
%% 5.1a)
% Here we start at A4 which is a frequency of 440 Hz and we determine the
% upper and lower frequencies of each octave by using the formula of
% 2^(N/12) where N is the number of notes between the first note and the
% second note. A negative N indicates moving down in pitch.

% Initialize A4 frequency
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

% Create Octave numbering to display in tables
octaves = (2:6)';

disp('Octave Frequencies (Hz)')
table(octaves,lower_freq,upper_freq)

disp('Octave Frequencies (Radians)')
table(octaves,lower_freq*2*pi,upper_freq*2*pi)

%% 5.1b)
% Compute the center frequencies of each octave (Octaves 2-6)
center_freq = zeros(5,1);

for i = 1:5
    center_freq(i) = (upper_freq(i)+lower_freq(i))/2;
end

disp('Octave Center Frequencies (Hz)')
table(octaves,center_freq)

%% 5.2a)
% We must determine a value 'beta' to normalize the frequency response so
% that the max value is equal to 1.
% Our strategy for determining beta will be performed by using MATLAB's max
% function to measure the peak value of the unscaled frequency response,
% and then have MATLAB compute the value that beta must be to normalize the
% peak. The code will look similar to the following.

% beta = 1/max(abs(hf_unscaled));
% hf_scaled = beta*hf_unscaled;

%% 5.2b)
% We must calculate w_c to be used in generating our BPFs. This is done by
% using the sampling frequency and the center frequencies (Hz) found
% previously.

% Sampling frequency (Hz)
fs = 8000;

% Calculate w_c values from center_freq
w_c = center_freq./fs.*2;

% As we design each individual BPF, we must determine the length (L)
% we need to use in order for the BPF to cover the entire bandwidth. This
% is done by using the frequency charts that we created in 5.1 (a) and (b).
% The length of each filter was determined by trial and error and measuring
% the passband (all frequencies with a magnitude above 0.5) after plugging
% in an experimental L until the correct passband was achieved.

%% Design of BPF for octave 2.
L2 = 197;                         % Bandwidth = 134.33Hz-55.45Hz
[h2,hf2] = create_BPF(L2,w_c(1));

% Create frequency vector for plot in 5.2c)
w2 = linspace(-pi,pi,L2);

%% Design of BPF for octave 3.
L3 = 125;                         % Bandwidth = 246.7538Hz-132.5442Hz
[h3,hf3] = create_BPF(L3,w_c(2));

% Create frequency vector for plot in 5.2c)
w3 = linspace(-pi,pi,L3);

%% Design of BPF for octave 4. 
L4 = 67;                         % Bandwidth = 500.89Hz-275.2158Hz
[h4,hf4] = create_BPF(L4,w_c(3));

% Create frequency vector for plot in 5.2c)
w4 = linspace(-pi,pi,L4);

%% Design of BPF for octave 5.
L5 = 35;                         % Bandwidth = 1018Hz-551.3127Hz
[h5,hf5] = create_BPF(L5,w_c(4));

% Create frequency vector for plot in 5.2c)
w5 = linspace(-pi,pi,L5);

%% Design of BPF for octave 6.
% For this filter, we didn't mind that the upper part of the passband
% surpassed the highest frequency in the octave as there won't be any notes
% inputted into the system that are higher than the 6th octave.

L6 = 19;                         % Bandwidth = 2065.2Hz-1109Hz
[h6,hf6] = create_BPF(L6,w_c(5));

% Create frequency vector for plot in 5.2c)
w6 = linspace(-pi,pi,L6);

%% 5.2c)
% Plot each of the BPFs on the same plot.
figure(1)
clf
plot(w2((end+1)/2:end),abs(hf2((end+1)/2:end)),'-x')
hold on
plot(w3((end+1)/2:end),abs(hf3((end+1)/2:end)),'-x')
plot(w4((end+1)/2:end),abs(hf4((end+1)/2:end)),'-x')
plot(w5((end+1)/2:end),abs(hf5((end+1)/2:end)),'-x')
plot(w6((end+1)/2:end),abs(hf6((end+1)/2:end)),'-x')
title('Frequency Response of BPFs')
xlabel('\omega (radians)')
ylabel('Normalized Magnitude')
xline(center_freq/fs*2*pi,'--')
yline(0.5,'--b')
legend('Octave 2','Octave 3','Octave 4','Octave 5','Octave 6','Center Frequencies')
hold off

% The plot generated here shows the frequency responses of each BPF, the
% center frequencies of each octave, and the passband level of 0.5 used to
% determine the width of each filter by modifying L.

%% 5.2d)
% In determining the width of each BPF by modifying L, there were some
% sacrifices that had to be made and some assumptions used. In the case of
% the sacrifices, there were moments when determining L that changing
% L by 1 would greatly modify the passband and cause it to either shrink or
% grow greatly on the upper and/or lower frequency side. In all of these
% cases, the L was chosen that made the passband more narrow. One of the
% assumptions that was used heavily was that the notes being passed into
% the system would exclusively be those within the octaves specified. This
% assumption was used for the upper and lower octave BPFs, as the passband
% could cover frequencies outside of the octaves in the lower
% or upper direction for octaves 2 and 6, respectively (it is assumed that 
% notes won't be in octave 1 or 7).

%% 5.3a)
% Generate a signal xx with the following characteristics:
% t = 0   --> 0.25, cos(2pi(220)t)
% t = 0.3 --> 0.55, cos(2pi(880)t)
% t = 0.6 --> 0.85, cos(2pi(440)t) + cos(2pi(1760)t)

% Sampling frequency
fs = 8000;

% Generate each of the time periods. "n" such as "t2n" or "t4n" indicate
% time where there is no signal.
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
yy2 = conv(xx,h2);
yy3 = conv(xx,h3);
yy4 = conv(xx,h4);
yy5 = conv(xx,h5);
yy6 = conv(xx,h6);


%% 5.3c)
% Generate time vectors for each filtered signal. Size of time vector will
% be length(xx) + length(h_scaledN) - 1.
tt2 = linspace((length(xx)-(length(xx) + L2 - 1))*1/fs,0.85-1/fs,length(xx) + L2 - 1);
tt3 = linspace((length(xx)-(length(xx) + L3 - 1))*1/fs,0.85-1/fs,length(xx) + L3 - 1);
tt4 = linspace((length(xx)-(length(xx) + L4 - 1))*1/fs,0.85-1/fs,length(xx) + L4 - 1);
tt5 = linspace((length(xx)-(length(xx) + L5 - 1))*1/fs,0.85-1/fs,length(xx) + L5 - 1);
tt6 = linspace((length(xx)-(length(xx) + L6 - 1))*1/fs,0.85-1/fs,length(xx) + L6 - 1);

% Plot each of the five outputs in a tiled single figure using subplot.
figure(2)
clf
title('Filtered Signal')

% Octave 2
subplot(5,1,1)
plot(tt2,yy2)
subtitle('BPF Octave #2')
xlabel('Time (s)')
ylabel('Amplitude (V)')

% Octave 3
subplot(5,1,2)
plot(tt3,yy3)
title('BPF Octave #3')
xlabel('Time (s)')
ylabel('Amplitude (V)')

% Octave 4
subplot(5,1,3)
plot(tt4,yy4)
title('BPF Octave #4')
xlabel('Time (s)')
ylabel('Amplitude (V)')

% Octave 5
subplot(5,1,4)
plot(tt5,yy5)
title('BPF Octave #5')
xlabel('Time (s)')
ylabel('Amplitude (V)')

% Octave 6
subplot(5,1,5)
plot(tt6,yy6)
title('BPF Octave #6')
xlabel('Time (s)')
ylabel('Amplitude (V)')

%% 5.3d)
% Examining the filtered signal and comparing it to the passband of each
% BPF, we can see that the filters are behaving as expected. For the first
% segment of the signal, the frequency is 220Hz. This is within octave #3
% which is correctly shown in the plot. For the second segment of the 
% signal, the frequency is 880Hz. This is within octave #5 which is shown
% correctly in the plot. For the third segment of the signal, there were 
% two frequencies present: 440Hz and 1760Hz. These two frequencies are 
% within octaves #4 and #6 which are shown correctly in the plot.

%% 5.3e)
% There is an observable transient at the start-up and end of each filter.
% The time duration of this transient is different for each of the BPFs.
% The measureable time is:
% Octave 2 BPF: 0.0255 seconds
% Octave 3 BPF: 0.0157 seconds
% Octave 4 BPF: 0.0065 seconds
% Octave 5 BPF: 0.0044 seconds
% Octave 6 BPF: 0.0018 seconds

%% 5.4a)
% We're looking to create function where the input signal is passed through
% a BPF and there is a score given which indicates the max 
% amplitude of the output of the BPF. The function for 
% octavescore is defined below.

%% 5.4b)
% The scoring system of the function octavescore gives the maximum output 
% of the filter over each 50 ms interval of the signal.

%% 5.4c)
% Octavescore function was adjusted so that the output signal
% accounts for the delay of the start-up of the BPF (delay is equal to
% (L-1)/2.

%% 5.4d)
% Function autodetect was created. This function indicates whether there 
% are notes corresponding to an octave within each time segment of 50 ms.

% NOTE: Scores of each octave should be combined into a 5 x N matrix before
% plugging into autodetect function.

% Testing our octavescore function and autodetect function on previously
% generated signal.
sc2 = octavescore(xx,h2,fs);
sc3 = octavescore(xx,h3,fs);
sc4 = octavescore(xx,h4,fs);
sc5 = octavescore(xx,h5,fs);
sc6 = octavescore(xx,h6,fs);

sc = [sc2;sc3;sc4;sc5;sc6];

notes = autodetect(sc);

% CODE WORKS FOR SIGNAL GENERATED PREVIOUSLY

%% Test code using incrementing notes (every other note).

% Total time duration of each note.
ttest = 0:1/fs:0.2-1/fs;

% Create empty vector for test notes' frequencies
test_notes = zeros(1,35);

% Calculate frequencies of each note (35 notes total)
for i = 0:34
    test_notes(i+1) = lower_freq(1)*2^(2*i/12);
end

% Create empty vector for test_signal
test_signal = zeros(1,length(ttest)*35);

% Create first note within test_signal
test_signal(1:length(ttest)) = cos(2*pi*test_notes(1)*ttest);

% Define remaining test signal with incrementing notes of 200 ms duration.
for i = 2:35
    test_signal(i*length(ttest)+1:(i+1)*length(ttest)) = cos(2*pi*test_notes(i)*ttest);
end

% Plug test_signal into each BPF using octavescore function
sc2test = octavescore(test_signal,h2,fs);
sc3test = octavescore(test_signal,h3,fs);
sc4test = octavescore(test_signal,h4,fs);
sc5test = octavescore(test_signal,h5,fs);
sc6test = octavescore(test_signal,h6,fs);

% Create 5 x N vector of octavescores
sctest = [sc2test;sc3test;sc4test;sc5test;sc6test];

% Plug 5 x N vector of octavescores into autodetect function
notes_test = autodetect(sctest);

% Test works pretty well. As notes land in the frequencies between BPF
% passbands, there is some error, but otherwise these functions work quite
% well.

%% 5.4.1 Testing

% Load labtest .mat file
labtest = load("labtest.mat");

% Extract signal xx from labtest .mat file
xxtest = labtest.xx;

% Plug signal xx into each BPF using octavescore function
sc2t = octavescore(xxtest,h2,fs);
sc3t = octavescore(xxtest,h3,fs);
sc4t = octavescore(xxtest,h4,fs);
sc5t = octavescore(xxtest,h5,fs);
sc6t = octavescore(xxtest,h6,fs);

% Concatenate octave scores into 5 x N.
sct = [sc2t;sc3t;sc4t;sc5t;sc6t];

% Plug 5 x N octave scores into autodetect function.
notes_t = autodetect(sct);

% Extract actual notes in labtest using labtest.notes
notes = labtest.notes;

% Calculate octaves of labtest.notes

% Initialize empty vectors to decrease computing power.
notes_freq = zeros(2,15);
notes_octave = zeros(2,15);

% For loop to calculate frequency of notes and octave.
for i = 1:15
    for j = 1:2
        if notes(j,i) ~= 0
            notes_freq(j,i) = 440*2.^((notes(j,i)-49)./12);
        else
            notes_freq(j,i) = 0;
        end

        for k = 1:5
            if notes_freq(j,i) <= upper_freq(k) && notes_freq(j,i) >= lower_freq(k)
                notes_octave(j,i) = k+1;
            end
        end
    end
end

% After examining the performance of the BPFs against the solution for the
% labtest.xx signal, there were 29 incorrect octave designations out of 375
% total, and so an error of 7.73%. It appeared that the notes that were the
% most troublesome were those that were on the boundaries between octaves
% as well as notes in the 2-4 octave range. We believe that because the
% transients in these ranges are longer, they affected the max amplitudes
% calculated by the octavescore function, and therefore the rest of the
% output. When two simultaneous notes were played in the same octave, the
% octaves next to that octave were shown as having notes being played in
% them as well. We believe this is due to the magnitude of the frequencies
% in that octave essentially being doubled which caused runoff into the
% surrounding octaves.

% Overall, we're very happy with how the BPFs performed with the labtest.xx
% data and that the error percentage was relatively low.

%% Functions

function [h,hf_scaled] = create_BPF(L,w_c)
%CREATE_BPF
% usage: [h,hf_scaled] = create_BPF(L, w_c)
% returns the impulse response of a BPF of input length L and center
% frequency w_c in both the frequency (hf_scaled) and time (h) domain.
% L = Length of BPF. This controls the width of the passband.
% w_c = center frequency of BPF. Calculated by
% center_freq(Hz)/sampling_freq/2.
% The output is the normalized impulse response of the BPF in the time
% domain.

% Create sample vector n from input L
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
hf_scaled = beta*hf_unscaled;           % hf_scaled

% Transform back to time vector
hf_scaled = ifftshift(hf_scaled);
h = ifft(hf_scaled);                    % h

end

function sc = octavescore(xx, hh, fs)
%OCTAVESCORE
% usage: sc = octavescore(xx, hh, fs)
% returns a score based on the max amplitude of the filtered output
% xx = input signal containing musical notes
% hh = impulse response of ONE bandpass filter
% fs = sampling rate
%
% The signal detection is done by filtering xx with a length-L
% BPF, hh, and then finding the maximum amplitude of the output
% within 50 millisecond segments.
% The score is a vector containing the maximum amplitudes
% of all the segments.

% Length of a segment based on sampling rate
seg = 0.05*fs;

% Length of filter
L = length(hh);

% Calculate number of segments
i = length(xx)/seg;
i = floor(i);           % Make sure i is an integer

% Initialize sc vector based on number of segments
sc = zeros(1,i);

% For loop moves through 
for j = 1:i

    xxf = fft(xx((j-1)*seg+1:j*seg),L);
    hhf = fft(hh,L);
    yyf = xxf.*hhf;
    yy = ifft(yyf,L);

    % Adjust yy for transient at beginning due to Hamming BPF
    % start-up delay.
    yy = yy((L-1)/2:end);
    sc(j) = max(yy);

end

end

function notes = autodetect(sc)
%AUTODETECT
% sc = 5 x N matrix with scores of all 5 octaves. N is equal to the number
% of 50 ms segments.
% This function will take the 5 x N matrix of scores and give a logical 1
% or zero whether there are notes present in an octave within that 50 ms
% segment.

notes = zeros(5,length(sc));

for i = 1:5
    for j = 1:length(sc)
        if abs(sc(i,j)) >= 0.5
            notes(i,j) = 1;
        else 
            notes(i,j) = 0;
        end

    end
end

end