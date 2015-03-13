% % ---------------------------------------------
% Name : Manan Vyas
% USCID: 7483-8632-00
% Email: mvyas@usc.edu
% EE519 : Speech Recognition : P5
% % ----------------------------------------------

% Setup
clc;
clear all;
close all;

load('final2014_p5.mat');
N = length(speech1_10k); % Length of speech
Fs = 10000; % Hz  Sampling Frequency

% Part (A) - Windowing
% Ambiguity in language here : It is said that it should be assumed that
% the signal center starts from n=0 ; i.e: at n=0 then the sample at n=0 
% being the 500th sample of the original sample and the window starts from
% this sample to the next 250 samples.

% The language of the question can also be interpreted and the center neing
% the 500th sample and the window extending from sample 375 to 625th
% sample.However, I have gone ahead with the above interpretation.
% Although, it doesnt make any difference as the input signal is periodic
% and both the chunks are similiar.
% Please discuss it with me if any issues/clarifications regards to this !

% Setup the Hamming window => 25msec => 250 samples;
wLen = 250;
window = hamming(wLen);
% Window the signal - considering the instant n=500 as the new center of
% the sppech signal i.e: n' = 0 (500) and hence applying the window from
% sample instants 500 to 750 for windowing
Speech = speech1_10k(500:749);
wSpeech = Speech.*window'; % windowed speech here
% Log magnitude plot of the DFT of the windowed signal
figure
hold on
logMag = (20*log(abs(fft(wSpeech,1024)))); % FFT -> log mag(FFT)
Mag = abs(fft(wSpeech,1024));
phase = angle(fft(wSpeech,1024));
Mag = Mag(1:512);
logMag = logMag(1:512);
phase = phase(1:512);
subplot(2,1,1); axis tight;
plot(wSpeech); title('Hamming Windowed Speech Signal'); xlabel('n ->');ylabel('Amplitude');
subplot(2,1,2)
plot(logMag);
xlabel('FFT Points'); ylabel('Magnitude (dB)'); axis tight;
title('Log Magnitude Spectrum of the Windowed signal over 512 FFT points');

% Part (B) - Peak Picking on the 512 points of the Log Magnitude spectrum
pickedPeaks = final_2014_peak_pick(Mag);
% We have the location of the picked peaks in pickedPeaks here
% Now finding the corresponding log mag values at those peaks the the
% frequency location corresponding to those FFT locations
for i=1:1:length(pickedPeaks)
   magAtPeaks(i) = Mag(pickedPeaks(i)); 
   phaseAtPeaks(i) = phase(pickedPeaks(i));
   frequencyAtPeaks(i) = (pickedPeaks(i));
end
L = length(pickedPeaks);
% Generate a 3D array containing this information
peakData(1,length(pickedPeaks),1) = Mag(length(pickedPeaks));
peakData(1,length(pickedPeaks),2) = phaseAtPeaks(length(pickedPeaks));
peakData(1,length(pickedPeaks),3) = frequencyAtPeaks(length(pickedPeaks));
frequencyAtPeaks = (frequencyAtPeaks).*(pi/1024);
% Got peak Data
% Plot the superimposed peaks on the original peaks
figure
hold on
plot(logMag);
scale = 1:(512):(512*60);
plot((frequencyAtPeaks).*(1025/pi),20*log(magAtPeaks),'rx');axis tight;
xlabel('FFT points'); ylabel('20*log(magnitude(fft))');
title('Log magitude FFT and magnitude of picked peaks superimposed');
legend('20*log|FFT(windowedSpeech|','log magnitude of peaks picked');
hold off;

% Part (C) : Reconstruction of the signal using Sinusiods and above data
% Now reconstruct using the original signal using sum of sines
Sn = generateSine(magAtPeaks, phaseAtPeaks, frequencyAtPeaks, length(speech1_10k));
figure
hold on;
axis tight;
subplot(2,1,1)
plot(Sn); xlabel('n ->'); ylabel('Amplitude'); title('Reconstructed speech signal');
subplot(2,1,2)
plot(speech1_10k,'r-'); xlabel('n ->'); ylabel('Amplitude'); title('Original speech signal');

% Part (D) : Reconstruction of the signal using Sinusiods and above data
% with phase offset 0
% Now reconstruct using the original signal using sum of sines; phase -> 0
phaseAtPeaks = zeros(1,length(pickedPeaks));
Sn2 = generateSine(magAtPeaks, phaseAtPeaks, frequencyAtPeaks, length(speech1_10k));
figure
hold on;
axis tight;
subplot(2,1,1)
plot(Sn./(max(Sn))); xlabel('n ->'); ylabel('Amplitude'); title('Reconstructed speech signal with phase offset 0');
subplot(2,1,2)
plot(speech1_10k,'r-'); xlabel('n ->'); ylabel('Amplitude'); title('Original speech signal');

% Part (E) : Comparision between reconstructed signals obtained from
% Part(C) and Part(D)
figure
hold on;
plot(Sn./(max(Sn)));
plot(Sn2./(max(Sn2)),'r-'); xlabel('n ->'); ylabel('Amplitude'); title('Comparision between part(C) and part(D)');
legend('Reconstructed signal','Reconstructed signal with phase offset = 0')



