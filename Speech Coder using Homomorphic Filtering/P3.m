% % ---------------------------------------------
% Name : Manan Vyas
% USCID: 7483-8632-00
% Email: mvyas@usc.edu
% EE519 : Speech Recognition : P3
% % ----------------------------------------------

% Setup
clc;
clear all;
close all;


Fs = 10000;   % Hz - Sampling frq
load('final2014_p3.mat'); % we have normalized Speech signal here
n = length(speech);
% Apply a 25ms Hamming window ...
% 10k samples/sec => 25msec = 250 samples => Length of Hamming window
wLen = 250;
window = hamming(wLen);
% Take first 250 samples of the Speech signal
Speech = speech(1:250);
wSpeech = Speech.*window' ; % Windowed Speech

% Computing Cepstrum
wSpeechFFT = fft(wSpeech,1024); % 1024 point FFT of the windowed Speech
wSpeechFFT_log = log(abs(wSpeechFFT)); % Log Mag spectrum of the FFT
wSpeechFFT_cepstrum = real(ifft(wSpeechFFT_log,1024)); % Defintion of Cepstrum

% Part (A) All plots
figure;
subplot(3,1,1); plot(Speech); title('250 samples of original Speech sample');
xlabel('n ->'); ylabel('Speech');

subplot(3,1,2); plot(wSpeech); 
title('Windowed version of above Speech signal using Hamming Window');
xlabel('n ->'); ylabel('s[n].W(hamming)[n]');

subplot(3,1,2); plot(wSpeech); 
title('Windowed version of above Speech signal using Hamming Window');
xlabel('n ->'); ylabel('s[n].W(hamming)[n]');

subplot(3,1,3); plot(wSpeechFFT_cepstrum(1:200)); axis([1 200 -1 1.5]);
title('Real cepstrum');
xlabel('Quefrency ->'); ylabel('Real ceps amplitude');

% Part (B) : Pitch Estimation
%Liftering

wSpeechFFT_cepstrum_coEff = wSpeechFFT_cepstrum(1:length(wSpeechFFT_cepstrum));
L = zeros(1,length(wSpeechFFT_cepstrum_coEff)); 
% Liftering Window
L(20:140) = 1;
% Low time lifted cepstrum
yOp = real(wSpeechFFT_cepstrum_coEff.*L);

%Finding peak in lifted cepstrum
[peak_val,peak_loc] = max(yOp);
pitch_period = peak_loc;
pitch_frq = (1/pitch_period)*Fs;

% Part (C) : Uniform Quantizer
pitchMin = 50;
pitchMax = 300;
B = 7 ; % Bit precision size
% The quantized value of the pitch frequency using a 7bit uniform quantizer
quantPitchFrqValue = quant(pitch_frq,B,pitchMin,pitchMax)
quantValue = quant(pitch_frq,B,pitchMin,pitchMax); 

% Part (D): Non uniform quantizer
% take the initial 28 co-efficients of the cepstrum co-eff
cepstrumCoEff28 = wSpeechFFT_cepstrum(1:28);

% Unifromly quantize these
for i=1:1:28
    cepstrumCoEff28_UniQ(i) = quant(cepstrumCoEff28(i), 7, min(cepstrumCoEff28),max(cepstrumCoEff28));
end
figure
subplot(2,1,1)
plot(cepstrumCoEff28,'--r');
hold on
stem(cepstrumCoEff28_UniQ);
title('Uniformly Quantized initial 28 cepstrum co-effs');
xlabel('n ->'); ylabel('c[n]'); axis([1 28 min(cepstrumCoEff28_UniQ) max(cepstrumCoEff28_UniQ)]);
legend('Original Cepstrum Coeffs','Uniformly Quantized Cepstrum Coeff(7bit)');
hold off

% Non Uniform Quantization
cepstrumCoEff28_NonUniQ = nonUniformQuant(cepstrumCoEff28);
subplot(2,1,2)
plot(cepstrumCoEff28,'--r');
hold on
stem(cepstrumCoEff28_NonUniQ);
title('Non-Uniformly Quantized initial 28 cepstrum co-effs');
xlabel('n ->'); ylabel('_c[n]'); axis([1 28 min(cepstrumCoEff28_NonUniQ) max(cepstrumCoEff28_NonUniQ)]);
legend('Original Cepstrum Coeffs','Non-Uniformly Quantized Cepstrum Coeff');
hold off

% Part (E): Minimum phase reconstruction
% Right sided Cepstral Lifter
Lifter(1) = 1;
Lifter(2:28) = 2;
% Lifter the non uniformly quantized values
cepstrumCoEff28_UniQ_Liftered = cepstrumCoEff28_UniQ.*Lifter;
minimumPhaseRecostruction = fft(cepstrumCoEff28_UniQ_Liftered,1024);

% Part (F) : Sampled log magnitude and phase and signal reconstruction
figure
% We need to map the x-axis value to reflect frequency in Hz. 
% To do that all you need to remember that the length of the
% FFT corresponds to the sampling rate Fs for continuous frequencies and 
% corresponds to 2? for discrete frequency. Therefore
% to get the positive and negative frequencies
scale = -(1024/2):(1024/2-1);
subplot (2,1,1)
plot(scale*Fs/1024,fftshift(20*log10(abs(minimumPhaseRecostruction)))); % to get the frquency scale in Hz
title('Log Magnitude plot of minimum phase recosntruction');xlabel('f(Hz)->');ylabel('20*log|fft(Ceps,1024)| db scale');
subplot(2,1,2);
plot(scale*Fs/1024,fftshift((angle(minimumPhaseRecostruction))));
title('Phase plot of minimum phase recosntruction');xlabel('f(Hz)->');ylabel('<fft(Ceps,1024)');
% Now sample the log magnitude and phase values
toBeSampled = log(abs(minimumPhaseRecostruction));
diff = [5000,-ones(1,100)];
index = 1;
newHarms = 1;
for i=1:1:512 %because the spectrum FFT is symmetrical
    index = index + 1;
    newfreq = i*Fs/1024;
    tempdiff = abs(quantValue - (i*Fs/(1024)));
    diff(index) = quantValue - (i*Fs/(1024));
    if(abs(diff(index)) > abs(diff(index-1)))
%        sampledLogMag(i) = 20*log(abs(minimumPhaseRecostruction((quantValue+diff(index-1))*1024/(Fs*(i)))));
        sampledLogMag(newHarms) = toBeSampled(i) ;%log(abs(minimumPhaseRecostruction(i)));
        sampledPhase(newHarms) = angle(minimumPhaseRecostruction(i));
        selectedHarmonicFrequencies(newHarms) = newfreq;
        selectedHarmonicFrequenciesIndex(newHarms) = i;
        quantValue = quantValue + quantPitchFrqValue;
        newHarms = newHarms + 1;
        diff = [5000,-ones(1,100)];
%         i = 1;
        if(i > 1024)
           break; 
        end
        index = 1;
    end
end
selectedHarmonicFrequencies
selectedHarmonicFrequenciesIndex
sampledLogMag
sampledPhase

% Part (G) : Resynthesis
reconLogMag = exp(sampledLogMag);
reconPhase = exp(sampledPhase);
reconstructedSpeech = sinewave(reconLogMag,selectedHarmonicFrequencies,reconPhase,1000);
% Normalize
reconstructedSpeech = reconstructedSpeech./max(reconstructedSpeech) ;
figure
subplot(2,1,1)
plot(reconstructedSpeech); title('Reconstructed Speech'); xlabel('n ->'); ylabel('y[n]');
subplot(2,1,2);
plot(speech);title('Original Speech'); xlabel('n ->'); ylabel('x[n]');
% MSE Calculation
sqDiff = (double(reconstructedSpeech) - double(speech)).^2;
MSE = sum(sqDiff)/(length(reconstructedSpeech))





























