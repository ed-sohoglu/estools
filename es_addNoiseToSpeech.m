function [mix,speech,ns] = es_addNoiseToSpeech(speech,SNR)

% Add spectrally matched noise to speech signal according to specified SNR
% speech = audio waveform of speech
% SNR = Signal-To-Noise-Ratio (dB)
% Ed Sohoglu 2021

% NB- ASSUMES MONO!!!!

% work out length of speech signal
nSamples = size(speech,1);

% make noise (spectrally matched to speech)
mag = abs(fft(speech));
ph = rand(nSamples,1)*2*pi-pi;
ns = real(ifft(mag.*exp(1i*ph)));

ns = ns*(rms(speech)/rms(ns)); % match rms of noise to rms of speech
ns = ns/db2mag(SNR); % apply SNR weighting

mix = speech+ns; % mix speech and noise

%snr(speech,ns) % check SNR