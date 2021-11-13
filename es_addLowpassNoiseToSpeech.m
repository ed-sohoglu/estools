function stim = es_addNoiseToSpeech(speech,fs,SNR)

% Add lowpass filtered noise to speech signal according to specified SNR
% speech = audio waveform of speech
% fs = sampling rate (Hz)
% SNR = Signal-To-Noise-Ratio (dB)
% Ed Sohoglu 2021

% work out length of speech signal
nSamples = size(speech,1);

% make noise
ns = sign(rand(nSamples,1)-0.5);
% lowpass filter at 4800 Hz
[blo,alo] = butter(2, 4800/fs);
ns = filter(blo,alo,ns);

ns = ns*(rms(speech)/rms(ns)); % match rms of noise to rms of speech
ns = ns/db2mag(SNR); % apply SNR weighting

stim = speech+ns; % mix speech and noise

% snr(speech,ns); % check SNR