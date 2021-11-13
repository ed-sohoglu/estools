function [S,to,fo,ampsrate] = es_spectrogram(sound_in,samprate,plot)
% Adapted by ES (2017) from Theunissen lab ModFilter code.
% Computes and plots gaussian spectrogram (S) from audio signal (sound_in). Also outputs time and frequency axes (to and fo)

% parameters for the spectrogram
fband = 32;                                % Width of the frequency band – we use 32 Hz for speech and 125 Hz for zebra finch song (standard deviation of the Gaussian filter)
nstd = 6;
twindow = 1000*nstd/(fband*2.0*pi);        % Window length in ms - 6 times the standard dev of the gaussian window
winLength = fix(twindow*samprate/1000.0);  % Window length in number of points
winLength = fix(winLength/2)*2;            % Enforce even window length
increment = fix(0.01*samprate);            % Sampling rate of spectrogram in number of points - set at 100 Hz
ampsrate = samprate/increment;             % Sampling rate of the amplitude envelope-->Needed for MPS
f_low=70;                                  % Lower frequency bounds to get average amplitude in spectrogram
f_high=5000;                               % Upper frequency bound to get average amplitude in spectrogram
DBNOISE = 32;                              % dB in Noise for the log compression - values below will be set to zero.
cmap = 'hot';

% Get spectrogram
[S, to, fo, pg] = GaussianSpectrum(sound_in, increment, winLength, samprate);

 % Convert spectral coefficients to dB magnitude
S = 20*log10(abs(S));

% Restrict frequency range
fo2keep = find(fo>=f_low & fo<=f_high);
S = S(fo2keep,:);
fo = fo(fo2keep);

if plot
    figure();
    maxB = max(max(S));
    minB = maxB-DBNOISE;
    imagesc(to*1000,fo,S); % to is in seconds
    axis xy;
    caxis('manual');
    caxis([minB maxB]);
    colormap(cmap);
end

function [s, to, fo, pg] = GaussianSpectrum(input, increment, winLength, samprate)
%
% Gaussian spectrum
% 	s = GaussianSpectrum(input, increment, winLength, samprate)
% 	Compute the gaussian spectrogram of an input signal with a gaussian
% 	window that is given by winLength. The standard deviation of that
% 	gaussian is 1/6th of winLength.
%	Each time frame is [winLength]-long and
%	starts [increment] samples after previous frame's start.
%	Only zero and the positive frequencies are returned.
%   to and fo are the time and frequency for each bin in s and Hz
%   pg is a rumming rms.

%%%%%%%%%%%%%%%%%%%%%%%
% Massage the input
%%%%%%%%%%%%%%%%%%%%%%%

% Enforce even winLength to have a symmetric window
if rem(winLength, 2) == 1
    winLength = winLength +1;
end

% Make input it into a row vector if it isn't
if size(input, 1) > 1,
	input = input';
end;

% Padd the input with zeros
pinput = zeros(1,length(input)+winLength);
pinput(winLength/2+1:winLength/2+length(input)) = input;
inputLength = length(pinput);

% The number of time points in the spectrogram
frameCount = floor((inputLength-winLength)/increment)+1;

% The window of the fft
fftLen = winLength;


%%%%%%%%%%%%%%%%%%%%%%%%
% Gaussian window 
%%%%%%%%%%%%%%%%%%%%%%%%
nstd = 6;                   % Number of standard deviations in one window.
wx2 = ((1:winLength)-((winLength+1)/2)).^2;
wvar = (winLength/nstd)^2;
ws = exp(-0.5*(wx2./wvar));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize output "s" 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if rem(fftLen, 2)   
    % winLength is odd
    s = zeros((fftLen+1)/2+1, frameCount);
else
    % winLength is even 
    s = zeros(fftLen/2+1, frameCount);
end

pg = zeros(1, frameCount);
for i=1:frameCount
    start = (i-1)*increment + 1;
    last = start + winLength - 1;
    f = zeros(fftLen, 1);
    f(1:winLength) = ws.*pinput(start:last);
    pg(i) = std(f(1:winLength));

    specslice = fft(f);
    if rem(fftLen, 2)   % winLength is odd
        s(:,i) = specslice(1:((fftLen+1)/2+1));
    else
        s(:,i) = specslice(1:(fftLen/2+1));
    end
    %s(:,i) = specslice(1:(fftLen/2+1));
end

% Assign frequency_label
if rem(fftLen, 2)   % winLength is odd
    select = 1:(fftLen+1)/2;
else
    select = 1:fftLen/2+1;
end
fo = (select-1)'*samprate/fftLen;

% assign time_label
to = ((1:size(s,2))-1)'.*(increment/samprate);
return