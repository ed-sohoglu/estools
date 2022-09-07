function [M,dwt,dwf] = es_MPS(S,to,fo,ampsrate,plot)
% Adapted by ES (2017) from Theunissen lab ModFilter code.
% Computes and plots modulation spectrum (M) from spectrogram (S). Also outputs temporal and spectral modulation axes (dwt and dwf)
% NB- ampsrate is the sampling rate of the amplitude envelope

% parameters for the modulation spectrum
DBNOISE = 32; % dB in Noise for the log compression - values below will be set to zero.
cmap = 'hot';

% work out size of time-frequency dimensions
nb = size(S,1); % Frequency
nt = length(to); % Time

% work out dB range of input spectrogram
maxB = max(max(S));
minB = maxB-DBNOISE;

% threshold spectrogram
%S(find(S <= minB)) = minB;

% subtract DC from spectrogram
S = S - mean(mean(S));

% calculate the 2D fft (obtains the modulation spectrum)
M = fft2(S);
M = fftshift(M);

% calculate amplitude of modulation spectrum
M = 20*log10(abs(M));

% calculate the labels for spectral frequencies in physical units
fstep = fo(2)-fo(1); % f_step is the separation between frequency bands
for ib=1:ceil((nb+1)/2)
    dwf(ib)= (ib-1)*(1/(fstep*nb));
    if (ib > 1)
        dwf(nb-ib+2)=-dwf(ib);
    end
end
dwf = fftshift(dwf).*1000;

% calculate the labels for temporal frequencies in physical units
for it=1:ceil((nt+1)/2)
    dwt(it) = (it-1)*(ampsrate/nt);
    if (it > 1 )
        dwt(nt-it+2) = -dwt(it);
    end
end
dwt = fftshift(dwt);

if plot
    figure();
    max_amp = max(max(M));
    min_amp = max_amp-DBNOISE;
    imagesc(dwt, dwf, M);
    caxis('manual');
    caxis([min_amp max_amp]);
    title('Amplitude Spectrum');
    axis xy;
    axis([-40 40 0 20]);
    colormap(cmap);
end

