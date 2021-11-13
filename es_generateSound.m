function output = generateSound (fs,duration,env,freq)
% generates tone with specified 'duration' (ms).
% 'fs' = sampling frequency (Hz), 'env' argument (if supplied) determines
% rise-fall time (ms) of envelope. Peak amplitude of sound is at full scale (i.e. 1)

inputLength = round(fs*(duration/1000));

% generate tone
t = [1:inputLength]';
sound = sin (2.*pi.*(freq/fs).*t);

if env > 0
    % apply ramps of duration 'env'
    envelopeLength = round(fs*(env/1000)) * 2;

    if envelopeLength > inputLength
        display('Error: ramps are longer than input sound. Sound wont be processed...');
        output = sound;
    else
        w = hann(envelopeLength);

        padLength = inputLength - envelopeLength;
        w = [w(1:envelopeLength/2); ones(padLength,1); w(envelopeLength/2+1:envelopeLength)];

        output = sound .* w;
    end
else
    % don't apply ramps
    output = sound;
end