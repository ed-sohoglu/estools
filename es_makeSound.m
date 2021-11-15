function output = es_makeSound(fs,duration,freq,env,type,AMrate,AMdepth)
% makes sound of 'type' tone or noise with specified 'duration' (ms) and carrier frequency (Hz).
% 'fs' = sampling frequency (Hz), 'env' determines rise-fall time (ms) of envelope.
% if 'type' = noise, 'bw' determines its bandwidth (in Hz)
% Peak amplitude of sound is at full scale (i.e. 1)

rand('state',sum(100 * clock));  % initialize random seed

if nargin < 4
    AMdepth = 0;
    AMrate = 0;
    type = 'tone';
    env = 0;
elseif nargin < 5
    AMdepth = 0;
    AMrate = 0;
    type = 'tone';
elseif nargin < 6
    AMdepth = 0;
    AMrate = 0;
end

inputLength = round(fs*(duration/1000));
t = [1:inputLength]';
    
switch type
    
    case 'tone'
        
        sound = sin(2.*pi.*(freq/fs).*t);
        
    case 'noise'
        
        % set variables for filter
        lf = freq/(2^(1/6));  % lowest frequency
        hf = freq*(2^(1/6));  % highest frequency
        lp = lf * duration/1000; % ls point in frequency domain
        hp = hf * duration/1000; % hf point in frequency domain
        N = length(t); % length of noise in samples
        
        sound = randn(N,1);             % Gausian noise
        
        filter = zeros(1, N);           % initializaiton by 0
        filter(1, lp : hp) = 1;         % filter design in real number
        filter(1, N - hp : N - lp) = 1; % filter design in imaginary number
        
        sound = fft(sound);                  % FFT
        sound = sound .* filter';            % filtering
        sound = ifft(sound);                 % inverse FFT
        sound = real(sound);
        
end

if AMrate % modulate carrier with sinewave AM
    sound = es_AM(fs,sound,AMrate,AMdepth);
end

if env
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

output = output / max(abs(output)); % -1 to 1 normalization