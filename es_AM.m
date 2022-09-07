function output = es_AM(fs,input,rate,depth)

len = length(input);
t = [1:len]';
modulator = 1 + sin(2.*pi.*(rate/fs).*t) * depth;
%modulator = modulator * 0.5;
output = input .* modulator;
