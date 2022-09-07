function output = es_rms(input,dim)
% RMS across data dimension specified by 'dim' argument 

output = sqrt(nanmean(input.^2,dim));