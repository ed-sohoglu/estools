function output = es_scale(input,lower,upper)
% Scales input to desired range
% e.g. y = scale([1 2 3],0,1). y will be [0 0.5 1]
% Ed Sohoglu 2013

mini = min(input(:));
maxi = max(input(:));

output = input * ((upper-lower)/(maxi-mini));

mini_new = min(output(:));

output = output + (lower-mini_new);
