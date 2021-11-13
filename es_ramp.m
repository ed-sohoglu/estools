function output = es_ramp(input,fs,duration)
% Apply window

[len_y,len_x] = size(input);
if len_y>len_x
    input = input';
    [len_y,len_x] = size(input);
    resized_flag = 1;
else
    resized_flag = 0;
end

inputLength = len_x;
winLength = round(fs*(duration/1000)*2);

if mod(winLength,2) ~= 0
    winLength = winLength + 1;
end
w = hann(winLength)';
w = repmat(w,size(input,1),1);

input(:,1:winLength/2) = input(:,1:winLength/2) .* w(:,1:winLength/2);
input(:,inputLength-(winLength/2)+1:inputLength) = input(:,inputLength-(winLength/2)+1:inputLength) .* w(:,winLength/2+1:winLength);

if resized_flag
    input = input';
end

output = input;
    