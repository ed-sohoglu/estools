function data = es_contrast(data,weights)
% Performs weighted contrast of nsubjects * nconditions data
% Output is a single column vector containing contrast for each subject

tmp = [];
for w=1:size(weights,1)
    tmp(:,w) = sum(data .* repmat(weights(w,:),size(data,1),1),2);
end

data = tmp;

end

