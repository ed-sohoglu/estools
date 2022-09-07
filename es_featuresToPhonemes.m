function phonemes = es_featuresToPhonemes(phonTable,features,type,whichFeatures)

feat2phMatrix = cell2mat(phonTable(2:end,whichFeatures)); % Excluding headers

[~,col2read] = intersect(phonTable(1,:),type);

for t=1:size(features,2)
    [~,row2read] = intersect(feat2phMatrix,features(:,t)','rows');
    if ~isempty(row2read)
        ph = phonTable(1+row2read,col2read);
        phonemes(:,t) = ph{1};
    else
        phonemes(:,t) = 0; 
    end
end