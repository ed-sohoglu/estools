function features = es_phonemesToFeatures(phonTable,phonemes,type,whichFeatures)

ph2featMatrix = phonTable(:,whichFeatures);

[~,col2read] = intersect(phonTable(1,:),type);

for t=1:size(phonemes,2)
    [~,row2read] = intersect(phonTable(:,col2read),phonemes(t));
    feaVec = cell2mat(ph2featMatrix(row2read,:));
    if ~isempty(feaVec)
        features(:,t) = feaVec;
    else
        features(:,t) = zeros(1,size(ph2featMatrix,2)); 
    end
end
