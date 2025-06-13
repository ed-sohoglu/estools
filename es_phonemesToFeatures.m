function features = es_phonemesToFeatures(phonTable,str,type,whichFeatures)

ph2featMatrix = phonTable(:,whichFeatures);
nFeatures = size(ph2featMatrix,2);

col2read = strcmp(phonTable(1,:),type);

for t=1:length(str)
    feaVec = cell2mat(ph2featMatrix(strcmp(phonTable(:,col2read),str(t)),:));
    if ~isempty(feaVec)
        features(:,t) = feaVec;
    else
        features(:,t) = zeros(1,nFeatures); 
    end
end
