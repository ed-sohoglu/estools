function strNew = es_phonConvert(lookupTableFullFileName,str,typeOld,typeNew)
% SAM-PA	CELEX	CPA	DISC
 
[~, ~, phonTable] = xlsread(lookupTableFullFileName);
[~,col2read] = intersect(phonTable(1,:),typeOld);
phonOld = phonTable(:,col2read);
[~,col2read] = intersect(phonTable(1,:),typeNew);
phonNew = phonTable(:,col2read);

strNew = [];
for seg=1:length(str)
    [~,ind] = intersect(phonOld,str(seg));
    strNew = [strNew phonNew{ind}];
end