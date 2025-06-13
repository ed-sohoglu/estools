function strNew = es_phonConvert(phonTable,str,typeOld,typeNew)
% Takes cell array of strings and does look-up operation in phonTable to
% convert between typeOld and typeNew transcription systems
% Cell array input can comprise multiple cells (for multiple tments
% comprising a word)
% Helpful to specify input as cell array because of SAMPA which encodes
% some tments with multiple characters e.g. U@
% Ed Sohoglu 2025
 
phonOld = phonTable(:,strcmp(phonTable(1,:),typeOld));
phonNew = phonTable(:,strcmp(phonTable(1,:),typeNew));

strNew = {};
for t=1:length(str)
    strNew{t} = phonNew{strcmp(phonOld,str(t))};
end