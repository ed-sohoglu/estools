modality = 'MEG';
filename = [];
clusterImageName = [];

D = spm_eeg_load(filename);

V = spm_vol(clusterImageName);
data_img = spm_read_vols(V);
data_img(~data_img) = NaN;
data_img = max(data_img,[],3);

Cind = D.selectchannels(modality);
[Cel, x, y] = spm_eeg_locate_channels(D, V.dim(1), Cind);

goodchan = D.indchantype(modality, 'GOOD');

%-Find channels corresponding to blobs
%--------------------------------------------------------------------------
chanind = [];
for i = 1:length(Cind)
    if ismember(Cind(i), goodchan)
        val = data_img(Cel(i, 1), Cel(i, 2));
        if ~isnan(val) && val~=0
            chanind = [chanind; Cind(i)];
        end
    end
end

D.chanlabels(chanind)
