clearvars

%% Choose what to plot

filename = '/imaging/es03/Loes_MEG/meg19_0008/PfmceMfspmeeg_run1_raw_trans1stdef.mat';
modality = 'MEGCOMB';
con2plot = {'Strong+M'};
times2plot = [.09 .1]; % in seconds

%% Load SPM file

D = spm_eeg_load(filename);

chanind = D.selectchannels(modality);
chanind_bad = D.badchannels;
if ~isempty(chanind_bad)
    chanind = setdiff(chanind,chanind_bad);
end
chanlabels = D.chanlabels(chanind);

%% Convert to fieldtrip

data4topo = D.fttimelock(chanind);
data4topo.avg = double(squeeze(data4topo.trial(D.indtrial(con2plot),:,:)));
data4topo.dimord = 'chan_time';
data4topo = rmfield(data4topo,'trial');
if strcmp(modality,'MEGCOMB') % Hack for plotting MEGCOMB
    % Swap MEGCOMB labels for MEG labels (so layout will be based on MEG sensors but data will be from MEGCOMB)
    data4topo.label = D.chanlabels(D.selectchannels('MEG')); 
end

%% Set topoplot config

cfg = [];
cfg.zlim = 'maxabs';
cfg.colormap = 'jet';
cfg.style = 'straight';
cfg.marker = 'off';
cfg.comment = 'xlim';

cfg.channel = data4topo.label;
cfg.xlim  = times2plot;
cfg.parameter = 'avg';
if strcmp(modality,'EEG')
    cfg.elec    = D.sensors('EEG');
else
    cfg.grad   = D.sensors('MEG');
end
cfg.layout = ft_prepare_layout(cfg);

%% Plot

figure;
ft_topoplotER(cfg,data4topo);
colorbar;