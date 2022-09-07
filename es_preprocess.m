function es_preprocess(step,prevStep,p,pathstem,subjects,dates,blocksin,blocksout,experimenterID,rawpathstem,badeeg,SID)

switch prevStep
    
    case 'CTF'
        prevStep = '*.meg4'; 
    case 'maxfilter'
        prevStep = '*ssst.fif';    
    case 'convert'   
        prevStep = 'spmeeg*.mat';
    case 'downsample'   
        prevStep = 'd*.mat';
    case 'epoch'     
        prevStep = 'e*.mat';
    case 'merge'    
        prevStep = 'c*.mat';
    case 'rereference'     
        prevStep = 'M*.mat';
    case 'TF_power'
        prevStep = 'tf*.mat';
    case 'TF_phase'
        prevStep = 'tph*.mat';
    case 'TF_rescale'
        prevStep = 'r*.mat';
    case 'filter'     
        prevStep = 'f*.mat';
    case 'baseline'   
        prevStep = 'b*.mat';
    case 'average'
        prevStep = 'm*.mat';
    case 'weight'
        prevStep = 'w*.mat';
    case 'combineplanar'
        prevStep = 'p*.mat';
    case 'adjust_sensor_positions'
        prevStep = 'a*.mat';
    case 'DSS'
        prevStep = 'n*.mat';
    case 'artefact'   
        prevStep = 'a*.mat';
    case 'blink'
        prevStep = 'clean*.mat';
    case 'image'
        prevStep = 'condition*.nii';
    case 'smooth'
        prevStep = 'sm*.img';
    case 'firstlevel'
        prevStep = 't*.img';       
end

switch step
    
    case 'erase'
        
        for s=1:length(subjects)
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);
            
            % search for input files
            files = dir(prevStep);
            
            for f=1:length(files)
                
                fprintf([ '\n\nProcessing ' files(f).name '...\n\n' ]);
                
                % main process
                delete(files(f).name);       
                if strfind(files(f).name,'spmeeg') % if spm8 file, delete .dat file (in addition to .mat file)     
                    datname = [strtok(files(f).name,'.') '.dat'];
                    delete(datname);
                end
                
            end % blocks
            
        end % subjects
        
        fprintf('\n\nData deleted!\n\n');
        
    case 'copy'
        
        % parameters
        outputstem = p.outputstem;
        
        for s=1:length(subjects)
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);
            
            % search for input files
            files = dir(prevStep);
            
            for f=1:length(files)
                
                fprintf([ '\n\nProcessing ' files(f).name '...\n\n' ]);
                
                % make output directory if it doesn't exist
                outputfullpath = [outputstem subjects{s}];
                if ~exist(outputfullpath,'dir')         
                    mkdir(outputfullpath);
                end
                
                % main process
                copyfile(files(f).name,outputfullpath);
                if strfind(files(f).name,'spmeeg') % if spm8 file, copy .dat file (in addition to .mat file)     
                    datname = [strtok(files(f).name,'.') '.dat'];
                    copyfile(datname,outputfullpath);
                end
                
            end % blocks
            
        end % subjects
        
        fprintf('\n\nData copied!\n\n');
        
    case 'definetrials'
        
        % parameters for SPM function
        fs = p.fs_new; % ASSUMES DOWNSAMPLED
        S.timewin = [p.preEpoch p.postEpoch];
        S.reviewtrials = 0;
        S.save = 0;
        
        % other parameters
        triggerChannelName = p.triggerChannelName;
        conditions = p.conditions;
        triggers = p.triggers;
        if isfield(p,'stimuli_list_fname'); stimuli_list_fname = [pathstem p.stimuli_list_fname]; end      
        
        for s=1:length(subjects)
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);

            % search for input files
            files = dir(prevStep);
            
            for f=1:length(files)
                
                fprintf([ '\n\nProcessing ' files(f).name '...\n\n' ]);
                
                % define trigger labels
                for c=1:length(conditions)
                    S.trialdef(c).conditionlabel = conditions{c};
                    S.trialdef(c).eventtype = triggerChannelName;
                    S.trialdef(c).eventvalue = triggers(c);
                end
                
                % set input file
                S.D = files(f).name;

                % main process
                trials = [];
                [trials.trl trials.labels] = spm_eeg_definetrial(S);
                if strcmp(SID{s},'016') && strcmp(blocksout{s}{f},'03')
                    if strfind(pathstem,'scene')
                        trials.trl = trials.trl(2:end,:);
                        trials.labels = trials.labels(2:end);
                    elseif strfind(pathstem,'change')
                        trials.trl = trials.trl(3:end,:);
                        trials.labels = trials.labels(3:end);
                    end
                end
                    
                fprintf('\n%d trigger events found...',length(trials.trl));
                
                % correct for trigger-to-sound delay
                if isfield(p,'delay')
                    trials.trl(:,1:2) = trials.trl(:,1:2) + round((p.delay/1000)*fs);
                end
                
                % incorporate event information for P6E1 (if stimulus list
                % filename supplied)
                if exist('stimuli_list_fname','var')
                    trials = es_get_events_MEG4(SID{s},sprintf('%03d',str2double(blocksout{s}{f})),trials,stimuli_list_fname,pathstem);
                    fprintf('\n%d trigger events remaining after matching triggers to stimuli list...',length(trials.trl));
                    
                    if ~isempty(find(ismember(trials.labels,'CA_REG_REG_Detected')==1))
                        ind = find(ismember(trials.labels,'NC_REG_REG_CorrectRejection')==1);
                        trials.labels(ind) = {'NC(CA)_REG_REG_CorrectRejection'};
                        ind = find(ismember(trials.labels,'NC_RAND_REG_CorrectRejection')==1);
                        trials.labels(ind) = {'NC(CA)_RAND_REG_CorrectRejection'};
                    elseif ~isempty(find(ismember(trials.labels,'CD_REG_REG_Detected')==1))
                        ind = find(ismember(trials.labels,'NC_REG_REG_CorrectRejection')==1);
                        trials.labels(ind) = {'NC(CD)_REG_REG_CorrectRejection'};
                        ind = find(ismember(trials.labels,'NC_RAND_REG_CorrectRejection')==1);
                        trials.labels(ind) = {'NC(CD)_RAND_REG_CorrectRejection'};
                    end
                    if ~isempty(find(ismember(trials.labels,'CA_REG_REG_Detected')==1))
                        ind = find(ismember(trials.labels,'NC_REG_REG_FalseAlarm')==1);
                        trials.labels(ind) = {'NC(CA)_REG_REG_FalseAlarm'};
                        ind = find(ismember(trials.labels,'NC_RAND_REG_FalseAlarm')==1);
                        trials.labels(ind) = {'NC(CA)_RAND_REG_FalseAlarm'};
                    elseif ~isempty(find(ismember(trials.labels,'CD_REG_REG_Detected')==1))
                        ind = find(ismember(trials.labels,'NC_REG_REG_FalseAlarm')==1);
                        trials.labels(ind) = {'NC(CD)_REG_REG_FalseAlarm'};
                        ind = find(ismember(trials.labels,'NC_RAND_REG_FalseAlarm')==1);
                        trials.labels(ind) = {'NC(CD)_RAND_REG_FalseAlarm'};
                    end
                else                   
                    if ~isempty(find(ismember(trials.labels,'CA_REG_REG')==1))
                        ind = find(ismember(trials.labels,'NC_REG_REG')==1);
                        trials.labels(ind) = {'NC(CA)_REG_REG'};
                        ind = find(ismember(trials.labels,'NC_RAND_REG')==1);
                        trials.labels(ind) = {'NC(CA)_RAND_REG'};
                    elseif ~isempty(find(ismember(trials.labels,'CD_REG_REG')==1))
                        ind = find(ismember(trials.labels,'NC_REG_REG')==1);
                        trials.labels(ind) = {'NC(CD)_REG_REG'};
                        ind = find(ismember(trials.labels,'NC_RAND_REG')==1);
                        trials.labels(ind) = {'NC(CD)_RAND_REG'};
                    end
                end
                               
                % save to file
                filename = strtok(files(f).name,'.');
                save(['trlDef_' filename],'trials');
                
            end % blocks
            
        end % subjects
        
        fprintf('\n\nTrials defined!\n\n');
        
    case 'definetrials_jp' % uses jp's trigger extraction function + includes ability to block out response triggers
        
        % parameters for jp's trigger function
        cfg.fs = p.fs; % ASSUMES NOT DOWNSAMPLED
        cfg.triggers = p.triggers;
        cfg.prestim = -p.preEpoch/1000; % jp's trigger function specifies times with positive numbers only
        cfg.poststim = p.postEpoch/1000;
        cfg.minduration = p.minduration/1000;
        cfg.maxduration = p.maxduration/1000;
        cfg.fromprevious = 0;
        
        % other parameters
        conditions = p.conditions;
        if isfield(p,'stimuli_list_fname')
            stimuli_list_fname = [pathstem p.stimuli_list_fname];
        end
        
        for s=1:length(subjects)
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);
            
            % search for input files
            files = dir(prevStep);
            
            for f=1:length(files)
                
                fprintf([ '\n\nProcessing ' files(f).name '...\n\n' ]);
                
                % main process
                
                % only read select triggers channels (i.e. to block out response triggers)
                B = ft_read_data(files(f).name,'chanindx',[307:310]); % trigger channels STI001 STI002 STI003 STI004 (can represent numbers 0-15. if higher numbers needed, add more channels)
                B(B>0) = 1; % values are 0 or 5, make it binary            
                % for each row, get the values encoded by that channel
                for k=1:size(B,1)
                    B(k,B(k,:)>0) = 2^(k-1);
                end      
                % add up these channels to get the new trigger channel
                sti101_new = sum(B,1);
                
                % get trials
                if strfind(files(f).name,'train')
                    cfg.triggers = p.triggers(4:9);
                    conditions = p.conditions(4:9);
                else
                    cfg.triggers = p.triggers(1:3);
                    conditions = p.conditions(1:3);
                end
                [trials.trl events] = jp_meg_gettrials(sti101_new, cfg);
                for c=1:length(conditions)
                    trials.labels([events.value]==cfg.triggers(c),1)=repmat(conditions(c),sum([events.value]==cfg.triggers(c)),1);
                end
                
                % correct for trigger-to-sound delay
                if isfield(p,'delay')
                    trials.trl(:,1:2) = trials.trl(:,1:2) + round((p.delay/1000)*cfg.fs);
                end
                
                % incorporate event information for P6E1 (if stimulus list
                % filename supplied)
                if exist('stimuli_list_fname','var')
                    trials = es_get_events_P6E1(subjects{s},files(f).name,trials,stimuli_list_fname);
                    fprintf('\n%d trigger events remaining after matching triggers to stimuli list...',length(trials.trl));
                end
                
                % save to file
                filename = strtok(files(f).name,'.');
                save(['trlDef_' filename],'trials');
                
            end % blocks
            
        end % subjects
        
        fprintf('\n\nTrials defined!\n\n');
     
      case 'artefact_ft'
        
        % parameters for fieldtrip function
        type = p.type;
        z_thresh = p.z;
        feedback = p.feedback;
        if isfield(p,'artefactchans') % look for artefacts in specified channels
            channels = p.artefactchans;
        elseif isfield(p,'montage_fname') % or else look in all EEG and EOG channels
            load([pathstem p.montage_fname]);
            channels = montage.labelorg;
            ind_eeg = find(cellfun(@isempty,strfind(channels,'EEG'))==0);
            ind_eog = find(cellfun(@isempty,strfind(channels,'EOG'))==0);
            channels = channels([ind_eeg ind_eog]);
        else
            error('Please specify channels to look for artefact in');
        end
        
        for s=1:length(subjects)
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);
            
            % search for input files
            files = dir(prevStep);
            filesTrlDef = dir('*trlDef.mat');

            for f=1:length(files)
                
                fprintf([ '\n\nProcessing ' files(f).name '...\n\n' ]);
                
                % loads trial structure containg trial def
                load(filesTrlDef(f).name);
                
                % exclude manually specified bad EEG/EOG channels
                channels_clean = {};
                if ~isempty(badeeg{s}) 
                    channels_clean = setdiff(channels,badeeg{s});
                end
                
                % load data
                cfg = [];
                cfg.dataset = files(f).name;
                if ~isempty(channels_clean)
                    cfg.channel = channels_clean;
                else
                    cfg.channel = channels;
                end
                cfg.continuous = 'yes';
                data = ft_preprocessing(cfg);
                
                % setup cfg structure
                cfg = [];
                cfg.trl = trials.trl;
                if cfg.trl(end,2) > length(data.time) % remove last trial if defined to be outside time range of data
                    fprintf('\nNot checking last trial since trial offset defined outside of data range\n');
                    cfg.trl(end,:) = [];
                end

                % cutoff and padding parameters
                if ~isempty(channels_clean)
                    cfg.artfctdef.zvalue.channel = channels_clean;
                else
                    cfg.artfctdef.zvalue.channel = channels;
                end
                if iscell(z_thresh)
                    cfg.artfctdef.zvalue.cutoff = z_thresh{s}(f);
                else
                    cfg.artfctdef.zvalue.cutoff = z_thresh;
                end
                cfg.artfctdef.zvalue.trlpadding = 0;
                cfg.artfctdef.zvalue.artpadding = 0;
                cfg.artfctdef.zvalue.fltpadding = 0.1;
                
                % algorithmic parameters
                if strcmp(type,'EOG')
                    cfg.artfctdef.zvalue.bpfilter   = 'yes';
                    cfg.artfctdef.zvalue.bpfilttype = 'but';
                    cfg.artfctdef.zvalue.bpfreq     = [1 15];
                    cfg.artfctdef.zvalue.bpfiltord  = 4;
                    cfg.artfctdef.zvalue.hilbert    = 'yes';
                elseif strcmp(type,'muscle')
                    cfg.artfctdef.zvalue.bpfilter   = 'yes';
                    cfg.artfctdef.zvalue.bpfilttype = 'but';
                    cfg.artfctdef.zvalue.bpfreq     = [90 120];
                    cfg.artfctdef.zvalue.bpfiltord  = 9;
                    cfg.artfctdef.zvalue.hilbert    = 'yes';
                    cfg.artfctdef.zvalue.boxcar  = 0.2;
                end
                
                % feedback
                cfg.artfctdef.zvalue.feedback = feedback;
                
                if strcmp(type,'EOG')
                    % detect artefacts
                    [cfg artefact_eog] = ft_artifact_zvalue(cfg,data);
                    % reject artefacts
                    cfg=[];
                    cfg.dataset = files(f).name;
                    cfg.trl = trials.trl;
                    cfg.artfctdef.reject = 'complete'; % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
                    cfg.artfctdef.eog.artifact = artefact_eog;
                    trials_clean = ft_rejectartifact(cfg);
                    trials.artefact_eog = setdiff(trials.trl(:,1),trials_clean.trl(:,1)); % stores onset of artefact trial
                    fprintf('\n%f percent of trials rejected',(length(trials.artefact_eog)/length(trials.trl))*100);
                elseif strcmp(type,'muscle')
                    % detect artefacts
                    [cfg artefact_muscle] = ft_artifact_zvalue(cfg,data);
                    % reject artefacts
                    cfg=[];
                    cfg.dataset = files(f).name;
                    cfg.trl = trials.trl;
                    cfg.artfctdef.reject = 'complete'; % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
                    cfg.artfctdef.muscle.artifact = artefact_muscle;
                    trials_clean = ft_rejectartifact(cfg);
                    trials.artefact_muscle = setdiff(trials.trl(:,1),trials_clean.trl(:,1)); % stores onset of artefact trial
                    fprintf('\n%f percent of trials rejected',(length(trials.artefact_muscle)/length(trials.trl))*100);
                end
                
                % save new (clean) trigger file
                save(filesTrlDef(f).name,'trials');
                
            end % blocks
            
        end % subjects
        
        fprintf('\n\nData artefact rejected!\n\n');
        
    case 'blink'
        
        % parameters for SPM function
        detect_method = 'Eyes';
        correct_method = 'Berg';
        
        for s=1:length(subjects)
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);
            
            % search for input files
            files = dir(prevStep);
            
            for f=1:length(files)
                
                fprintf([ '\n\nProcessing ' files(f).name '...\n\n' ]);
                
                % copy file
                S = [];
                S.D = spm_eeg_load(files(f).name);
                S.newname = ['clean_' files(f).name];
                Dclean = spm_eeg_copy(S);
                
                % main process 1 (detect spatial confounds)
                S = [];
                S.D = Dclean;
                S.method = detect_method;
                Dclean = spm_eeg_spatial_confounds(S);
                
                % main process 2 (correct spatial confounds)
                S = [];
                S.D = Dclean;
                S.correction = correct_method;
                Dclean = spm_eeg_correct_sensor_data(S);
                Dclean.save;
                
            end % blocks
            
        end % subjects
        
        fprintf('\n\nEye blinks removed!\n\n');
        
    case 'convert+epoch' % convert data into SPM format
        
        % parameters for SPM function
        S.mode = 'epoched';
        if isfield(p,'channels') % use specified channels
            S.channels = p.channels;
        elseif isfield(p,'montage_fname') % or else use channels specified in montage file
            load([pathstem p.montage_fname]);
            S.channels = montage.labelorg;
        else % or else convert all channels
            S.channels = 'all';
        end
        
        for s=1:length(subjects)
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);
            
            % search for input files
            files = dir(prevStep);
            
            for f=1:length(files)
                
                fprintf([ '\n\nProcessing ' files(f).name '...\n\n' ]);
                
                % set input file
                S.dataset = files(f).name;
                
                % load trial structure
                load(['trlDef_' files(f).name]); % loads trial structure
                
                % set trial definition
                S.trl = trials.trl;
                S.conditionlabels = trials.labels;
                
                % main process
                D = spm_eeg_convert(S);
                
                % set manually specified bad EEG/EOG channels
                if ~isempty(badeeg{s})
                    
                    cids = D.indchannel(badeeg{s});
                    if ~isempty(cids)
                        D = D.badchannels(cids,1);
                        D.save;
                    end
                    
                end
                
                % set bad trials (if EOG artefacts present)
                if isfield(trials,'artefact_eog')
                   
                   [discard ind] = intersect(trials.trl(:,1,1),trials.artefact_eog);
                   D = D.reject(ind,1);
                   D.save;
                    
                end
                
                % set bad trials (if muscle artefacts present)
                if isfield(trials,'artefact_muscle')
                   
                   [discard ind] = intersect(trials.trl(:,1,1),trials.artefact_muscle);
                   D = D.reject(ind,1);
                   D.save;
                    
                end
                
                % add additional event info (if any)
                if isfield(trials,'events_custom')
                   
                   D.events_custom = trials.events_custom; 
                   D.save;
                    
                end
                
            end % blocks
            
        end % subjects
        
        fprintf('\n\nData converted!\n\n');
        
    case 'convert' % convert data into SPM format
        
        % parameters for SPM function
        S.mode = 'continuous';
        if isfield(p,'channels') % use specified channels
           S.channels = p.channels;
        elseif isfield(p,'montage_fname') % or else use channels specified in montage file
           load([pathstem p.montage_fname]);
           S.channels = montage.labelorg;
        else % or else convert all channels
            S.channels = 'all';
        end
        
        for s=1:length(subjects)
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);
            
            % search for input files
            files = dir(prevStep);
            
            for f=1:length(files)
                
                fprintf([ '\n\nProcessing ' files(f).name '...\n\n' ]);
                
                % set input file
                S.dataset = files(f).name;
                
                % set output file
                S.outfile = [outputstem 'spm8_' strtok(files(f).name,'.')];
                
                % main process
                D = spm_eeg_convert(S);
                
                % set manually specified bad EEG/EOG channels
                if ~isempty(badeeg{s})
                    
                    cids = D.indchannel(badeeg{s});
                    if ~isempty(cids)
                        D = D.badchannels(cids,1);
                        D.save;
                    end
                    
                end
                
            end % blocks
            
        end % subjects
        
        fprintf('\n\nData converted!\n\n');
        
    case 'convert_CTF' % convert data into SPM format
        
        % parameters for SPM function
        S.mode = 'continuous';
        S.saveorigheader = 1;
        if isfield(p,'channels') % use specified channels
           S.channels = p.channels;
        elseif isfield(p,'montage_fname') % or else use channels specified in montage file
           load([pathstem p.montage_fname]);
           S.channels = montage.labelorg;
        else % or else convert all channels
            S.channels = 'all';
        end
         
        for s=1:length(subjects)
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            for b=1:length(blocksin{s})
                
                fprintf([ '\n\nProcessing ' blocksin{s}{b} '...\n\n' ]);
                
                % change to input directory
                filePath = [rawpathstem subjects{s} '_' experimenterID{s} '_' dates{s} '_' blocksin{s}{b} '.ds/'];
                cd(filePath);
                
                % search for input files
                file = dir(prevStep);      
                
                % set input file
                S.dataset = [filePath file.name];
                
                % set output file
                outputstem = [pathstem subjects{s} '/'];
                if ~exist(outputstem,'dir'); mkdir(outputstem); end
                S.outfile = [outputstem 'spmeeg_' blocksout{s}{b}];
                
                % main process
                D = spm_eeg_convert(S);
                
                % set manually specified bad EEG/EOG channels
                if ~isempty(badeeg{s})
                    
                    cids = D.indchannel(badeeg{s});
                    if ~isempty(cids)
                        D = D.badchannels(cids,1);
                        D.save;
                    end
                    
                end
                
                % automatically detect bad channels for CTF system
                thresh = 3; % criterion is 3 standard deviations from mean of channel standard deviations
                modality = D.modality;
                if ~iscell(modality); modality = {modality}; end
                for mod=1:length(modality)
                    cids = D.selectchannels(modality{mod});
                    stdev = [];
                    for i=1:length(cids)
                        dataSelected = D(cids(i),:);
                        stdev(i) = std(dataSelected,1);
                    end
                    stdev = stdev-mean(stdev); % demean
                    stdev2 = std(stdev,1);
                    cids = cids(find(stdev>stdev2*thresh));
                    if ~isempty(cids)
                        D = D.badchannels(cids,1);
                        D.save;
                    end
                end
                
            end % blocks
            
        end % subjects
        
        fprintf('\n\nData converted!\n\n');
        
    case 'epoch'
        
        % parameters for SPM function
        S.bc = 0;
        
        for s=1:length(subjects)
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);
            
            % search for input files
            files = dir(prevStep);
            
            for f=1:length(files)
                
                fprintf([ '\n\nProcessing ' files(f).name '...\n\n' ]);
                
                % set input file
                S.D = files(f).name;
                
                % load trial structure
                load(['trlDef_' files(f).name]); % loads trial structure
                
                % set trial definition
                S.trl = trials.trl;
                S.conditionlabels = trials.labels;
                
                % main process
                D = spm_eeg_epochs(S);
                
                % set bad trials (if EOG artefacts present)
                if isfield(trials,'artefact_eog')
                   
                   [discard ind] = intersect(trials.trl(:,1,1),trials.artefact_eog);
                   D = D.reject(ind,1);
                   D.save;
                    
                end
                
                % set bad trials (if muscle artefacts present)
                if isfield(trials,'artefact_muscle')
                   
                   [discard ind] = intersect(trials.trl(:,1,1),trials.artefact_muscle);
                   D = D.reject(ind,1);
                   D.save;
                    
                end
                
                % add additional event info (if any)
                if isfield(trials,'events_custom')
                   
                   D.events_custom = trials.events_custom; 
                   D.save;
                    
                end
                
            end % blocks
            
        end % subjects
        
        fprintf('\n\nData epoched!\n\n');
        
    case 'downsample' % downsample data
        
        % parameters for SPM function
        S.fsample_new = p.fs_new; 
        
        for s=1:length(subjects)
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);
            
            % search for input files
            files = dir(prevStep);
            
            for f=1:length(files)
                
                fprintf([ '\n\nProcessing ' files(f).name '...\n\n' ]);
                
                % set input file
                S.D = files(f).name;
                
                % main process
                spm_eeg_downsample(S);
                
            end
            
        end % subjects
        
        fprintf('\n\nData downsampled!\n\n');
        
    case 'filter' % filter data
        
        % parameters for SPM function
        S.type = 'but';
        S.order = 5;
        S.band = p.filter;
        S.freq = p.freq;
        S.dir = 'twopass';
        
        for s=1:length(subjects)
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);
            
            % search for input files
            files = dir(prevStep);
            
            for f=1:length(files)
                
                fprintf([ '\n\nProcessing ' files(f).name '...\n\n' ]);
                
                % set input file
                S.D = files(f).name;
                
                % main process
                spm_eeg_filter(S);
                
            end
            
        end % subjects
        
        fprintf('\n\nData filtered!\n\n');
              
    case 'merge' % merge data
        
        % parameters for SPM function
        S.recode = 'same';
        
        for s=1:length(subjects)
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);
            
            % search for input files
            files = dir(prevStep);
            
            % set input files for merging (note: spm_eeg_merge requires file
            % names as a character array)
            mergePrimary = repmat(' ',1,50); % this is to ensure that the first file in list has the sensor locations etc. (because otherwise that information won't be retained in merged file)
            mergeSecondary = repmat(' ',1,50);
            countPrimary = 0;
            countSecondary = 0;
            for f=1:length(files)              
                
                fprintf([ '\n\nProcessing ' files(f).name '...\n\n' ]);
                
                D = spm_eeg_load(files(f).name);
                if ~isempty(D.fiducials)
                    countPrimary = countPrimary + 1;
                    mergePrimary(countPrimary,1:length(files(f).name)) = files(f).name;
                else
                    countSecondary = countSecondary + 1;
                    mergeSecondary(countSecondary,1:length(files(f).name)) = files(f).name;
                end
            
            end
            
            if countPrimary && countSecondary
                files2merge = [mergePrimary; mergeSecondary];
            else
                files2merge = mergePrimary;
            end
            
            S.D = files2merge;

            % main process
            spm_eeg_merge(S);
            
        end % subjects
        
        fprintf('\n\nData merged!\n\n');
    
     case 'rereference'
        
        % parameters for SPM function
        S.keepothers = 'no';
        
        if isfield(p,'montage_fname')
            load([pathstem p.montage_fname]); % load montage
        else
            error('Please supply montage filename');
        end
        
        for s=1:length(subjects)
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);
            
            % search for input files
            files = dir(prevStep);
            
            for f=1:length(files)
                
                fprintf([ '\n\nProcessing ' files(f).name '...\n\n' ]);
                
                % find out bad channels for rereferencing
                D = spm_eeg_load(files(f).name);
                bad = D.badchannels;
                
                montage_new = montage; % get custom template montage
                
                montage_new.tra(bad,:) = 0; % leave bad channels untouched
                montage_new.tra(:,bad) = 0;
                for i=1:length(bad)
                    montage_new.tra(bad(i),bad(i)) = 1;
                end
                
                good = setxor(D.selectchannels('EEG'),bad); % rereference (excluding bad channels)
                montage_new.tra(good,good) = -1/length(good);
                for i=1:length(good)
                    montage_new.tra(good(i),good(i)) = (length(good)-1)/length(good);
                end
                
                S.montage = montage_new; % set new montage
                
                % set input file
                S.D = files(f).name;
                
                % main process
                spm_eeg_montage(S);
                
            end
            
        end % subjects
        
        fprintf('\n\nData (EEG) rereferenced!\n\n');    
        
    case 'baseline' % baseline correct data
        
        % parameters for SPM function
        S.timewin = [p.preBase p.postBase];
        
        for s=1:length(subjects)
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);
            
            % search for input files
            files = dir(prevStep);
            
            for f=1:length(files)
                
                fprintf([ '\n\nProcessing ' files(f).name '...\n\n' ]);
                
                % set input file
                S.D = files(f).name;
                
                % main process
                spm_eeg_bc(S);
                
            end
            
        end % subjects
        
        fprintf('\n\nData baseline corrected!\n\n');
    
    case 'baseline_custom' % baseline correct data (using baseline data from another file). ASSUMES MERGED FILE!
        
        % parameters
        time = [p.preBase p.postBase];
        baselineFile = p.baselineFile; % full name of merged file (inc. directory) with baseline data
        % e.g. '/Users/edizsohoglu/Documents/projects/MEG4/data_MEG/preprocess_scene/*/cefdspmeeg_01.mat'
        
        for s=1:length(subjects)
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);
            
            % search for input files
            files = dir(prevStep);
                
            fprintf([ '\n\nProcessing ' files.name '...\n\n' ]);
            
            % load input file
            D = spm_eeg_load(files.name);
            D_new = D.clone(['b2' fname(D)], size(D));
            D_new(:,:,:) = D(:,:,:);
            
            % load baseline file
            D_baseline = spm_eeg_load(strrep(baselineFile,'*',subjects{s}));
            
            if D_new.ntrials ~= D_baseline.ntrials
                error('Number of trials in target file and baseline file do not match!');
            end
            
            % main process
            startSample = D.indsample(time(1)/1000);
            endSample = D.indsample(time(2)/1000);
            for i = 1:D_new.ntrials
                tmp = mean(D_baseline(:, startSample:endSample, i), 2);
                D_new(:, :, i) = D_new(:, :, i) - repmat(tmp, 1, D_new.nsamples);
            end
            D_new.save;
            
        end % subjects
        
        fprintf('\n\nData baseline corrected!\n\n');
        
    case 'artefact' % detect artefacts
        
        % parameters for SPM function
        S.badchanthresh = p.badchanthresh; 
        S.methods.fun = 'threshchan';
        S.methods.settings.threshold = p.thresh;
        S.methods.channels = p.artefactchans;
        
        for s=1:length(subjects)
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);
            
            % search for input files
            files = dir(prevStep);
            
            for f=1:length(files)
                
                fprintf([ '\n\nProcessing ' files(f).name '...\n\n' ]);
                
                % set input file
                S.D = files(f).name;
                
                % main process
                spm_eeg_artefact(S);
                
            end
            
        end % subjects
        
        fprintf('\n\nData artefact rejected!\n\n');
        
    case 'average' % average data
        
        % parameters for SPM function
        if p.robust == 1; 
            S.robust.savew = 1; % Robust averaging
            S.robust.bycondition = p.bycondition;
            S.robust.ks = 2;
            S.robust.removebad = 1;
        else
            S.robust = 0; % No robust averaging
        end
        S.circularise = 1; % gives phase locking value (PLV) when averaging phase values of TF data
        
        for s=1:length(subjects)
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);
            
            % search for input files
            files = dir(prevStep);
            
            for f=1:length(files)
                
                fprintf([ '\n\nProcessing ' files(f).name '...\n\n' ]);
                
                % load input file
                S.D = spm_eeg_load(files(f).name);
                
                % main process
                spm_eeg_average(S);
                
            end
            
        end % subjects
        
        fprintf('\n\nData averaged!\n\n');
        
    case 'weight' % Compute contrast on averaged data (can perform groups of contrasts such that each group gets written to a different file- if parameters specified as cell arrays)
                      
        % parameters for SPM function
        S.c = p.contrast_weights; % contrast matrix (one row per contrast)
        S.label = p.contrast_labels; % cell array of contrast labels (one row per contrast)
        S.WeightAve = 0;        
        
        for s=1:length(subjects)
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);
            
            % search for input files
            files = dir(prevStep);
            
            for f=1:length(files)
                
                fprintf([ '\n\nProcessing ' files(f).name '...\n\n' ]);
                
                % load input file
                S.D = spm_eeg_load(files(f).name);
                
                if strcmp(S.label,'All') % average across all conditions
                    S.c = repmat(1/S.D.ntrials,1,S.D.ntrials);
                end
                
                % main process
                spm_eeg_contrast(S);
                
            end
            
        end % subjects
            
        
        fprintf('\n\nData contrasted!\n\n');
        
    case 'sort' % sort conditions according to specified order
        
        % parameters for SPM function
        S.condlist = p.conditions;
        S.save = 1;
        
        for s=1:length(subjects) % for multiple subjects
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);
            
            % search for input files
            files = dir(prevStep);
            
            for f=1:length(files)
                
                fprintf([ '\n\nProcessing ' files(f).name '...\n\n' ]);
                
                % set input files
                S.D = files(f).name;
                
                % main process
                spm_eeg_sort_conditions(S);
                
            end
            
        end
        
        fprintf('\n\nConditions sorted!\n\n');
        
    case 'adjust_sensor_positions' % currently set to average sensor positions across trials (can also be used to discard trials with excessive head movements, see spm_eeg_megheadloc)
        
        % parameters for SPM function
        S.rejectbetween = 0;
        S.rejectwithin = 0;
        S.losttrack = 'preserve';
        S.correctsens = 1;
        S.save = 0;
        
        for s=1:length(subjects) % for multiple subjects
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);
            
            % search for input files
            files = dir(prevStep);
            
            for f=1:length(files)
                
                fprintf([ '\n\nProcessing ' files(f).name '...\n\n' ]);
                
                % set input files
                D = spm_eeg_load(files(f).name);
                S.D = D;
                
                % main process
                D = spm_eeg_megheadloc(S);
                
                % save new file
                D_new = D.clone(['a' fname(D)], size(D));
                D_new(:,:,:) = D(:,:,:);
                D_new.save;
                
            end
            
        end
        
        fprintf('\n\nSensor positions fixed!\n\n');
        
    case 'grand_average' % assumes merged files (i.e. one per subject)
        
        % parameters for SPM function
        S.weighted = 0;
        
        % set input files for averaging
        subjects2average = [];
        for s=1:length(subjects) % for multiple subjects
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);
            
            % search for input files
            files = dir(prevStep);
            
            fprintf([ '\n\nProcessing ' files.name '...\n\n' ]);
            
            subjects2average = [subjects2average; [filePath '/' files.name]];
            
        end
        
        % set input files
        S.D = subjects2average;
        
        % setup output filename (have to do this despite what grandmean()
        % documentation says!)
        S.outfile = [pathstem 'g' files.name];
        
        % main process
        spm_eeg_grandmean(S);
        
        fprintf('\n\nData grand averaged!\n\n');
        
    case 'TF' % perform time-frequency analysis
        
        % parameters for SPM function
        S.frequencies = p.freqs;
        S.method = p.method;
        if strcmp('morlet',p.method)
            S.settings.ncycles = p.ncycles;
            S.settings.subsample = 12;
            S.phase = p.phase;
        elseif strcmp('mtmconvol',p.method)          
            S.settings.timeres = p.timeres;
            S.settings.timestep = p.timestep;
            S.settings.freqres = p.freqres;
        elseif strcmp('morletFT',p.method)
            S.settings.ncycles = p.ncycles;
            S.settings.timestep = 50;
            S.phase = p.phase;
        end
        if isfield(p,'tf_chans')
            S.channels = p.tf_chans;
        end

        for s=1:length(subjects)
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);
            
            % search for input files
            files = dir(prevStep);
            
            for f=1:length(files)
                
                fprintf([ '\n\nProcessing ' files(f).name '...\n\n' ]);
                
                % set input file
                S.D = files(f).name;
                
                % main process
                spm_eeg_tf(S);
                
            end
            
        end % subjects
        
        fprintf('\n\nData subjected to time-frequency analysis!\n\n');
        
    case 'TF_rescale' % perform baseline correction for time-frequency analysis
        
        % parameters for SPM function
        S.method = 'rel';
        S.timewin = [p.preBase_tf p.postBase_tf];
        
        for s=1:length(subjects)
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);
            
            % search for input files
            files = dir(prevStep);
            
            for f=1:length(files)
                
                fprintf([ '\n\nProcessing ' files(f).name '...\n\n' ]);
                
                % set input file
                S.D = files(f).name;
                
                % main process
                spm_eeg_tf_rescale(S);
                
            end
            
        end % subjects
        
        fprintf('\n\nTime-frequency analysis baseline corrected!\n\n');
        
    case 'combineplanar' % combine MEGPLANAR sensor pairs using RMS (if no MEGPLANAR sensors present,the resulting file is identical to the input file but filename is prepended with 'p')
        
        correct = p.correctPlanar;
        if correct
            pre = p.preBase; % specify baseline period for baseline correction of RMSed planar gradiometer data
            post = p.postBase;
        end
        
        for s=1:length(subjects)
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);
            
            % search for input files
            files = dir(prevStep);
            
            for f=1:length(files)
                
                fprintf([ '\n\nProcessing ' files(f).name '...\n\n' ]);
                
                % set input file
                S.D = files(f).name;
                
                % main process
                D = spm_eeg_load(S.D);
                D_new = D.clone(['p' fname(D)], size(D));
                
                cind = D.selectchannels('MEGPLANAR');
                if ~isempty(cind)
                    if strcmp(D.transformtype,'time') % for time domain data
                        D(cind,:,:) = es_combineplanar3(D(cind,:,:));
                        D_new(:,:,:) = D(:,:,:);
                        
                        if correct % baseline correction
                            startSample = D_new.indsample(pre/1000);
                            endSample = D_new.indsample(post/1000);
                            D_new(:,:,:) = D_new(:,:,:) - repmat(mean(D_new(:,startSample:endSample,:),2),[1 D.nsamples 1]);
                        end
                    else % for time-frequency domain data
                        D(cind,:,:,:) = es_combineplanar3(D(cind,:,:,:));
                        D_new(:,:,:,:) = D(:,:,:,:);
                        
                        if correct % baseline correction
                            startSample = D_new.indsample(pre/1000);
                            endSample = D_new.indsample(post/1000);
                            D_new(:,:,:,:) = D_new(:,:,:,:) - repmat(mean(D_new(:,:,startSample:endSample,:),2),[1 1 D.nsamples 1]);
                        end
                    end
                end
                
                D_new.save;
                
            end
            
        end % subjects
        
        fprintf('\n\nMEGPLANAR data combined!\n\n');
                   
    case 'image' % assumes merged files (i.e. one per subject)
        
        % parameters for SPM function
        S.mode = p.imagetype;

        % other parameters
        modality = p.mod;
 
        for m=1:length(modality) % for multiple modalities
            
            fprintf([ '\n\nCurrent imaging modality = ' modality{m} '...\n\n' ]);
            
            S.channels = modality{m}; % set imaging modality
            
            for s=1:length(subjects) % for multiple subjects
                
                fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
                
                % change to input directory
                filePath = [pathstem subjects{s}];
                cd(filePath);
                
                % search for input files
                files = dir(prevStep);
                
                fprintf([ '\n\nProcessing ' files.name '...\n\n' ]);
                
                % set input file
                S.D = [pathstem subjects{s} '/' files.name];
                
                % load file (required for subsequent steps)
                D = spm_eeg_load(S.D);
                
                % main process
                spm_eeg_convert2images(S);

                % move created folder to modality-specific folder
                copyfile(strtok(files.name,'.'),S.channels);
                rmdir(strtok(files.name,'.'),'s');
                
            end
            
        end
        
        fprintf('\n\nData converted to image files!\n\n');
        
    case 'smooth'
        
        % parameters for SPM function
        smooth = [p.xSmooth p.ySmooth p.zSmooth];
        
        % other parameters
        modality = p.mod;
        
        for m=1:length(modality)
            
            fprintf([ '\n\nCurrent imaging modality = ' modality{m} '...\n\n' ]);
            
            for s=1:length(subjects) % for multiple subjects
                
                fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
                
                % change to input directory (image folder level)
                filePath = [pathstem subjects{s} '/' modality{m}];
                cd(filePath);                
                
                files = dir('condition*.nii');
                
                for f=1:length(files)
                        
                    fprintf([ '\n\nProcessing ' files(f).name '...\n\n' ]);
                    
                    % set input file
                    inputFile = files(f).name;
                    
                    % main process
                    spm_smooth(inputFile,['sm_' inputFile],smooth);
%                     Vi = spm_vol({inputFile ['sm_' inputFile]});
%                     Vo = spm_vol(['sm_' inputFile]);
%                     spm_imcalc(Vi,Vo,'((i1+eps).*i2)./(i1+eps)'); % reinsert NaNs for voxels outside space-time volume

                end
                
            end
            
        end
        
        fprintf('\n\nImage files smoothed!\n\n');

    case 'image_1D' % combine MEGPLANAR sensor pairs using RMS (if no MEGPLANAR sensors present,the resulting file is identical to the input file but filename is prepended with 'p')
        
        % other parameters
        modality = p.mod;
        
        for s=1:length(subjects)
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);
            
            % search for input files
            files = dir(prevStep);
            
            for f=1:length(files)
                
                fprintf([ '\n\nProcessing ' files(f).name '...\n\n' ]);
                
                % load input file
                D = spm_eeg_load(files(f).name);
                conditions = D.conditions;
                
                for m=1:length(modality)
                    
                    fprintf([ '\n\nCurrent imaging modality = ' modality{m} '...\n\n' ]);
                    
                    for c=1:length(conditions)
                        
                        % main process
                        foldername = [pathstem subjects{s} '/' modality{m} '/1D_type_' conditions{c} '/'];
                        if ~exist(foldername,'dir')
                            mkdir(foldername);
                        end
                        filename = [foldername sprintf('trial%04d.img', c)];
                        
                        chaninds = D.selectchannels(modality{m});
                        coninds = D.indtrial(conditions{c});
                        data = D(chaninds,:,coninds);
                        data = sqrt(nanmean(data.^2,1))';
                        time = D.time*1000;
                        
                        N     = nifti;
                        dat   = file_array(filename, [length(data), 1], 'FLOAT64-LE');
                        N.dat = dat;
                        N.mat = [...
                            diff(time(1:2))   0               0        time(1);...
                            0                 1               0        0;...
                            0                 0               1        0;...
                            0                 0               0        1];
                        N.mat(1,4) = N.mat(1,4) - N.mat(1,1);
                        N.mat_intent = 'Aligned';
                        create(N);
                        
                        N.dat(:, :) = data;
                        
                    end
                    
                end
                
            end
            
        end % subjects
        
        fprintf('\n\nData converted to 1D image files!\n\n');
        
    case 'mask' % make mask for image files so that stats are done only on time window of interest
                % takes as input one unmasked image file (from any subject
                % and condition)
        
        % parameters for SPM function
        S.timewin = [p.preImageMask p.postImageMask]; % time window outside which to mask
        
        % other parameters
        modality = p.mod;
        
        for m=1:length(modality)
            
            fprintf([ '\n\nCurrent imaging modality = ' modality{m} '...\n\n' ]);
            
            % change to input directory (image folder level)
            filePath = [pathstem subjects{1} '/' modality{m}];
            cd(filePath);
                      
            % search for input files
            files = dir(prevStep);
            files = files(1); % Use first file found (assume all images [e.g. from different conditions] have same format)
            
            fprintf([ '\n\nProcessing ' files.name '...\n\n' ]);
            
            % set input (and output) files
            S.image = files.name;
            S.outfile = [pathstem sprintf([modality{m} '_mask_%d_%dms.nii'],S.timewin(1),S.timewin(2))];
            
            % main process
            if ~isempty(strfind(S.image,'tf')) || ~isempty(strfind(S.image,'tph')) || ~isempty(strfind(S.image,'power')) || ~isempty(strfind(S.image,'phase'))
                fprintf('\nMaking time-frequency mask...\n');
                spm_eeg_mask_TF_es(S); % Use custom function if TF data
            else
                fprintf('\nMaking time-domain mask...\n');
                spm_eeg_mask(S); % else use normal SPM function for time-domain data
            end
            
        end
        
        fprintf('\n\nMask for image files made!\n\n');
        
    case 'firstlevel'
        
        % parameters for SPM function
        windows = p.windows;
        
        % other parameters
        modality = p.mod;
             
        for w=1:size(windows,1)
            
            S.window = windows(w,1:2);
            fprintf(sprintf('\n\nWindow %d-%d...\n\n',S.window(1),S.window(2)));
            
            for s=1:length(subjects) % for multiple subjects
                
                fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
                
                for m=1:length(modality)
                    
                    folders = dir([pathstem subjects{s} '/' modality{m} '/type*']);
                    
                    for fl=1:length(folders)
                        
                        fprintf([ '\n\nCurrent condition = ' folders(fl).name '...\n\n' ]);
                        
                        % change to input directory
                        filePath = [pathstem subjects{s} '/' modality{m} '/' folders(fl).name];
                        cd(filePath);
                        
                        % search for input file
                        file = dir(prevStep);
                        
                        fprintf([ '\n\nProcessing ' file.name '...\n\n' ]);
                        
                        % set input file
                        S.images = file.name;
                        S.Pout = filePath;
                        
                        % main process
                        spm_eeg_firstlevel(S);
                        tmpfile = dir('*con*.img');
                        movefile(tmpfile.name,sprintf(['t%d_%d_' file.name],S.window(1),S.window(2)));
                        tmpfile = dir('*con*.hdr');
                        file.name(end-3:end) = '.hdr';
                        movefile(tmpfile.name,sprintf(['t%d_%d_' file.name],S.window(1),S.window(2)));
                        
                    end
                    
                end
                
            end
            
        end
        
        fprintf('\n\nImage files windowed!\n\n');
        
    case 'DSS'
        
        modality = p.mod;
        if isfield(p,'file2dss'); file2dss= p.file2dss; else file2dss = []; end % if DSS using components from another file
        if isfield(p,'comp2dss'); comp2dss = p.comp2dss; else error('Please specify "comp2dss"'); end % number of components to keep/remove
        if isfield(p,'time2dss'); time2dss = p.time2dss; else time2dss = []; end
        if isfield(p,'biasWeights2dss'); W_bias = p.biasWeights2dss; else W_bias = []; end
        if isfield(p,'baselineWeights2dss'); W_baseline = p.baselineWeights2dss; else W_baseline = []; end
        if isfield(p,'cond2dss'); cond2dss= p.cond2dss; else cond2dss = []; end % if DSS using components from another file
        comp2plot = [1:10]; % components to plot
        
        for s=1:length(subjects)
            
            fprintf([ '\n\nCurrent subject = ' subjects{s} '...\n\n' ]);
            
            % create folder for saving DSS outputs
            dssFolderName = [pathstem 'DSS2/'];
            if ~exist(dssFolderName,'dir'); mkdir(dssFolderName); end
            
            % change to input directory
            filePath = [pathstem subjects{s}];
            cd(filePath);
            
            % search for input files
            files = dir(prevStep);
            
            fprintf([ '\n\nProcessing ' files.name '...\n\n' ]);
            
            % set input file
            S.D = files.name;
            
            % load spm file
            D = spm_eeg_load(S.D);
            D_new = D.clone(['n2' fname(D)], size(D));
            D_new(:,:,:) = D(:,:,:);
            
            % get trial info
            if isempty(cond2dss)
                conditions = D_new.condlist;
            else
                conditions = cond2dss; 
            end
            
            for m=1:length(modality)
                
                chanind = setdiff(D_new.selectchannels(modality{m}),D_new.badchannels); % select only good channels
                if ~isempty(time2dss);
                    timeind = unique(D_new.indsample([time2dss(1):time2dss(2)]/1000));
                else
                    timeind = D_new.indsample(D.time);
                end
                timeaxis = D_new.time;                
                
                % this part for computing DSS matrices
                data_current = {}; data_current_trimmed = {}; bias = {}; trialsind_current = {}; trialsind_current_good = {};
                cov0 = zeros(length(chanind),length(chanind));
                cov1 = zeros(length(chanind),length(chanind));
                for c=1:length(conditions)
                    
                    if isempty(strfind(conditions{c},'*'))
                        trialsind_current{c} = D.indtrial(conditions{c});
                    else
                        [discard conind] = regexp(D.conditions,conditions{c},'match');
                        trialsind_current{c} = find(~cellfun(@isempty,conind));
                    end
                    
                    % get data from current condition in loop
                    data_current{c} = D_new(chanind,:,trialsind_current{c}); % channels * time * trials
                    data_current{c} = permute(data_current{c},[2 1 3]); % rearrange data as time * channels * trials

                    % discard outliers
                    ntrials = size(data_current{c},3);
                    trialsind = nt_find_outlier_trials(nt_demean2(data_current{c}),3); % first-pass artefact rejection
                    trialsind2 = nt_find_outlier_trials(nt_demean2(data_current{c}(:,:,trialsind)),3); % second-pass artefact rejection
                    trialsind = trialsind(trialsind2);
                    data_current{c} = data_current{c}(:,:,trialsind);
                    trialsind_current_good{c} = trialsind_current{c}(trialsind);
                    fprintf('\nRejecting %2.0f percent of trials. %2.0f remaining.\n',((ntrials-length(trialsind))/ntrials)*100,length(trialsind));
                    
                    % zero out time samples outside 'time2dss' range (if specified)
                    data_current_trimmed{c} = zeros(size(data_current{c}));
                    data_current_trimmed{c}(timeind,:,:) = data_current{c}(timeind,:,:);
                    
                    % prepare biasing function (in this case, average across trials)
                    bias{c} = mean(data_current_trimmed{c},3);
                    
                    % compute covariances
                    cov0 = cov0 + nt_cov(data_current{c});
                    if isempty(file2dss)
                        cov1 = cov1 + nt_cov(bias{c});
                    else
                        load(strrep(file2dss,'*',[subjects{s} '_' modality{m}]),'cov1');
                        cov1 = cov1 + cov1;
                    end
                    
                end
                
                % compute DSS matrix that maximizes repeatability
                [todss,pwr0,pwr1] = nt_dss0(cov0,cov1);
                
                if ~isempty(W_bias) % apply contrast to subspace of repeatable components
                    
                    % get repeatable components
                    z_tmp = {};
                    for c=1:length(conditions)
                        z_tmp{c} = nt_mmat(data_current{c},todss);
                    end
                    
                    % keep subspace of repeatable components
                    NKEEP = 15;
                    for c=1:length(conditions)
                        z_tmp{c} = z_tmp{c}(:,1:NKEEP,:);
                    end
                    
                    contrast_baseline = zeros(size(mean(z_tmp{1},3)));
                    for c=1:length(conditions)                                      
                        contrast_baseline = contrast_baseline + mean(z_tmp{c},3) .* W_baseline(c);                        
                    end
                    
                    contrast_bias = zeros(size(mean(z_tmp{1},3)));
                    for c=1:length(conditions)                                      
                        contrast_bias = contrast_bias + mean(z_tmp{c},3) .* W_bias(c);                        
                    end
                    
                    % compute c0 covariance
                    cov0 = zeros(NKEEP,NKEEP);
                    for c=1:length(conditions)                                      
                        cov0 = cov0 + nt_cov(contrast_baseline);                        
                    end
                    
                    % compute c1 covariance
                    cov1 = zeros(NKEEP,NKEEP);
                    for c=1:length(conditions)                                      
                        cov1 = cov1 + nt_cov(contrast_bias);                           
                    end
                    
                    % compute 2nd order DSS matrix
                    [todss2,pwr0,pwr1] = nt_dss0(cov0,cov1);
                    
                end
                
                % this part for projecting back to sensor space
                z = {};
                if ~isempty(W_bias); cov = zeros(size(todss2,1),size(todss,1));
                else; cov = zeros(size(todss,1),size(todss,1)); end
                for c=1:length(conditions)
                    
                    % keep selected DSS components
                    if ~isempty(W_bias)
                        z{c} = nt_mmat(z_tmp{c},todss2);
                    else
                        z{c} = nt_mmat(data_current{c},todss);
                    end
                    cov = cov + nt_xcov(z{c},data_current{c});
                    
                end              
                data_dss = {};
                for c=1:length(conditions)
                    data_dss{c} = nt_mmat(z{c}(:,comp2dss,:),cov(comp2dss,:));
                end
                
                % this part for visualizing components
                
                % plot timecourse (for each condition separately)
                figure;
                for c=1:length(conditions)                                                         
                    for i=1:length(comp2plot)                        
                        subplot(1,length(comp2plot),i);
                        plot(timeaxis,mean(z{c}(:,comp2plot(i),:),3)); hold all
                        xlim([timeaxis(1) timeaxis(end)]);
                    end                                      
                end
                tightfig;
                figname = sprintf([dssFolderName '%s_%s_Components'],[strtok(S.D,'.') '_' subjects{s}],modality{m});
                saveas(gcf,figname);
                close;
                                
                % plot topos
                figure;
                dataTopo.freq = 1;
                dataTopo.label = D_new.chanlabels(chanind); % {1 x N}
                dataTopo.dimord = 'chan_freq';
                cfg = [];
                cfg.layout='CTF275.lay';
                cfg.zlim = 'maxabs';
                cfg.marker = 'off';
                cfg.comment = 'no';
                for i=1:length(comp2plot)
                    subplot(2,length(comp2plot),i+length(comp2plot));
                    dataTopo.powspctrm = cov(comp2plot(i),:)'; % topo vector as [N x 1]
                    ft_topoplotER(cfg,dataTopo);                
                end
                tightfig;
                figname = sprintf([dssFolderName '%s_%s_Topos'],[strtok(S.D,'.') '_' subjects{s}],modality{m});
                saveas(gcf,figname);
                close;
                
                % plot RMS of projected and original data
                figure(300);
                for i=1:length(comp2plot) % compare projections when retaining (cumulatively) different number of components
                    subplot(length(comp2plot)+2,1,i);
                    data2plot = nt_mmat(z{c}(:,1:comp2plot(i),:),cov(1:comp2plot(i),:));
                    plot(D_new.time,es_rms(mean(data2plot,3),2),'LineWidth',2); title(sprintf('%d components',comp2plot(i)));
                    xlim([timeaxis(1) timeaxis(end)]);
                end
                subplot(length(comp2plot)+2,1,length(comp2plot)+1);
                data2plot = nt_mmat(z{c}(:,:,:),cov(:,:));
                plot(timeaxis,es_rms(mean(data2plot,3),2),'LineWidth',2); hold all; title('All components');
                xlim([timeaxis(1) timeaxis(end)]);
                subplot(length(comp2plot)+2,1,length(comp2plot)+2);
                plot(timeaxis,es_rms(mean(data_current{c},3),2),'LineWidth',2); hold all; title('Original');
                xlim([timeaxis(1) timeaxis(end)]);
                tightfig;
                figname = sprintf([dssFolderName '/%s_%s_Projected'],[strtok(S.D,'.') '_' subjects{s}],modality{m});
                saveas(gcf,figname);
                close;
                
                % plot power of components
                figure;
                bar(cumsum(pwr1)/sum(pwr1));
                figname = sprintf([dssFolderName '%s_%s_Power'],[strtok(S.D,'.') '_' subjects{s}],modality{m});
                saveas(gcf,figname);
                close;
                
                % save DSSed data
                for c=1:length(conditions)
                    data_dss{c} = permute(data_dss{c},[2 1 3]); % rearrange back to channels * time * trials
                    D_new(chanind,:,trialsind_current_good{c}) = data_dss{c}; % replace data with DSSed data
                    trialsind_bad = setdiff(trialsind_current{c},trialsind_current_good{c});
                    if ~isempty(trialsind_bad); D_new = D_new.badtrials(trialsind_bad,1); end % set bad trials in SPM data
                end
                D_new.save;
                
                % save .mat file of DSS variables and options
                filename = sprintf([dssFolderName '%s_%s_Variables'],[strtok(S.D,'.') '_' subjects{s}],modality{m});
                if ~isempty(W_bias)
                    save(filename,'cov0','cov1','todss','pwr0','pwr1','z','cov','bias','comp2dss','time2dss','file2dss','W_baseline','W_bias','todss2','z_tmp');
                else
                    save(filename,'cov0','cov1','todss','pwr0','pwr1','z','cov','bias','comp2dss','time2dss','file2dss');
                end
                
            end
            
        end % subjects
        
        fprintf('\n\nData DSSed!\n\n');
        
    otherwise
        
        fprintf('\n\nProcess not found. Please try again!\n\n');
        
end