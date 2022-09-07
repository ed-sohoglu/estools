function h = es_plotSummary(data,conditions,type,within)

% Plots means and sems of data in nSubjects x Nconditions format
% 'conditions' = cell array of factors/levels e.g. conditions{1} = {'Attend' 'Unattend'}; conditions{2} = {'Pred' 'Unpred'};
% 'type' = 'bar' or 'line'
% 'within' = 1 to remove between subject variance
% NB currently only works for <= 2 factor designs
% Ed Sohoglu 2020

if nargin<4
    within = 0;
elseif nargin<3
    type = 'line';
elseif nargin<2
    conditions = {};
end

nSubjects = size(data,1);
for i=1:numel(conditions)
    design(i) = numel(conditions{i});
end

if within
    data = es_removeBetween(data);
end

means = mean(data,1);
sems = std(data,0,1) / sqrt(nSubjects);

means = es_reshape(means,design);
sems = es_reshape(sems,design);

if strcmp(type,'bar')
    h = bar(means); % bar plots first dimension of matrix along ticks (second dimension is nested within ticks and linked with legend)
    %% This part for adding error bars to grouped bar graph
    hold on
    nGroups = size(means, 1);
    nBars = size(means, 2);
    % Calculating the width for each bar group
    groupWidth = min(0.8, nBars/(nBars + 1.5));
    for i = 1:nBars
        x = (1:nGroups) - groupWidth/2 + (2*i-1) * groupWidth / (2*nBars);
        errorbar(x, means(:,i), zeros(nGroups,1), sems(:,i), '.','Color','k');
    end
    hold off
    %%
    if numel(design) > 1
        legend(conditions{2});
    end
    set(gca,'xtick',1:design(1));
    set(gca,'xticklabel',conditions{1});
    set(gca,'ylim',[min(means(:))-max(abs(sems(:))) max(means(:))+max(abs(sems(:)))]);
elseif strcmp(type,'line')
    means = means';
    sems = sems';
    h = errorbar(means,sems,'linewidth',2);
    if numel(design) > 1
        legend(conditions{1});
        xlim([0 design(2)+1]);
        set(gca,'xtick',1:numel(conditions{2}));
        set(gca,'xticklabel',conditions{2});
    else
        xlim([0 design(1)+1]);
        set(gca,'xtick',1:numel(conditions{1}));
        set(gca,'xticklabel',conditions{1});
    end
end