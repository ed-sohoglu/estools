function es_plotIndividual(data,conditions,colours)

% Plots data for individual observations as well as means and sems
% conditions = cell array of condition names
% colours = cell array and can be labels e.g. 'r', 'b' or vectors like [1 .5 0]
% Ed Sohoglu 2021

if nargin<3
    for c=1:numel(conditions)
        colours{c} = [1 .5 0]; % by default use orange
    end
end

nSubjects = size(data,1);
nConditions = size(data,2);

figure(gcf); hold all

% plot individual observations
for c=1:nConditions
    scatter(c+(rand(nSubjects,1)-.5)*.4, data(:,c), 10, colours{c}, 'filled');
    scatter(c+(rand(nSubjects,1)-.5)*.4, data(:,c), 10, colours{c}, 'filled');
    scatter(c+(rand(nSubjects,1)-.5)*.4, data(:,c), 10, colours{c}, 'filled');
end

% plot mean and sem
for c=1:nConditions
    mn = mean(data(:,c),1);
    sem = std(data(:,c),0,1) / sqrt(nSubjects);
    
%     % add horizontal line showing mean
%     xData = linspace(c-.1,c+.1);
%     yData = ones(size(xData)).*mn;
%     line(xData,yData,'Color',colours{c},'LineWidth',2);
    
    % bars and errorbars
    bar(c,mn,'FaceAlpha',0,'EdgeColor',colours{c},'LineWidth',2);
    errorbar(c,mn,sem,'LineStyle','none','Color',colours{c},'LineWidth',2,'CapSize',12);
end

% labels and aesthetics
xlim([0 nConditions+1]);
set(gca,'xtick',1:nConditions);
set(gca,'FontSize',14);
if nargin>1
    set(gca,'xticklabel',conditions);
end