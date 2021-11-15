function es_plotIndividual(data,conditions)

% Plots data for individual observations as well as means
% Ed Sohoglu 2021

nSubjects = size(data,1);
nConditions = size(data,2);

figure(gcf); hold all
for c=1:nConditions
    scatter(c+(rand(nSubjects,1)-.5)*.4, data(:,c), 10, [1 .5 0], 'filled');
    scatter(c+(rand(nSubjects,1)-.5)*.4, data(:,c), 10, [1 .5 0], 'filled');
    scatter(c+(rand(nSubjects,1)-.5)*.4, data(:,c), 10, [1 .5 0], 'filled');
end
for c=1:nConditions
    xData = linspace(c-.25,c+.25);
    yData = ones(size(xData)).*mean(data(:,c),1);
    line(xData,yData,'Color','k','LineWidth',2);
end
xlim([0 nConditions+1]);
set(gca,'xtick',1:nConditions);
set(gca,'FontSize',14);
if nargin>1
    set(gca,'xticklabel',conditions);
end