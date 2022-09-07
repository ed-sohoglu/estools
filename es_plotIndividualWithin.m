function es_plotIndividualWithin(data,conditions)

% Plots data for individual observations as well as means
% Ed Sohoglu 2021

nSubjects = size(data,1);
nConditions = size(data,2);

figure(gcf); hold all
plot(data');
plot(nanmean(data,1),'k','LineWidth',4);
xlim([0 nConditions+1]);
set(gca,'xtick',1:nConditions);
set(gca,'FontSize',14);
if nargin>1
    set(gca,'xticklabel',conditions);
end