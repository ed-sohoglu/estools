function es_plotTimeSeries(data,timeAxis,masks,maskColours,lineStyles,er,masksUnc)
% es_plotTimeSeries(data,timeAxis,masks,maskColours,lineStyles)
%
% Plots time-series data and labels significant timepoints for multiple contrasts with horizontal
% coloured lines (one for each contrasts).
% data = n timeseries X t timepoints
% masks = cell array of masks (vectors with t timepoints) denoting significant time samples
% maskColours = cell array of colours for masks  e.g. {'r' 'b'}
% lineStyles = cell array of line styles, one for each time-series (e.g. {'b-' 'r-'})
% er = error bars as n timeseries X t timepoints
% masksUnc = same as 'masks' but denoting uncorrected significance
% NB- masks need to contain 1 for significant time samples and NaN for non
% significant samples (won't work if zeros)
%
% Ed Sohoglu 2020

if nargin < 7
     masksUnc = [];
end

if nargin < 6
     er = [];
end
    
for i=1:size(data,1)
    if isempty(er)
        plot(timeAxis,data(i,:),lineStyles{i},'LineWidth',2); hold all
    else
        boundedline(timeAxis,data(i,:),er(i,:),lineStyles{i},'alpha'); hold all
    end
end

limits = ylim;
ytickSpacing = mean(diff(get(gca,'ytick')));
for m=1:length(masks)
    if ~isempty(masksUnc)
        line(timeAxis,limits(1)+masksUnc{m}*ytickSpacing*m*.2,'LineWidth',1,'Color',maskColours{m});
    end
    line(timeAxis,limits(1)+masks{m}*ytickSpacing*m*.2,'LineWidth',4,'Color',maskColours{m});
end

