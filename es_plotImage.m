function h = es_plotImage(xaxis,yaxis,data)

h = imagesc(data); % first dimension of data will be yaxis, second will be xaxis
set(gca,'YDir','normal');
polarmap;
if length(yaxis) > 10
    set(gca,'ytick',1:10:length(yaxis)); set(gca,'yticklabel',round(yaxis(1:10:end)));
else
    set(gca,'ytick',1:length(yaxis)); set(gca,'yticklabel',round(yaxis(1:end)));
end
if length(xaxis) > 10
    set(gca,'xtick',1:10:length(xaxis)); set(gca,'xticklabel',round(xaxis(1:10:end)));
else
    set(gca,'xtick',1:length(xaxis)); set(gca,'xticklabel',round(xaxis(1:end)));
end
