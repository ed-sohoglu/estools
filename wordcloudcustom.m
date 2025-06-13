function wordcloudcustom(key_words,freq)
% Ed's custom word cloud function, adapted from code written by Ikko Kimura
% (Matlab Central)

Max = 40; % Range of FontSize
Min = 6;
Rotate = 0;

% Color = 'parula';
% switch Color
%         case 'hsv'
%     colors = colormap(hsv(length(key_words)));
%         case 'jet'
%     colors = colormap(jet(length(key_words)));
%         case 'parula'
%     colors = colormap(parula(length(key_words)));
%         case 'hot'
%     colors = colormap(hot(length(key_words)));
%         case 'cool'
%     colors = colormap(cool(length(key_words)));
%         case 'spring'
%     colors = colormap(spring(length(key_words)));
%         case 'summer'
%     colors = colormap(summer(length(key_words)));
%         case 'autumn'
%     colors = colormap(autumn(length(key_words)));
%         case 'winter'
%     colors = colormap(winter(length(key_words)));
% end

% or just grey?
colors = zeros(numel(key_words),3)+[0.2510 0.2510 0.2510];
%colors(1,:) = [1 0 0]; % highlight first item

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sizes = (freq * (Max-Min))+Min; % Linearly Hmmmm...
% Initialize
pos(:,1) = [1:numel(key_words)]';
pos(:,2) = [1:numel(key_words)]';
rng(6);
pos(:,1) = pos(randperm(numel(key_words)),1);
rng(7);
pos(:,2) = pos(randperm(numel(key_words)),2);
% plot the data
xlim([min(pos(:,1))-10 max(pos(:,1))+10]); % x axis range 
ylim([min(pos(:,2))-3 max(pos(:,2))+3]); % y axis range
hold on
for i=1:length(key_words)
    if ~isnan(sizes(i))
        h = text(pos(i,1),pos(i,2),char(key_words{i}),'FontSize',sizes(i),'Color',colors(i,:),'HorizontalAlignment','center','FontName','Helvetica');
        % rotate text every other time if prompted
        if Rotate==1
            if mod(i,2)==0
                set(h,'Rotation',270);
            end
        end
    end
end
hold off
axis off
end