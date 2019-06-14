fig=figure;

set(gcf,'Units','centimeters')
%figure('Units', 'pixels', ...
%    'Position', [100 100 500 375]);
%hold on;
scaling=1.5;
width=scaling*8;height=scaling*5.5;fontsize=20;
set(gcf,'Position',[0,0,width,height]);
set(gcf,'Papersize',[width,height])
set(gcf,'paperposition',[0,0,width,height])
set(gcf, 'renderer', 'painters');
set(gcf, 'PaperPositionMode', 'auto');


% figure,
% 
% set(gcf,'Units','centimeters')
% %figure('Units', 'pixels', ...
% %    'Position', [100 100 500 375]);
% %hold on;
% %scaling=1.5;
% scaling=1.5;
% width=scaling*9;height=scaling*6;fontsize=20;
% set(gcf,'Position',[0,0,width,height]);
% set(gcf,'Papersize',[width,height])
% set(gcf,'paperposition',[0,0,width,height])
% set(gcf, 'renderer', 'painters');
% set(gcf, 'PaperPositionMode', 'auto');
% set(gca,'FontSize',15);