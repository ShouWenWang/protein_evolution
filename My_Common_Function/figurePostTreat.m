scale2=1;
if exist('f1','var')
set(f1,'linewidth',1.5,'MarkerSize',4) % original, 1.5;  4
end

if exist('f2','var')
set(f2,'linewidth',1.5,'MarkerSize',4)
end

if exist('a1','var')
set(a1,'FontSize',22*scale2,'Interpreter','latex')
end

if exist('a2','var')
set(a2,'FontSize',22*scale2,'Interpreter','latex')
end

if exist('a3','var')
set(a3,'FontSize',22*scale2,'Interpreter','latex')
end
%default northeast.  can change 

if exist('h1','var')
set(h1,'interpreter','latex','fontsize',20*scale2)
end

 set(gca,'FontSize', 18*scale2);
 if exist('gca0','var')
     set(gca0,'FontSize',18*scale2);
 end

%standard output form

%% for Matlab_F1V1
if exist('fig_name','var')
if exist('my_path','var')
   fig_path=strcat(my_path,'/',fig_name);
   print('-depsc2','-painters',fig_path);
else 
   print('-depsc2','-painters',fig_name);
end
end
%print -depsc2 -painters ./figure/force.eps
% 
% set(gca,'YDir','normal')
% print(fig_name,'-r400','-djpeg');

%%  useful plot commands

%set(gca, 'XTick' , [1 100 10000]);
%control the position to show XTick,  (it shows its default value at this
%position

%set(gca, 'YTick' , [1 100 10000]);
%b=bar(X,Prob,0.9,'Facecolor','k','EdgeColor','k');
% axis tight

%% two axis

% figureParameter
% left_color = [124 0 0]/255;
% right_color = [0 128 0 ]/255;
% set(fig,'defaultAxesColorOrder',[left_color; right_color]);
% 
% x=1:10;
% y = x.^2;
% yyaxis left
% bar(x,y);
% 
% z = 10./x;
% yyaxis right
% f1=plot(x,z);
% fig_name='./figure/test.eps';
% figurePostTreat

%% bar plot: change bar width and color
% b=bar(X,Prob,1,'Facecolor','k','EdgeColor','k');