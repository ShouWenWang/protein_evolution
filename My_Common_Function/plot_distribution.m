function plot_distribution(data,plot_range,bin_size,name)
[Prob_0,X]=hist(data,bin_size);
Prob=Prob_0./(sum(Prob_0));

index=Prob~=0;
Prob=Prob(index);
X=X(index);

% index1=abs(X)<0.03;
% index2=abs(X)>=0.03 & Prob~=0;
% Prob1=Prob(index1);
% X1=X(index1);
% Prob2=Prob(index2);
% X2=X(index2);


figureParameter
%f1=plot(X1,Prob1,'*k',X2,Prob2,'or');
%f1=plot(X,Prob,'*b');
b=bar(X,Prob,1,'Facecolor','k','EdgeColor','k');
%b.CData(1,:) = [0 0.8 0.8];
ax=gca;
ax.YAxis.Exponent = -2;
xlim([plot_range(1) plot_range(2)]);
a1=xlabel(name);
%axis tight
%set(gca,'EdgeColor','k');
a2=ylabel('Probability');
%ylim([0 0.08]);
%set(gca, 'YTick' , [0 0.01 0.02]);
%set(gca,'yticks',[0 0.01 0.02]);
fig_name='./figure/distribution.eps';
figurePostTreat