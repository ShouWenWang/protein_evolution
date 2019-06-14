function [ratio_us_1,ratio_us_2,ratio_rama_1,ratio_rama_2,ratio_inv_1,ratio_inv_2]=sequence_mutation_Math_double_select_new(bias1,bias2)
%% a key part here is we use ICA to separate the two constraints

%% random mutation

red = [255 0 0]/255;
dark_green = [0 180 0]/255;
light_green=[0 255 0]/255;

%% initializaton 
count=100; %L
scenarios=1;  % 1,  the state is 0 or 1; 0, the state is +1 or -1;

%% generate the two constraints
% first constraint
sector_size=20;
weight1=randn(count,1);
factor1=floor(20*randn(sector_size,1))+1;
 
weight1((1:length(factor1)))=factor1;
%gamma1=shift1*sqrt(sum(Delta_1.^2));

% second constraint
weight2=randn(count,1);
factor2=floor(20*randn(sector_size,1))+1;

%weight2(count-length(factor2)+1:1:count)=factor2;
weight2(20+1:1:20+length(factor2))=factor2;
%gamma2=shift2*sqrt(sum(Delta_2.^2));
save ./Data/two_Delta_first_20_overlap_0 weight1 weight2

figure,bar(1:count,[weight1,weight2])


% Delta_1=abs(Delta_1);
% Delta_2=abs(Delta_2);
% 
%  load two_Delta_first_20_overlap_20;
% % load two_Delta_first_50;
% load paper_two_sectors

Delta_1=weight1;
Delta_2=weight2;

%% random expectation
%rand_expec=sqrt(2/(pi*count_0))*sum(abs(Delta_0))/sqrt(sum(Delta_0.*Delta_0));

%%

over_lap=sum(Delta_1.*Delta_2)/sqrt(sum(Delta_1.*Delta_1)*sum(Delta_2.*Delta_2));
% 
% figureParameter
% bar(1:count,[Delta_1,Delta_2])
% axis tight
% fig_name='./figure/two_sectors.eps';
% figurePostTreat

%M=20000000;
M=50000;
% M=50000;
sequence_array=zeros(M,count);
sum_K1=zeros(1,M);
sum_K2=zeros(1,M);

if scenarios==1
    
    for j=1:M  % 0 or 1
    sequence=floor(rand(count,1)+0.5);
    sum_K1(j)=sum(sequence.*Delta_1);
    sum_K2(j)=sum(sequence.*Delta_2);
    sequence_array(j,:)=sequence';
    end

else 
    
    for j=1:M %Ising model, -1 or +1
    sequence=floor(rand(count,1)-0.5)*2+1;
    sum_K1(j)=sum(sequence.*Delta_1);
    sum_K2(j)=sum(sequence.*Delta_2);
    sequence_array(j,:)=sequence';
    end

end


%% selection
d1=0.3*sqrt(var(sum_K1)); %d=0.3*std is used in the paper
gamma1=mean(sum_K1)-bias1*sqrt(var(sum_K1));
d2=0.3*sqrt(var(sum_K2));
gamma2=mean(sum_K2)-bias2*sqrt(var(sum_K2));

kappa1=1/d1^2;
kappa2=1/d2^2;
weighting_fun_1=@(x) sqrt(kappa1/(2*pi))*exp(-0.5*kappa1*(x-gamma1).^2);
weighting_fun_2=@(x) sqrt(kappa2/(2*pi))*exp(-0.5*kappa2*(x-gamma2).^2);
Weight=weighting_fun_1(sum_K1').*weighting_fun_2(sum_K2');

%% covariance matrix 
pseudo_ratio=0;
[seq_corr,mean_sequence]=covariance_matrix(sequence_array,Weight,pseudo_ratio);

if plot_fig
    figureParameter
    f1=pcolor(1:count,1:count,2*seq_corr);
    colorbar;
    set(gca, 'clim', [-0.05 0.05]);
    a1=xlabel('$x$');
    a2=ylabel('$y$');
    fig_name='./figure/mode_correlation.eps';
    figurePostTreat

    figureParameter
    f1=plot(1:count,mean_sequence,'*k');
    a1=xlabel('$l$');
    %a2=ylabel('$\langle S^*\rangle$');
    fig_name='./figure/mean-beta.eps';
    figurePostTreat
end

%% reweighting 1


q0=0.5;
phi_conserv=-log((1-mean_sequence)./mean_sequence)+log((1-q0)/q0);

if  plot_fig
    figureParameter
    f1=plot(1:count,phi_conserv,'*k');
    a1=xlabel('$l$');
    a2=ylabel('$\phi_l$');
    %a2=ylabel('$1/\sigma_l^2$');
    %h1=legend('$\phi_l$');
    axis tight
    %ylim([0 6]);
    fig_name='./figure/weight.eps';
    figurePostTreat
end

reweight_seq_corr=abs(seq_corr).*abs(phi_conserv*phi_conserv');
%reweight_seq_corr=(seq_corr).*(phi_conserv*phi_conserv');
% for j=1:count
%     reweight_seq_corr(j,j)=0;
% end

select_data=reweight_seq_corr;
[NormVector_corr,orderEigValue_corr]=orderedEigSystem(select_data,0);

[S,A,U,ll,info]=icaML(NormVector_corr(:,1:2)',2);
S(1,:)=S(1,:)/sqrt(sum(S(1,:).*S(1,:)));
S(2,:)=S(2,:)/sqrt(sum(S(2,:).*S(2,:)));
NormVector_corr(:,1:2)=S';

Delta_10=Delta_1;
Delta_20=Delta_2;
Delta_weight_1=Delta_10./(sqrt(sum(Delta_10.*Delta_10)));
Delta_weight_2=Delta_20./(sqrt(sum(Delta_20.*Delta_20)));

vector=sqrt(abs(NormVector_corr));
for j=1:count
vector(:,j)=vector(:,j)/sqrt(sum(vector(:,j).*vector(:,j)));
end

ratio_rama_1=(abs(Delta_weight_1')*vector);
ratio_rama_2=(abs(Delta_weight_2')*vector);





%% reweighting-2: inverse matrix


reweight_seq_corr=inv(seq_corr);
for j=1:count
    for k=j:count
       if j==k
       reweight_seq_corr(j,j)=0; % this works equally well compared with more elaborated choice  
       end
    end
end
%reweight_seq_corr=abs(reweight_seq_corr);


select_data=reweight_seq_corr;
[NormVector_corr,orderEigValue_corr]=orderedEigSystem(select_data,0);

[S,A,U,ll,info]=icaML(NormVector_corr(:,1:2)',2);
S(1,:)=S(1,:)/sqrt(sum(S(1,:).*S(1,:)));
S(2,:)=S(2,:)/sqrt(sum(S(2,:).*S(2,:)));
NormVector_corr(:,1:2)=S';

Delta_10=Delta_1;
Delta_20=Delta_2;
Delta_weight_1=Delta_10./(sqrt(sum(Delta_10.*Delta_10)));
Delta_weight_2=Delta_20./(sqrt(sum(Delta_20.*Delta_20)));
vector=abs(NormVector_corr);

ratio_inv_1=(abs(Delta_weight_1')*vector);
ratio_inv_2=(abs(Delta_weight_2')*vector);

%%  reweighting-3
reweight_seq_corr=zeros(count,count);
for j=1:count
    for k=j:count
       if j==k
       reweight_seq_corr(j,j)=0; % this works equally well compared with more elaborated choice  
       else
       reweight_seq_corr(j,k)=seq_corr(j,k)/(seq_corr(j,j)*seq_corr(k,k));
       %reweight_seq_corr(j,k)=seq_corr(j,k);
       reweight_seq_corr(k,j)=reweight_seq_corr(j,k);
       end
    end
end
%reweight_seq_corr=abs(reweight_seq_corr);

% figureParameter
% f1=pcolor(1:count,1:count,0.2*reweight_seq_corr);
% colorbar;
% xlabel('Residue');
% ylabel('Residue');
% set(gca, 'clim', [-0.05 0.05]);
% %set(gca, 'clim', [-1 1]);
% %a1=xlabel('$x$');
% %a2=ylabel('$y$');
% fig_name='./figure/mode_correlation.eps';
% figurePostTreat


select_data=reweight_seq_corr;
[NormVector_corr,orderEigValue_corr]=orderedEigSystem(select_data,0);

[S,A,U,ll,info]=icaML(NormVector_corr(:,end-1:end)',2);
S(1,:)=S(1,:)/sqrt(sum(S(1,:).*S(1,:)));
S(2,:)=S(2,:)/sqrt(sum(S(2,:).*S(2,:)));
NormVector_corr(:,end-1:end)=S';

Delta_10=Delta_1;
Delta_20=Delta_2;
Delta_weight_1=Delta_10./(sqrt(sum(Delta_10.*Delta_10)));
Delta_weight_2=Delta_20./(sqrt(sum(Delta_20.*Delta_20)));
vector=abs(NormVector_corr);

ratio_us_1=(abs(Delta_weight_1')*vector);
ratio_us_2=(abs(Delta_weight_2')*vector);


%% eigenmode analysis
if  plot_fig
    
max_N=count;
selection_index=1:max_N;%count-max_N+1:count;
%select_data=seq_corr;
select_data=reweight_seq_corr(selection_index,selection_index);
Delta_1=Delta_1(selection_index);
Delta_2=Delta_2(selection_index);


[NormVector_corr,orderEigValue_corr]=orderedEigSystem(select_data,0);
%% ICA for our reweighting
% [S,A,U,ll,info]=icaML(NormVector_corr(:,end-1:end)',2);
% S(1,:)=S(1,:)/sqrt(sum(S(1,:).*S(1,:)));
% S(2,:)=S(2,:)/sqrt(sum(S(2,:).*S(2,:)));
% NormVector_corr(:,end-1:end)=S';
%% ICA for Rama reweighting or inverse matrix
% [S,A,U,ll,info]=icaML(NormVector_corr(:,1:2)',2);
% S(1,:)=S(1,:)/sqrt(sum(S(1,:).*S(1,:)));
% S(2,:)=S(2,:)/sqrt(sum(S(2,:).*S(2,:)));
% NormVector_corr(:,1:2)=S';

figureParameter
f1=plot(1:max_N,orderEigValue_corr(:,1),'*r');
a1=xlabel('Mode index: $j$');
%a2=ylabel('Eigenvalue: $\lambda_j$');
axis tight
%ylim([0 0.6]);
fig_name='./figure/corr_eig.eps';
figurePostTreat


Delta_10=Delta_1;
Delta_20=Delta_2;


Delta_weight_1=Delta_10./(sqrt(sum(Delta_10.*Delta_10)));
Delta_weight_2=Delta_20./(sqrt(sum(Delta_20.*Delta_20)));


%% for the Rama's method
vector=abs(NormVector_corr);
% vector=sqrt(abs(NormVector_corr));
% for j=1:count
% vector(:,j)=vector(:,j)/sqrt(sum(vector(:,j).*vector(:,j)));
% end

ratio_1=(abs(Delta_weight_1')*vector);
ratio_2=(abs(Delta_weight_2')*vector);


    
figureParameter
bar(1:max_N,ratio_1,1,'Facecolor','k','EdgeColor','k');
a1=xlabel('Mode index');
%a2=ylabel(' Recovery');
%xlim([0 76]);
axis tight;
ylim([0 1.03]);
fig_name='./figure/corr_ratio_1.eps';
figurePostTreat

figureParameter
bar(1:max_N,ratio_2,1,'Facecolor','b','EdgeColor','b');
a1=xlabel('Mode index');
%a2=ylabel(' Recovery');
%xlim([0 76]);
axis tight;
ylim([0 1.03]);
fig_name='./figure/corr_ratio_2.eps';
figurePostTreat


figureParameter
bar(1:max_N,[abs(NormVector_corr(:,end)),abs(Delta_weight_2)]);
%set(gca,'XTICKlabel',index2);
xlabel('Residue');
h1=legend('$|\nu_l^{(L-1)}|$','$|\Delta_l^{(2)}|$');
set(h1,'location','north')
axis tight
%ylim([0 1]);
fig_name='./figure/v1.eps';
figurePostTreat


figureParameter
bar(1:max_N,[abs(NormVector_corr(:,end-1)),abs(Delta_weight_1)]);
%set(gca,'XTICKlabel',index2);
xlabel('Residue');
h1=legend('$\nu_l^{(L)}$','$\Delta_l^{(1)}\;$');
set(h1,'location','north')
axis tight
%ylim([0 1]);
fig_name='./figure/v2.eps';
figurePostTreat





figureParameter
f1=plot(1:max_N,mean_sequence,'*','Color',light_green);
%xlim([0 1]);
a1=xlabel('Residue');
%set(gca,'XTICKlabel',index2);
%h1=legend('$\nu_{1}^l$');%-\Delta_6^l\;$');
%set(h1,'location','northeast')
axis tight
ylim([0 1]);
fig_name='./figure/Delta.eps';
figurePostTreat

figureParameter
bar(1:max_N,Delta_2,'b');
%xlim([0 1]);
a1=xlabel('Residue');
%set(gca,'XTICKlabel',index2);
%h1=legend('$\nu_{1}^l$');%-\Delta_6^l\;$');
%set(h1,'location','northeast')
axis tight
ylim([-5 11]);
fig_name='./figure/Delta.eps';
figurePostTreat

figureParameter
f1=plot(1:count,phi_conserv,'*');
axis tight
ylim([0 1.02]);
a1=xlabel('Residue');
%a2=ylabel('$\langle S_l\rangle_*$');
hold on;

plot(1:max_N,Delta_1/5,'or');
%a2=ylabel('$\nu_l^{(N)}$');
xlabel('Residue');
%set(gca,'XTICKlabel',index2);
%h1=legend('$\phi_l$','$\Delta_6^l\;$');
%set(h1,'location','northeast')
axis tight
ylim([-3 3]);
fig_name='./figure/v_end.eps';
figurePostTreat
hold off

% 


figureParameter
set(fig,'defaultAxesColorOrder',[dark_green; red]);
yyaxis left
plot(1:max_N,phi_conserv,'*');
a1=xlabel('Residue');
ylim([-2 2]);
yyaxis right
plot(1:count,mean_sequence,'o');
axis tight
ylim([0 1]);
fig_name='./figure/v_end.eps';
figurePostTreat



% figureParameter
% bar(1:max_N,[NormVector_corr(:,1),Delta_weight]);
% %set(gca,'XTICKlabel',index2);
% h1=legend('$\nu_{1}^l$','$\Delta_l$');%-\Delta_6^l\;$');
% set(h1,'location','northeast')
% xlim([0 max_N+1]);
% fig_name='./figure/corr_mode2.eps';
% figurePostTreat
% 
% 
% 
% 
% figureParameter
% bar(1:count,[Delta_weight,phi]);
% xlabel('$l$');
% %set(gca,'XTICKlabel',index2);
% h1=legend('$\Delta_l$','$\phi_l$');%-\Delta_6^l\;$');
% set(h1,'location','north')
% xlim([0 count]);
% fig_name='./figure/corr_mode2.eps';
% figurePostTreat
 
            
end