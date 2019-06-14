function [ratio_us,ratio_rama,ratio_inv,ratio_PCA,ratio_rama_1,ratio_Cocco]=sequence_mutation_Math(bias,Delta_0,M,method,plot_fig,kappa) %(scale,shift,M,kappa)

%% this is the main file used in the paper for simulating the selection and data analysis
%% last update: 2019-2-26

%[count,Delta,mean_sequence,seq_corr]
%% generating random Delta

% % scale=5;
% % Delta_0=1*randn(20,1);
% % %factor=scale*floor(20*(randn(95,1)+0.5))+1;
% % %factor=(4:1:10);
% % Delta_0(1:10)=scale*randn(10,1);
% % save FigS13_Delta;
% 
% load FigS13_Delta;
%  load Delta_first_1; 
%% The Delta for the conformational change of two PDZ structures 1be9 and 1BFE, used in our paper

%load Delta_for_conformational_change;
% Delta_0=Delta_lambda;
 

%% The Delta for the conformational change associated with the normal mode of structure 1bfe
% load generated_data/20190228_single_point_mutation_from_eigenmodes_1bfe.mat;
% Delta_0=eigen_array_many(3,:)'; 

% load paper_figure_Comparison_Delta_first_50

%% Initialization

%M=500000;
scenario=1;  % 0, for -1,1;  1 for 0, 1;  others
%plot_fig=0;

% method, two-entry variable [a,b], how to perform selection: 
% .  [a=1,b=1], reweight each sequence using the quadratic free energy,x^2
% .  [a=1,b=0], reweight each sequence using the quartic free energy, x^4
%    [a=0,b=1], Window-selection to get some sequences, then each sequence has the same weight
%    [a=0,b=0], Threshold-selection to get some sequences, then each one has the same weight      
%method=[1,1]; %reweight each sequence with the quadratic free energy

remove_conservation=0; % no remove by default

if ~exist('M')
    M=50000; %the default value
end

if exist('kappa') % it is OK not to pass kappa
    d=1/sqrt(kappa); % assume quardatic selection;  d is the selection width;
end


%% 
addpath('./functions')

if ~exist('figure','dir')
    mkdir figure;
end

if ~exist('Data','dir')
    mkdir hh
end
%%

count_0=length(Delta_0);

%% random expectation of the recovery
rand_expec=sqrt(2/(pi*count_0))*sum(abs(Delta_0))/sqrt(sum(Delta_0.*Delta_0));


%% generating sequences, different scenarios
sequence_array=zeros(M,count_0);
sum_K=zeros(M,1);

if scenario==1 % beta=0,1
    
    for j=1:M
    sequence=floor(rand(count_0,1)+0.5);
    sum_K(j)=sum(sequence.*Delta_0);
    sequence_array(j,:)=sequence;
    end

else
    if scenario==0  % beta=-1,1

        for j=1:M
        sequence=floor(rand(count_0,1)-0.5)*2+1;
        sum_K(j)=sum(sequence.*Delta_0);
        sequence_array(j,:)=sequence;
        end
    
    else   % beta=0,1,2,..,scenario (non-uniform)
        
        for j=1:M
        sequence=round(rand(count_0,1)*scenario)/scenario-0.5;
        sum_K(j)=sum(sequence.*Delta_0);
        sequence_array(j,:)=sequence;
        end
        
    end
        
end
   
% figureParameter
% hist(sum_K,100,'k');
% xlim([0 800]);
%  fig_name='./figure/sum_K.eps';
%  figurePostTreat

if plot_fig
 plot_distribution(sum_K,[-150 150],100,'$\sum_lS_l\Delta_l$');
end

%% selection and covariance
if ~exist('kappa')
     d=0.3*sqrt(var(sum_K)); %d=0.3*std is used in the paper
end
%bias=0; %bias=3.5 is used in the paper
gamma=mean(sum_K)+bias*sqrt(var(sum_K))-0;

[seq_corr_0,mean_sequence_0]=selection_covariance(sequence_array,sum_K,d,gamma,method);
count_0=length(mean_sequence_0);

%% remove the conserved sites, not run by default

if remove_conservation
    index=mean_sequence_0<1.02; % no remove
    %index=mean_sequence_0<0.98; %  remove
    mean_sequence=mean_sequence_0(index);
    sequence_array_new=sequence_array(:,index);
    Delta=Delta_0(index);
    count=sum(index);

    [seq_corr,mean_sequence]=selection_covariance(sequence_array_new,sum_K,d,gamma,weight_or_not);

else 
    seq_corr=seq_corr_0;
    mean_sequence=mean_sequence_0;
    Delta=Delta_0;
    count=count_0;
end
    
%% the variance at each site
var_l=zeros(1,count);
for j=1:count
   var_l(j)=seq_corr(j,j); 
end


if plot_fig
    
figureParameter
plot(1:count,1./var_l,'*r');
axis tight
%ylim([0 9]);
a1=xlabel('Residue');
%ylim([-5.3*10^4 3*10^4]);
a2=ylabel('$ 1/\sigma_l^2$');
fig_name='./figure/sigma-beta-2.eps';
figurePostTreat

end


%% SCA, using the conservation phi to reweight the matrix, with square root of eigenvector as predictor

q0=0.5;
epsilon=10^(-20); % used to avoid singularity
phi_conserv=-log((1-mean_sequence)./(epsilon+mean_sequence))+log((1-q0)/q0);

% figureParameter
% f1=plot(1:count,phi_conserv,'*k');
% a1=xlabel('$l$');
% a2=ylabel('$\phi_l$');
% %a2=ylabel('$1/\sigma_l^2$');
% %h1=legend('$\phi_l$');
% axis tight
% %ylim([0 6]);
% fig_name='./figure/weight.eps';
% figurePostTreat

seq_corr_SCA=zeros(count,count);
for j=1:count
   for k=1:count
      seq_corr_SCA(j,k)=abs(seq_corr(j,k)).*abs(phi_conserv(j)*phi_conserv(k));
       
   end
end

select_data=seq_corr_SCA;
[NormVector_SCA,orderEigValue_SCA]=orderedEigSystem(select_data,0);
Delta1=Delta;
Delta_weight=Delta1./(sqrt(sum(Delta1.*Delta1)));

vector=sqrt(abs(NormVector_SCA));
for j=1:count
vector(:,j)=vector(:,j)/sqrt(sum(vector(:,j).*vector(:,j)));
end
ratio_rama=(abs(Delta_weight')*abs(vector))';

vector_1=abs(NormVector_SCA);
ratio_rama_1=(abs(Delta_weight')*abs(vector_1))';


if plot_fig
    figureParameter
    f1=plot(1:count,ratio_rama,'+k');
    xlim([0.5,100.5]);
    ylim([0,1]);
    fig_name='./figure/correct_ratio_SCA.eps';
    figurePostTreat
end


%% ICOD: inverse the matrix, and setting the diagonal terms to be zero

seq_corr_ICOD=inv(seq_corr);
for j=1:count
    for k=j:count
       if j==k
       seq_corr_ICOD(j,j)=0; % this works equally well compared with more elaborated choice  
       end
    end
end

select_data=seq_corr_ICOD;
[NormVector_ICOD,orderEigValue_ICOD]=orderedEigSystem(select_data,0);

Delta1=Delta;
Delta_weight=Delta1./(sqrt(sum(Delta1.*Delta1)));
ratio_inv=(abs(Delta_weight')*abs(NormVector_ICOD))';
%improved_prediction_analysis(Delta,NormVector_ICOD(:,1));

% 
if plot_fig
figureParameter
f1=plot(1:count,ratio_inv,'+r');
%xlim([0.5,100.5]);
xlabel('Eigenvector');
ylabel('Recovery');
ylim([0,1]);
fig_name='./figure/correct_ratio_ICOD.eps';
figurePostTreat
end
%%  Our reweighting method, divide the covariance matrix Cij by the variance of at i and j respectively. Also,  set the off-diagonal terms to zero
    
reweight_seq_corr=zeros(count,count);
for j=1:count
    for k=j:count
       if j==k
       reweight_seq_corr(j,j)=0; % this works equally well compared with more elaborated choice  
       else
       reweight_seq_corr(j,k)=seq_corr(j,k)/(seq_corr(j,j)*seq_corr(k,k));
       reweight_seq_corr(k,j)=reweight_seq_corr(j,k);
       end
    end
end


select_data=reweight_seq_corr;
[NormVector_RW,orderEigValue_RW]=orderedEigSystem(select_data,0);
Delta1=Delta;
Delta_weight=Delta1./(sqrt(sum(Delta1.*Delta1)));
ratio_us=(abs(Delta_weight')*abs(NormVector_RW))';
%improved_prediction_analysis(Delta,NormVector_RW(:,end));


%% the naive method:PCA

select_data=seq_corr;
[NormVector_PCA,orderEigValue_PCA]=orderedEigSystem(select_data,0);
Delta1=Delta;
Delta_weight=Delta1./(sqrt(sum(Delta1.*Delta1)));
ratio_PCA=(abs(Delta_weight')*abs(NormVector_PCA))';


%% test of Eq.(9) from Cocco11
conserv_weight=sqrt(var_l).*sqrt(var_l');
rescaled_Corr=seq_corr./conserv_weight;

max_N=count;
select_data=rescaled_Corr;
[NormVector_cocco,orderEigValue_cocco]=orderedEigSystem(select_data,0);

predict_Cocco=zeros(max_N,max_N);
for j=1:max_N
predict_Cocco(:,j)=sqrt(max_N*(1/orderEigValue_cocco(j,1)-1))*NormVector_cocco(:,j)./sqrt(var_l)';
end

Delta1=Delta;
Delta_weight=Delta1./(sqrt(sum(Delta1.*Delta1)));
predict_Cocco=predict_Cocco./(sqrt(sum(predict_Cocco.*predict_Cocco)));
ratio_Cocco=(abs(Delta_weight')*abs(predict_Cocco))';

%% confirm the prediction about entries of the inverse matrix
% inv_corr=inv(seq_corr);
% kappa=0.09;
% esti_inv=zeros(count,count);
% for j=1:count
%     for k=1:count
%           esti_inv(j,k)=kappa*Delta1(j)*Delta1(k);
%        if k==j
%           esti_inv(j,k)=kappa*Delta1(j)*Delta1(k)+1/mean_sequence(j)+1/(1-mean_sequence(j));  
%        end
%         
%     end
% end
% 
% diag_index=zeros(count,2);
% for j=1:count
%     diag_index(j,1)=j;diag_index(j,2)=j;
% end
% 
% 
% figureParameter
% f1=plot(1:count,esti_inv(20,:),'-r',1:count,inv_corr(20,:),':k')
% xlim([0 count+1]);
% a1=xlabel('$\Delta$');
% %a2=ylabel('Distribution');
% fig_name='./figure/predicted_inv.eps';
% figurePostTreat


%% eigenmode analysis for a given method, method exploration
if plot_fig
max_N=count;

%select_data=reweight_seq_corr;
select_data=seq_corr_ICOD;

[NormVector_corr,orderEigValue_corr]=orderedEigSystem(select_data,0);

figureParameter
image(1:count,1:count,10*abs(select_data));
colorbar;
%set(gca, 'clim', [-0.25 0.25]);
%set(gca, 'clim', [-1 1]);
a1=xlabel('$x$');
a2=ylabel('$y$');
fig_name='./figure/corr1.jpg';
figurePostTreat
%set(gca,'YDir','normal')
%print(fig_name,'-r100','-djpeg');


figureParameter
f1=plot(1:max_N,orderEigValue_corr(:,1),'*r');
a1=xlabel('Mode index: $j$');
%a2=ylabel('Eigenvalue: $\lambda_j$');
%axis tight
xlim([1 count]);
fig_name='./figure/corr_eig.eps';
figurePostTreat

Delta1=Delta;
Delta_weight=Delta1./(sqrt(sum(Delta1.*Delta1)));
ratio=(abs(Delta_weight')*abs(NormVector_corr))';


% %% for the Rama's method
%  vector=sqrt(abs(NormVector_corr));
% for j=1:count
% vector(:,j)=vector(:,j)/sqrt(sum(vector(:,j).*vector(:,j)));
% end
% ratio=(abs(Delta_weight')*abs(vector))';
% 
figureParameter
bar(1:count,ratio,1,'Facecolor','k','EdgeColor','k');
a1=xlabel('Mode index');
%a2=ylabel(' Recovery');
%xlim([0 76]);
axis tight;
ylim([0 1.03]);
fig_name='./figure/corr_vec.eps';
figurePostTreat


figureParameter
bar(1:max_N,[NormVector_corr(:,1),1.2*Delta_weight]);
%bar(1:max_N,[abs(vector(:,1)),1.5*abs(Delta_weight)]);
%set(gca,'XTICKlabel',index2);
xlabel('Residue');
h1=legend('$\nu_l^{(1)}$','$\Delta_l$');%-\Delta_6^l\;$');
set(h1,'location','northeast')
axis tight
ylim([-1 1]);
fig_name='./figure/v1.eps';
figurePostTreat

            
end
