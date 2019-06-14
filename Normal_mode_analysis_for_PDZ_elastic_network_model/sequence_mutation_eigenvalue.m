function sequence_mutation_eigenvalue

clear all

%% important parameters
PDB_name_1='1be9.pdb'; reference_index_1=301;
%PDB_name_1='1BFE.pdb'; reference_index_1=306;
cut_index_low=309;cut_index_high=398;
cutoff_distance=7.5;
factor=0.8; % perturbation factor, close to 1 is a small perturbation
M=40000; % number of samples
%M=100;
mut_prob=0.2;
%% start

[Atom_type,coords]=read_PDB_data_C_beta(PDB_name_1,reference_index_1,cut_index_low,cut_index_high);

% get the positions of the C beta atoms
N=length(Atom_type);
beta_atom_index=zeros(2,1); % dynamically expand
count=0;
for i=1:N
   if Atom_type(i)==0
       count=count+1;
       beta_atom_index(count)=i;
   end
end



%% original eigenvalues
[NormVector,orderEigValue,coord_normVector,spring,Hessen_2d]=normal_mode_computation_Cbeta_mutation(coords,Atom_type,cutoff_distance,factor);
figure,plot(1:length(orderEigValue(:,1)),orderEigValue(:,1),'*r');
xlim([1 50]);
%s=size(NormVector);
%%
soft_mode=NormVector(:,7:27);
soft_mode_eig=orderEigValue(7:27,1);

%% generate a random sequence

Atom_type_mutant=Atom_type;

eigen_array_many=zeros(6,M);
vector_array_many=zeros(6,M);
%vector_array_many_original=zeros(M,s(1),6);

sequence_array=zeros(count,M);
for j=1:M
sequence=floor(rand(count,1)+mut_prob);
Atom_type_mutant(beta_atom_index)=sequence; 
[NormVector_mutant,orderEigValue_mutant,coord_normVector_mutant,spring,Hessen_2d]=normal_mode_computation_Cbeta_mutation(coords,Atom_type_mutant,cutoff_distance,factor);
eigen_array_many(:,j)=orderEigValue_mutant(7:12,1)-orderEigValue(7:12,1);
vector_array_many(:,j)=sum(NormVector_mutant(:,7:12).*NormVector(:,7:12));
%vector_array_many_original(j,:,:)=NormVector_mutant(:,7:12);
sequence_array(:,j)=sequence;
end





eigen_array3=eigen_array_many(3,:);
plot_distribution(eigen_array3,[-0.015 0],40,'$\lambda_3-\lambda_3^{(0)}$')

eigen_array4=eigen_array_many(4,:);
plot_distribution(eigen_array4,[-0.02 0],40,'$\lambda_4-\lambda_4^{(0)}$')

eigen_array6=eigen_array_many(6,:);
plot_distribution(eigen_array6,[-0.03 0],40,'$\lambda_6-\lambda_6^{(0)}$')


% [Prob,X]=hist(eigen_array,100);
% Prob=Prob./(sum(Prob));
% 
% a1=0.0255;
% b1=-0.008619;
% c1=0.002534;
% Prob_guass=a1*exp(-((X-b1)/c1).^2);
% 
% figureParameter
% f1=plot(X,Prob,'*b',X,Prob_guass,'-r');
% %xlim([0.065 0.08]);
% a1=xlabel('$\lambda_3-\lambda_3^{(0)}$');
% a2=ylabel('Distribution');
% fig_name='./figure/eign_distribution.eps';
% figurePostTreat


%% select the eiganvalue range, and compute the correlation
%index_center=( eigen_array>0.081 & eigen_array<0.0815 );
%index_center=( eigen_array>0.073 & eigen_array<0.074 );
%index_center=( eigen_array<0.0773 );
%index_center=( eigen_array<1 ); % for all sequences
%d=0.0035;
d=0;
index_center3=( eigen_array3>-0.009-d & eigen_array3<-0.008-d );
%index_center3=( eigen_array3>-1 & eigen_array3<0 );
%index_center3=( eigen_array3>-0.0117 & eigen_array3<-0.0112 );

%index_center3_1=( eigen_array3>-0.009 & eigen_array3<-0.008 );

%index_center3_2=( eigen_array3>-0.009+d & eigen_array3<-0.0080+d );
%index_center4=( eigen_array4>-0.0116 & eigen_array4<-0.0103 );
%index_center6=( eigen_array6>-0.0161 & eigen_array6<-0.0141 );
%index_center=(  eigen_array<-0.0152 );
%index_center=( eigen_array<-0.012 ); % 

%%
index_center=index_center3;

% M_new=0;  % final data size
% select_N=100; 
% while M_new<8000  %final size
% %index_center=index_center3(1:300) |index_center3_1(1:300)  | index_center3_2(1:300) ;%| index_center4 ;%& index_center6;
% index_center=index_center0(1:select_N);
% %index_center=index_center3(1:select_N);
% M_new=sum(index_center);
% select_N=floor(select_N*1.1);
% end

%%
%index_center=index_center3_1 & index_center4;

sequence_array_new=sequence_array(:,index_center);

M_new=sum(index_center);

sum_rule_tot=sum(sequence_array,1)/(count);

sum_rule=sum(sequence_array_new,1)/(count);



plot_distribution(sum_rule,[0 1],10,'$K/L$');



[Prob,X]=hist(sum_rule,10);
Prob=Prob./(sum(Prob));

[Prob1,X1]=hist(sum_rule_tot,10);
Prob1=Prob1./(sum(Prob1));


figureParameter
f1=plot(X,Prob,'-*k',X1,Prob1,'-or');
%f1=plot(X,Prob,'-*b');
%b=bar(X,Prob);
%b.CData(1,:) = [0 0.8 0.8];
xlim([0 1  ]);
a1=xlabel('$K/L$');
%ylim([0 0.12]);
a2=ylabel('Distribution');
fig_name='./figure/distribution.eps';
figurePostTreat


%%



%figure,hist(sum_rule,20)

mean_sequence=mean(sequence_array_new,2);




%figure,hist(mean_sequence);
plot_distribution(mean_sequence,[0 1],2,'$\langle \beta_l\rangle_*$');

figureParameter
bar(1:count,mean_sequence,'r');
fig_name='./figure/corr_mode.eps';
a1=ylabel('$\langle \beta_l\rangle*$');
%set(h1,'location','south')
xlim([0 80]);
%set(gca,'YTICKLABEL',[-0.1 0 0.1 0.2 0.3 0.4 0.5]);
figurePostTreat




%[hot_index_x,hot_index_y,hotspot_data]=search_hotspot(mean_sequence,0.5,[0.12 1]);


%% covariance analysis
seq_corr=zeros(count,count);
for k=1:M_new
%     for i=1:count
%         for j=1:count
%            
%            seq_corr(i,j)=seq_corr(i,j)+(sequence_array_new(i,k)-mean_sequence(i))*(sequence_array_new(j,k)-mean_sequence(j));
% 
%         end
%     end

       seq_corr=seq_corr+(sequence_array_new(:,k)-mean_sequence)*(sequence_array_new(:,k)-mean_sequence)';

       
end

seq_corr=seq_corr./M_new;


%%


% X0=abs(mean_sequence-0.5)+1;
% 
% Y0=X0*X0';
% Y0=Y0./max(max(abs(Y0)));
% 
% %test_M=0.04*randn(count,count);
% 
% %seq_corr0=seq_corr;
% 
% seq_corr=seq_corr0.*Y0;

figureParameter
f1=pcolor(1:count,1:count,seq_corr);
colorbar;
set(gca, 'clim', [-0.15 0.15]);
%set(gca, 'clim', [-1 1]);
%a1=xlabel('$x$');
%a2=ylabel('$y$');
fig_name='./figure/mode_correlation.eps';
figurePostTreat



%% analyze the correlation matrix
% histogram without the self-correlation
seq_corr_off_diag=seq_corr;
seq_diag=zeros(count,1);
for i=1:count
    seq_corr_off_diag(i,i)=0;
    seq_diag(i)=seq_corr(i,i);
end

corr_vector=zeros(count*count,1);
for i=1:count
  corr_vector((i-1)*count+1:i*count)=seq_corr_off_diag(:,i);
end

new_corr_vector=corr_vector(abs(corr_vector)>=0 );%& abs(corr_vector)<0.24);
    
%figure, hist(new_corr_vector,100);
plot_distribution(new_corr_vector,[-0.1 0.05],40,'$C_{ll''}^-$');





%% check the correlation formula

% 
% 
% 
% new_Corr=Delta_lambda3'*Delta_lambda3;
% 
% for j=1:count
%     new_Corr(j,j)=0;
% end
% 
% new_Corr_org=seq_corr_off_diag./max(max(abs(seq_corr_off_diag)));
% 
% 
% 
% 
% %[index2,hotspot_data_org]=search_hotspot(new_Corr_org,0,[0.2 1.1]);
% 
% var_sequence=mean_sequence.*(1-mean_sequence);
% weighting=sqrt(var_sequence)*sqrt(var_sequence)';%var_sequence*var_sequence';%
% new_Corr2=new_Corr.*weighting;
% new_Corr2=new_Corr2./max(max(abs(new_Corr2)));
% 
% 
% figureParameter
% bar(1:length(index2),var_sequence(index2),'r');
% fig_name='./figure/corr_mode_var.eps';
% set(gca,'XTICKlabel',index2);
% a1=ylabel('$\sigma_l^*$');
% %set(h1,'location','south')
% xlim([0 length(index2)+1]);
% %set(gca,'YTICKLABEL',[-0.1 0 0.1 0.2 0.3 0.4 0.5]);
% figurePostTreat
% 
% figureParameter
% bar(1:length(index2),[-new_Corr_org(11,index2)',new_Corr2(11,index2)']);
% set(gca,'XTICKlabel',index2);
% h1=legend('$C_{ll''}$','$\tilde{C}_{ll''}$');
% set(h1,'location','northwest')
% xlim([0 length(index2)+1]);
% %ylim([0.15 0.6]);
% fig_name='./figure/corr_mode.eps';
% figurePostTreat
% 

%% inverse matrix
% inv_corr=inv(seq_corr);
% 
% 
% figureParameter
% f1=pcolor(1:count,1:count,inv_corr);
% colorbar;
% set(gca, 'clim', [-20 20]);
% %set(gca, 'clim', [-1 1]);
% %a1=xlabel('$x$');
% %a2=ylabel('$y$');
% fig_name='./figure/mode_correlation.eps';
% figurePostTreat
% 
% [NormVector_inv_corr,orderEigValue_inv_corr]=orderedEigSystem(inv_corr,0);
% 
% 
% figureParameter
% f1=plot(1:count,orderEigValue_inv_corr(:,1),'or');
% a1=xlabel('Mode: $k$');
% a2=ylabel('Eigenvalue: $\lambda_k^c$');
% xlim([0 80]);
% fig_name='./figure/corr_eig.eps';
% figurePostTreat
% 
% figureParameter
% f1=plot(1:count,orderEigValue_inv_corr(:,1),'or');
% a1=xlabel('Mode: $k$');
% a2=ylabel('Eigenvalue: $\lambda_k^c$');
% xlim([0 80]);
% fig_name='./figure/corr_eig.eps';
% figurePostTreat
% 
% 
% 
% Delta_lambda=Delta_lambda3;
% Delta_lambda=Delta_lambda./(sqrt(sum(Delta_lambda.*Delta_lambda)));
% 
% figureParameter
% bar(1:count,[abs(NormVector_inv_corr(:,1)),-Delta_lambda']);
% %set(gca,'XTICKlabel',index2);
% h1=legend('$\nu_1^l$','$\Delta_3^l$');%-\Delta_6^l\;$');
% set(h1,'location','north')
% ylim([0 0.6]);
% fig_name='./figure/corr_mode2.eps';
% figurePostTreat

%% reweighting of the covariance matrix 
% according to our theory
%reweight_seq_corr=seq_corr./(seq_diag*seq_diag');

% Ranganathan's approach
phi_conserv=log((1-mean_sequence)./mean_sequence);

figureParameter
bar(1:count,1./seq_diag,'r');
a1=xlabel('$l$');
%a2=ylabel('$\phi_l$');
a2=ylabel('$1/\sigma_l^2$');
%h1=legend('$\phi_l$','$1/\sigma_l^2$');
xlim([0 80]);
fig_name='./figure/weight.eps';
figurePostTreat


reweight_seq_corr=abs(seq_corr).*(phi_conserv*phi_conserv');

%%  eigenvalue analysis of Corr
%select_data=seq_corr_off_diag;






%select_data=seq_corr;
select_data=reweight_seq_corr;

[NormVector_corr,orderEigValue_corr]=orderedEigSystem(select_data,0);


figureParameter
f1=plot(1:count,orderEigValue_corr(:,1),'or');
a1=xlabel('Mode: $k$');
a2=ylabel('Eigenvalue: $\lambda_k^c$');
xlim([0 80]);
%ylim([0.01 4]);
fig_name='./figure/corr_eig.eps';
figurePostTreat


Delta_lambda=Delta_lambda3;
Delta_lambda=Delta_lambda./(sqrt(sum(Delta_lambda.*Delta_lambda)));

phi=phi_conserv./(sqrt(sum(phi_conserv.*phi_conserv)));

figureParameter
bar(1:count,[Delta_lambda',phi]);
%a2=ylabel('$-\Delta_3^l$');
a1=xlabel('Site index');
%set(gca,'XTICKlabel',index2);
h1=legend('$\Delta_3^l$','$\phi_l$');%-\Delta_6^l\;$');
set(h1,'location','southeast')
%ylim([0 0.6]);
fig_name='./figure/corr_mode2.eps';
figurePostTreat




ratio=(abs(Delta_lambda)*abs(NormVector_corr))';

figureParameter
bar(1:count,ratio,'r');
a1=xlabel('Mode: $k$');
a2=ylabel('Recovery of $\Delta_3^l$');
xlim([0 80]);
fig_name='./figure/corr_vec.eps';
figurePostTreat


%% 

mode_index=1;

%X1=X0./sqrt(sum(X0.*X0));

Delta_lambda_plus=Delta_lambda3;%-Delta_lambda4;%-Delta_lambda6;%+Delta_lambda6;
Delta_lambda_plus=Delta_lambda_plus./sqrt(sum(Delta_lambda_plus.*Delta_lambda_plus));

Delta_lambda_mis=Delta_lambda4-Delta_lambda3;
%Delta_lambda_mis=Delta_lambda_mis./sqrt(sum(Delta_lambda_mis.*Delta_lambda_mis));

new_Delta_lambda1=Delta_lambda_plus+Delta_lambda_mis;
new_Delta_lambda1=new_Delta_lambda1./sqrt(sum(new_Delta_lambda1.*new_Delta_lambda1));

Delta_lambda1=Delta_lambda1./(sqrt(sum(Delta_lambda1.*Delta_lambda1)));

Delta_lambda2=Delta_lambda2./(sqrt(sum(Delta_lambda2.*Delta_lambda2)));

Delta_lambda3=Delta_lambda3./(sqrt(sum(Delta_lambda3.*Delta_lambda3)));

figureParameter
bar(1:count,[abs(NormVector_corr(:,1)),-Delta_lambda3']);
%set(gca,'XTICKlabel',index2);
h1=legend('$\nu_1^l$','$\Delta_3^l$');%-\Delta_6^l\;$');
set(h1,'location','north')
ylim([0 0.6]);
fig_name='./figure/corr_mode2.eps';
figurePostTreat

vect1=NormVector_corr(:,1)-NormVector_corr(:,2);

vect1=vect1./sqrt(sum(vect1.*vect1));

figureParameter
bar(1:count,[-new_Delta_lambda1',-Delta_lambda1']);
%set(gca,'XTICKlabel',index2);
h1=legend('$\nu_1^l$','$\Delta_3^l$');
set(h1,'location','north')
%ylim([0 0.5]);
fig_name='./figure/corr_mode2.eps';
figurePostTreat



figureParameter
bar(1:length(index2),[NormVector_corr(index2,mode_index),X1(index2),-Delta_lambda(index2)']);
set(gca,'XTICKlabel',index2);
h1=legend('$\nu_1^l$','$X_0$','$\Delta_4^l$');
set(h1,'location','north')
%ylim([0 0.5]);
fig_name='./figure/corr_mode1.eps';
figurePostTreat

%%

[hot_index_2,hotspot_data_2]=search_hotspot(new_Corr2,0,[0.3 1.1]);





figure,plot(1:length(hot_index_org),-hotspot_data_org(2,:),'*r',1:length(hot_index_2),hotspot_data_2(3,:),'+b')


%%

new_Corr=new_Corr./max(max(abs(new_Corr)));



figureParameter
f1=pcolor(1:count,1:count,new_Corr2);
colorbar;
set(gca, 'clim', [-1 1]);
%set(gca, 'clim', [-1 1]);
%a1=xlabel('$x$');
%a2=ylabel('$y$');
fig_name='./figure/mode_correlation.eps';
figurePostTreat




%%

temp=zeros(1,count);
for j=1:count
[temp_v,max_ind]=max(abs(select_data(j,:)));
temp(j)=select_data(j,max_ind);
end
temp=temp./(sqrt(sum(temp.*temp)));


temp0=sum(select_data);
temp0=temp0./(sqrt(sum(temp0.*temp0)));

Delta_lambda=Delta_lambda./(sqrt(sum(Delta_lambda.*Delta_lambda)));



X0=mean_sequence-0.5;
X0=X0./(sqrt(sum(X0.*X0)));

figureParameter
b=bar(1:count,[X0,-Delta_lambda']);
set(b(1),'FaceColor','r');
set(b(2),'FaceColor','g');
fig_name='./figure/corr_mode.eps';
h1=legend('$\langle \beta_l\rangle-0.5\;$','$\Delta_4^l$');
set(h1,'location','north')
ylim([0 0.5]);
%set(gca,'YTICKLABEL',[-0.1 0 0.1 0.2 0.3 0.4 0.5]);
figurePostTreat


figureParameter
bar(1:count,[NormVector_corr(:,mode_index),-temp']);
fig_name='./figure/corr_mode.eps';
h1=legend('$\nu_1^l$','$S_l^{(1)}$');
set(h1,'location','north')
ylim([-0.1 0.5]);
%set(gca,'YTICKLABEL',[-0.1 0 0.1 0.2 0.3 0.4 0.5]);
figurePostTreat



[hot_index_x,hot_index_y,hotspot_data]=search_hotspot(Delta_lambda,0,[0.15 1]);

index2=sort(hot_index_x);

figureParameter
bar(1:length(index2),[NormVector_corr(index2,mode_index),-temp(index2)',-temp0(index2)',-Delta_lambda(index2)']);
set(gca,'XTICKlabel',index2);
h1=legend('$\nu_1^l$','$S_l^{(1)}$','$S_l^{(2)}$','$\Delta_4^l$');
set(h1,'location','north')
%ylim([0.15 0.6]);
fig_name='./figure/corr_mode.eps';
figurePostTreat


%% hotspot analysis, show high-correlation values 

[hot_index_x,hot_index_y,hotspot_data]=search_hotspot(seq_corr_off_diag,0,[0.015 1]);

hotspot_index=hot_index_x;
%hotspot_index=[8 9 11 12 15 28 58 64];
Amino_acid_loc=visualizing_hotspot(coords,beta_atom_index,hotspot_index,spring,reference_index_1,Atom_type);


figureParameter
f1=pcolor(1:length(hot_index_x),1:length(hot_index_y),hotspot_data);
colorbar;
set(gca, 'clim', [-0.1 0.1]);
a1=xlabel('$x$');
a2=ylabel('$y$');
fig_name='./figure/mode_correlation.eps';
figurePostTreat



%% % identify hot spot pairs
% hot_spot_pair=zeros(hot_spot_count,2);
% hot_spot_pair(:,1)=real(hot_spot_index_temp);
% hot_spot_pair(:,2)=imag(hot_spot_index_temp);
% 
% %% 
% hot_spot_index=[11 15 28 58 64]';
% sequence_array_new(hot_spot_index,:);
% 
% hot_spot=sequence_array_new(hot_spot_index,:)-mean_sequence;
% 
% hot_spot_mean=mean(hot_spot,2);
% 
% hot_spot_N=length(hot_spot_index);
% hot_spot_corr=zeros(hot_spot_N,hot_spot_N);
% 
% for i=1:hot_spot_N
%     for j=1:hot_spot_N
%         hot_spot_corr(i,j)=mean(hot_spot(i,:).*hot_spot(j,:));
%     end
% end
% 

