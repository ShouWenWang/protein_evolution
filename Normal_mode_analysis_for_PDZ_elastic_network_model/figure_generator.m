

%%

%PDB_name_1='1be9.pdb'; reference_index_1=301;
PDB_name_1='1BFE.pdb'; reference_index_1=306;
cut_index_low=309;cut_index_high=398;


cutoff_distance=7.5;
index=3; % second softest mode that is not zero

coords=read_PDB_data_C_alpha(PDB_name_1,reference_index_1,cut_index_low,cut_index_high);

figure;plot3(coords(:,1),coords(:,2),coords(:,3),'b.')

[NormVector,orderEigValue,coord_normVector,spring]=normal_mode_computation_Calpha(coords,cutoff_distance);

for i=index
    make_movie(i,NormVector,coords,spring);
end

plot_eigenmode(index,NormVector,coords,spring);
plot_mode_correlation(index,NormVector);

%% Cbeta ANM

%PDB_name_1='1be9.pdb'; reference_index_1=301;
PDB_name_1='1BFE.pdb'; reference_index_1=306;
cut_index_low=309;cut_index_high=398;
cutoff_distance=7.5;

index=9; % second softest mode that is not zero


[Atom_type,coords]=read_PDB_data_C_beta(PDB_name_1,reference_index_1,cut_index_low,cut_index_high);


[NormVector,orderEigValue,coord_normVector,spring]=normal_mode_computation_Cbeta(coords,Atom_type,cutoff_distance);
 
make_movie(index,NormVector,coords,spring);
plot_eigenmode(index,NormVector,coords,spring);
plot_mode_correlation(index,NormVector);

%% mutation of the C beta atoms



%%
[NormVector,orderEigValue,coord_normVector,spring]=normal_mode_computation_Cbeta(coords,Atom_type,cutoff_distance);
 
 index=9;
 plot_eigenmode(index,NormVector,coords,spring);
plot_mode_correlation(index,NormVector);

 plot_eigenmode(index,NormVector_mutant,coords,spring);
plot_mode_correlation(index,NormVector_mutant);

%%  hot spot for conformational change

[Norm_Involvement,displace_vector]=conformational_change_Cbeta;

[NormVector,orderEigValue,coord_normVector,spring,]=normal_mode_computation_Cbeta_mutation(coords,Atom_type,cutoff_distance,factor);

%%  study the PCA 



figureParameter
f1=pcolor(1:count,1:count,seq_corr);
colorbar;
set(gca, 'clim', [-0.1 0.1]);
a1=xlabel('$x$');
a2=ylabel('$y$');
fig_name='./figure/mode_correlation.eps';
figurePostTreat

[NormVector_corr,orderEigValue_corr]=orderedEigSystem(seq_corr);
figureParameter,bar(NormVector_corr(:,1))

C=zeros(count,count);

for k=1:1
for i=1:count
    for j=1:count
C(i,j)=C(i,j)+NormVector_corr(i,k)*NormVector_corr(j,k);
    end
end
end

figureParameter
f1=pcolor(1:count,1:count,C);
colorbar;
set(gca, 'clim', [-0.1 0.1]);
a1=xlabel('$x$');
a2=ylabel('$y$');
fig_name='./figure/mode_correlation.eps';
figurePostTreat


%% combinatorial
N=100;
y0=zeros(N,1);
for k=1:N
y0(k)=nchoosek(N,k)/nchoosek(N,N/2);
end


x0=(1:N)./N;
figure,plot(x0,log(y0))


%%
count=80; %L
M1=3;M2=5;
%N=20;
L=count;

for k=(M1+M2+1):L-2

tot=nchoosek(L-2,k-M1-M2)+nchoosek(L-2,k-M1)+nchoosek(L-2,k-M2)+nchoosek(L-2,k);
p12=nchoosek(L-2,k-M1-M2)/tot;
p13=(nchoosek(L-3,k-M1-1)+nchoosek(L-3,k-M1-M2-1))/tot;
p23=(nchoosek(L-3,k-M2-1)+nchoosek(L-3,k-M1-M2-1))/tot;
aver1=(nchoosek(L-2,k-M1-M2)+nchoosek(L-2,k-M1))/tot;
aver2=(nchoosek(L-2,k-M1-M2)+nchoosek(L-2,k-M2))/tot;
aver3=(nchoosek(L-2,k-M1-M2)*(k-M1-M2)+nchoosek(L-2,k)*k+nchoosek(L-2,k-M2)*(k-M2)+nchoosek(L-2,k-M1)*(k-M1))/(tot*(L-2));
var1(k)=aver1*(1-aver1);
var2(k)=aver2*(1-aver2);
var3(k)=aver3*(1-aver3);


covar12(k)=p12-aver1*aver2;
covar13(k)=p13-aver1*aver3;
covar23(k)=p23-aver2*aver3;

aver1_array(k)=aver1;
aver2_array(k)=aver2;
aver3_array(k)=aver3;


end

x0=(M1+M2+1):L-2;
%figure,plot(x0,-covar12(x0),'-r',x0,sqrt(var1(x0).*var2(x0)),'-*k');

figureParameter
f1=plot(x0/L,aver1_array(x0)/M1,'-r',x0/L,aver2_array(x0)/M2,'-*k',x0/L,aver3_array(x0),'-ob');
a1=xlabel('$\gamma_1/L$');
h1=legend('$\langle \beta_1\rangle/n_1$','$\langle \beta_2\rangle/n_2$','$\langle \beta_3\rangle/n_3$');
set(h1,'location','northwest');
fig_name='./figure/mode_correlation.eps';
figurePostTreat



figureParameter
f1=plot(x0/L,covar12(x0),'-r',x0/L,covar13(x0),'-*k',x0/L,covar23(x0),'-ob');
a1=xlabel('$\gamma_1/L$');
h1=legend('$C_{12}$','$C_{13}$','$C_{23}$');
set(h1,'location','southeast');
fig_name='./figure/mode_correlation.eps';
figurePostTreat



figure,plot(x0,-covar12(x0),'-r',x0,3*var1(x0).*var2(x0)*M1*M2/L,'-*k');

figure,plot(x0,-covar13(x0),'-r',x0,3*var1(x0).*var3(x0)*M1/L,'-*k');


figure,plot(x0,-covar23(x0),'-r',x0,3*var2(x0).*var3(x0)*M2/L,'-*k');

figure,plot(x0,covar23(x0)./covar12(x0),'-r');
xlim([30 60]);

%% beta=-1,1
%%
%clear all;
L=1000; %L
M1=2:2:100;
M2=10*(1:5);
%N=20;

N=length(M1);

%gamma_0=0:2:L-M1-M2-4;
Covar_array=zeros(N,3,3);
mean_seq_array=zeros(N,3);

for l=1:length(M2)
for k=1:N
    
    [Covar,mean_seq]=sequence_mutation_analytics_plus_minus_one(L,M1(k),M2(l),0);
     Covar_array(k,:,:)=Covar;
     mean_seq_array(k,:)=mean_seq;

end
  plot_data{l,1}=M1.*M2(l)/L^2;%/(M1.^2+M2(l)^2+L-2);
  plot_data{l,2}=Covar_array(:,1,2);
  plot_data{l,3}=M1./(M1.^2+M2(l)^2+L-2);
  plot_data{l,4}=Covar_array(:,1,3);
  
end


figure, plot(plot_data{1,1},plot_data{1,2},'+r',plot_data{2,1},plot_data{2,2},'ok',plot_data{3,1},plot_data{3,2},'sb',plot_data{4,1},plot_data{4,2},'xg',plot_data{5,1},plot_data{5,2},'-cyan');


a1=xlabel('$n_1$');
a2=ylabel('$C_{1,2}$');
h1=legend('$n_2=6$','$n_2=10$','$n_2=20$','$n_2=30$','$n_2=50$');
set(h1,'location','northeast');
fig_name='./figure/mode_correlation.eps';
figurePostTreat

%figure,plot(x0,-covar12(x0),'-r',x0,sqrt(var1(x0).*var2(x0)),'-*k');

figureParameter
f1=plot(gamma_0/L,atanh(mean_seq_array(:,1))/M1,'-r',gamma_0/L,atanh(mean_seq_array(:,2))/M2,'*k',gamma_0/L,atanh(mean_seq_array(:,3)),'ob');
a1=xlabel('$\gamma/L$');
h1=legend('$f(\langle \alpha_1\rangle)\;\;$','$f(\langle \alpha_2\rangle)$','$f(\langle \alpha_3\rangle)$');
set(h1,'location','northwest');
legend boxoff
fig_name='./figure/mode_correlation.eps';
figurePostTreat

%% partition function approach
L=100;
Delta=zeros(1,L);
sector=3*(1:10);
Delta(1:length(sector))=sector;
im=sqrt(-1);
kappa=100;
gamma=0;

F=@(h) h*gamma-log(cosh(h*Delta));
G=@(h) exp(-F(h));


q = integral(@free_energy,-100,100);


%% analytical solution
L=1000;
M1=6; M2=10;
lambda0=100;
constr_mean=@(h) -lambda0+tanh(M1*h)*M1+tanh(M2*h)*M2+(L-1)*tanh(h);

constr_cov=@(h) -lambda0+tanh(M1*h)*M1+tanh(M2*h)*M2+(L-1)*tanh(h)+2*M1./sinh(2*M1*h)+2*M2./sinh(2*M2*h);

Free_energy=@(h) -lambda0*h+log(cosh(M1*h))+log(cosh(M2*h))+(L-1)*log(cosh(h));
Free_energy_cov=@(h) -lambda0*h+log(cosh(M1*h))+log(cosh(M2*h))+(L-1)*log(cosh(h))-log(tanh(M1*h))-log(tanh(M2*h));


h=0.001:0.001:0.1;
figure,plot(h,constr_mean(h),'-r',h,constr_cov(h),'-b')
%ylim([-2 2]);
figure,plot(h,Free_energy(h),'-k',h,Free_energy_cov(h),'-g')
ylim([-2 2]);


h_mean=fzero(constr_mean,0);
h_var=fzero(constr_cov,h_mean);

mean_seq=[tanh(M1*h_mean) tanh(M2*h_mean) tanh(h_mean)];


cov_12=tanh(M1*h_var)*tanh(M2*h_var)*exp(-Free_energy_cov(h_var)+Free_energy(h_mean))-tanh(M1*h_mean)*tanh(M2*h_mean);



%% two constraint
clear all; close all;
scale1=1;%[1 0.5 2 5 10];
shift1=[0 0 1 2];
%shift1=[0 0 1];
scale2=1;%[1 0.5 5 10 20];
shift2=[0 1 2 1];
%shift2=[0 1 2];
for m=1:length(shift1)
for l=1:length(scale1)
[count,weight1,weight2,mean_sequence,seq_corr,gamma1,gamma2]=sequence_mutation_Math_double_select(scale1(l),scale2(l),shift1(m),shift2(m));
data{l,m,1}=count;
data{l,m,2}=weight1;
data{l,m,3}=weight2;
data{l,m,4}=mean_sequence;
data{l,m,5}=seq_corr;
data{l,m,6}=gamma1;
data{l,m,7}=gamma2;
end
end


for m=1:length(shift1)
for l=1:length(scale1)
    weight1=data{l,m,2};
    weight2=data{l,m,3};
   
    mean_sequence=data{l,m,4};
    sigma=sqrt(1-mean_sequence.^2);
    
    f1=weight1.*sigma;
    f2=weight2.*sigma;
    
    seq_corr=data{l,m,5};
    gamma1=data{l,m,6};
    gamma2=data{l,m,7};
 
    C11=zeros(1,1);C12=zeros(1,1);C22=zeros(1,1);
    yy_s=zeros(1,1);
    s=0;
    for i=1:length(weight1)
        for j=i+1:length(weight1)
            if abs(seq_corr(i,j))>0.01;
           s=s+1;
           C11(s)=f1(i)*f1(j);
           C12(s)=f1(i)*f2(j);
           C22(s)=f2(i)*f2(j);
           yy_s(s)=seq_corr(i,j);%mean_sequence(i);%
            end
        end
    end

  
    R=[[sum(weight1.*weight1)  sum(weight1.*weight2)];[ sum(weight2.*weight1) sum(weight2.*weight2)]];
    inv_R=inv(R);
    
    A=[[sum(weight1.*weight1.*sigma)  sum(weight1.*weight2.*sigma)];[ sum(weight2.*weight1.*sigma) sum(weight2.*weight2.*sigma)]];
    inv_A=inv(A);
   
    plot_data{l,m,1}=C11*inv_A(1,1)+2*C12*inv_A(1,2)+C22*inv_A(2,2);
    plot_data{l,m,2}=yy_s;
    
    plot_data{l,m,3}=1-mean_sequence.^2;
    plot_data{l,m,4}=weight1*inv_R(1,1)*gamma1+weight1*inv_R(1,2)*gamma2+ weight2*inv_R(2,1)*gamma1 +weight2*inv_R(2,2)*gamma2 ;
    plot_data{l,m,5}=mean_sequence;
   
end
end

%save data_universal_covar_two_selection

clear f1 f2

m=1;

figureParameter
f1=plot(plot_data{1,m,1},plot_data{1,m,2},'*r',plot_data{2,m,1},plot_data{2,m,2},'xb');%,plot_data{3,m,1},plot_data{3,m,2},'ok',plot_data{4,m,1},plot_data{4,m,2},'sg',plot_data{5,m,1},plot_data{5,m,2},'+cyan');
ylim([-0.3 0]);
xlim([0 0.2]);
%a1=xlabel('$\sigma_l\sigma_{l''}\sum_{j,j''}\Delta_j^l\Delta_{j''}^{l''}A^{-1)_{jj''}$');
a2=ylabel('$C_{ll''}$');
%h1=legend('$\bar{\gamma}=0$','$\bar{\gamma}=0.5$','$\bar{\gamma}=1$','$\bar{\gamma}=1.5$','$\bar{\gamma}=2$');
%set(h1,'location','southeast')
fig_name='./figure/Corr_collapse.eps';
figurePostTreat


figureParameter
f1=plot(plot_data{1,m,4},plot_data{1,m,5},'*r',plot_data{2,m,4},plot_data{2,m,5},'xb',plot_data{3,m,4},plot_data{3,m,5},'ok',plot_data{4,m,4},plot_data{4,m,5},'sg',plot_data{5,m,4},plot_data{5,m,5},'+cyan');
ylim([0 1]);
xlim([0 1]);
a1=xlabel('$\sum_{jj''} \Delta_j^lR^{-1}_{jj''}\gamma_j{''}$');
a2=ylabel('$\langle \alpha_l\rangle$');
%h1=legend('$\bar{\gamma}=0$','$\bar{\gamma}=0.5$','$\bar{\gamma}=1$','$\bar{\gamma}=1.5$','$\bar{\gamma}=2$');
%set(h1,'location','southeast')
fig_name='./figure/mean_collapse.eps';
figurePostTreat



l=5;

clear f1 f2

figureParameter
f1=plot(plot_data{l,1,1},plot_data{l,1,2},'*r',plot_data{l,2,1},plot_data{l,2,2},'xb',plot_data{l,3,1},plot_data{l,3,2},'ok',plot_data{l,4,1},plot_data{l,4,2},'sg');
ylim([-0.4 0]);
xlim([0 0.3]);
a1=xlabel('$\sigma_l\sigma_{l''}\sum_{j,j''}\Delta_j^lA^{-1)_{jj''}\Delta_{j''}^{l''}$');
a2=ylabel('$C_{ll''}$');
h1=legend('$\bar{\gamma}_1=0,\bar{\gamma}_2=0\;\;$','$\bar{\gamma}_1=0,\bar{\gamma}_2=1$','$\bar{\gamma}_1=1,\bar{\gamma}_2=2$','$\bar{\gamma}_1=2,\bar{\gamma}_2=1$');
%set(h1,'location','southeast')
fig_name='./figure/Corr_collapse.eps';
figurePostTreat


figureParameter
f1=plot(plot_data{l,1,4},plot_data{l,1,5},'*r',plot_data{l,2,4},plot_data{l,2,5},'xb',plot_data{l,3,4},plot_data{l,3,5},'ok',plot_data{l,4,4},plot_data{l,4,5},'sg');
ylim([0 1]);
xlim([0 1]);
a1=xlabel('$\sum_{jj''} \Delta_j^lR^{-1}_{jj''}\gamma_{j''}$');
a2=ylabel('$\langle \alpha_l\rangle$');
%h1=legend('$\bar{\gamma}=0$','$\bar{\gamma}=0.5$','$\bar{\gamma}=1$','$\bar{\gamma}=1.5$','$\bar{\gamma}=2$');
%set(h1,'location','southeast')
fig_name='./figure/mean_collapse.eps';
figurePostTreat



%% one constraint
clear all; close all;
scale=5;
%scale=[1 5 10 20];
%shift=[0 0.5 1 1.5];
%shift=[0 0.25 0.5 1];
shift=[0 3.5 8];
for m=1:length(shift)
for l=1:length(scale)
[factor,count,weight,mean_sequence,seq_corr,cover_ll]=sequence_mutation_Math_single_select(scale(l),shift(m));
data{l,m,1}=factor;
data{l,m,2}=count;
data{l,m,3}=weight;
data{l,m,4}=mean_sequence;
data{l,m,5}=seq_corr;
data{l,m,6}=cover_ll;
end
end



%load data_universal_covar_more
%load data_universal_covar

for m=1:length(shift)
for l=1:length(scale)
    factor=data{l,m,1};
    weight=data{l,m,3};
    mean_sequence=data{l,m,4};
    seq_corr=data{l,m,5};
    cover_ll=data{l,m,6};
    %cover_ll=1-mean_sequence.^2;
    f=weight.^2.*cover_ll';
    tot_f=sum(f);
    xx_w=zeros(1,1);
    yy_s=zeros(1,1);
    zz_s=zeros(1,1);
    s=0;
    for i=1:length(weight)
        for j=i+1:length(weight)
            if abs(seq_corr(i,j))>0.0001
           s=s+1;
           %xx_w(s)=shift(m)*weight(i)/sqrt(sum(weight.^2));
           temp=f(i)*f(j)/(weight(i)*weight(j));%+(f(i)+f(j))*mean_sequence(i)*mean_sequence(j);
           xx_w(s)=temp/tot_f;
           yy_s(s)=seq_corr(i,j);
           zz_s(s)=(1-mean_sequence(i)^2)*(1-mean_sequence(j)^2);
            end
        end
    end

    plot_data{l,m,1}=xx_w;%/sum(weight.^2);  % prediction of covariance matrix
    plot_data{l,m,2}=yy_s;  % actual value of covariance matrix
    plot_data{l,m,3}=cover_ll;  % variance
    %plot_data{l,m,4}=shift(m)*weight/sqrt(sum(weight.^2)); % prediction of shift
    plot_data{l,m,4}=shift(m)*weight/sum(weight.^2);
    plot_data{l,m,5}=mean_sequence; % actual value of shift
    
    
end
end

figureParameter
f1=pcolor(1:count,1:count,seq_corr);
colorbar;
set(gca, 'clim', [-0.05 0.05]);
xlabel('Residue');
ylabel('Residue');
%set(gca, 'clim', [-1 1]);
%a1=xlabel('$x$');
%a2=ylabel('$y$');
fig_name='./figure/mode_correlation.eps';
figurePostTreat

% m=3;
% figure, plot3(plot_data{1,m,1},plot_data{1,m,3},plot_data{1,m,2},'*r',plot_data{2,m,1},plot_data{2,m,3},plot_data{2,m,2},'xb',plot_data{3,m,1},plot_data{3,m,3},plot_data{3,m,2},'ok',plot_data{4,m,1},plot_data{4,m,3},plot_data{4,m,2},'sg')
% xlabel('Delta1*Delta2');
% ylabel('Simga1*Simga2');

k=1;
figureParameter
f1=plot(-plot_data{k,1,1},plot_data{k,1,2},'ok',-plot_data{k,2,1},plot_data{k,2,2},'xb',-plot_data{k,3,1},plot_data{k,3,2},'*g',-plot_data{1,4,1},plot_data{1,4,2},'.r');%,plot_data{k,5,1},plot_data{k,5,2},'*cyan');
ylim([-0.06 0.06]);
%xlim([0 0.05]);
%axis tight
a1=xlabel('$\Delta_l\Delta_{l''}\sigma_l^2\sigma_{l''}^2/\sum_l \Delta_l^2\sigma_l^2$');
%a1=xlabel('$\Delta_l\Delta_{l''}/\sum_l \Delta_2^l$');
%a2=ylabel('$C_{ll''}$');
set(gca, 'YTick' , [-0.05 0 0.05 ]);
%h1=legend('$\gamma=0$','$\gamma=0.25$','$\gamma=0.5$','$\gamma=1$');
%set(h1,'location','northeast')
fig_name='./figure/coavr_collapse.eps';
figurePostTreat




figureParameter
f1=plot(plot_data{k,1,4},plot_data{k,1,5},'ok',plot_data{k,2,4},plot_data{k,2,5},'xb',plot_data{k,3,4},plot_data{k,3,5},'*g',plot_data{1,4,4},plot_data{1,4,5},'.r');%,plot_data{k,5,4},plot_data{k,5,5},'*cyan');
axis tight
ylim([0 1]);
xlim([-0.5 0.5]);
a1=xlabel('$\gamma \Delta_l/\sum_l \Delta_l^2$');
%a2=ylabel('$\langle S_l\rangle_*$');
h1=legend('$\delta\bar{E}^*=0$','$\delta\bar{E}^*=0.25$','$\delta\bar{E}^*=0.5$','$\delta\bar{E}^*=1$');
%set(h1,'location','southeast')
fig_name='./figure/mean_collapse.eps';
figurePostTreat





figureParameter
f1=plot(xx_w,yy_s,'*b');
xlim([0 0.01]);
ylim([-0.01 0]);
a1=xlabel('$\Delta_l\Delta_{l''}\sigma_l\sigma_{l''}/\sum_l \Delta_l^2\sigma_l$');
%a1=xlabel('$\Delta_l\Delta_{l''}/\sum_l \Delta_2^l$');
a2=ylabel('$C_{ll''}$');
%h1=legend('$\bar{\gamma}=0$','$\bar{\gamma}=0.5$','$\bar{\gamma}=1$','$\bar{\gamma}=1.5$','$\bar{\gamma}=2$');
%set(h1,'location','southwest')
fig_name='./figure/coavr_collapse.eps';
figurePostTreat


m=1;
figureParameter
f1=plot(tanh(xx_w),yy_s,'*r');
ylim([0 1]);
a1=xlabel('$\gamma \Delta_l/\sum_l \Delta_2^l$');
a2=ylabel('$\langle \alpha_l\rangle$');
%a1=xlabel('$\Delta_l\Delta_{l''}\sigma_l\sigma_{l''}/\sum_l \Delta_2^l\sigma_l$');
fig_name='./figure/var.eps';
figurePostTreat



m=1;
figureParameter
f1=plot(plot_data{1,m,1},plot_data{1,m,3},'*r',plot_data{2,m,1},plot_data{2,m,3},'xb',plot_data{3,m,1},plot_data{3,m,3},'ok',plot_data{4,m,1},plot_data{4,m,3},'sg');
%ylim([-1 1]);
a1=xlabel('$l$');
a2=ylabel('$\sigma_{l}$');
fig_name='./figure/var.eps';
figurePostTreat

m=1;
figureParameter
f1=plot(plot_data{1,m,1},plot_data{1,m,2},'*r',plot_data{2,m,1},plot_data{2,m,2},'xb',plot_data{3,m,1},plot_data{3,m,2},'ok',plot_data{4,m,1},plot_data{4,m,2},'sg');
ylim([-0.3 0]);
a1=xlabel('$\Delta_l\Delta_{l''}$');
a2=ylabel('$C_{ll''}$');
h1=legend('$\Delta^*=10$','$\Delta^*=50$','$\Delta^*=100$','$\Delta^*=200$');
set(h1,'location','northeast')
fig_name='./figure/coavr_collapse.eps';
figurePostTreat
% 
% k=2;
% figureParameter
% f1=plot(1:count,plot_data{k,1,3},'*r',1:count,plot_data{k,2,3},'xb',1:count,plot_data{k,3,3},'ok',1:count,plot_data{1,4,3},'sg',1:count,plot_data{1,5,3},'*cyan');
% %ylim([-1 1]);
% a1=xlabel('$l$');
% a2=ylabel('$\sigma_l$');
% h1=legend('$\bar{\gamma}=0$','$\bar{\gamma}=0.5$','$\bar{\gamma}=1$','$\bar{\gamma}=1.5$','$\bar{\gamma}=2$');
% set(h1,'location','southwest')
% fig_name='./figure/coavr_collapse.eps';
% figurePostTreat



% k=2;
% figureParameter
% f1=plot(shift,plot_data{k,1,3},'*r',1:count,plot_data{k,2,3},'xb',1:count,plot_data{k,3,3},'ok',1:count,plot_data{1,4,3},'sg',1:count,plot_data{1,5,3},'*cyan');
% %ylim([-1 1]);
% a1=xlabel('$$');
% a2=ylabel('$\sigma_l$');
% h1=legend('$\bar{\gamma}=0$','$\bar{\gamma}=0.5$','$\bar{\gamma}=1$','$\bar{\gamma}=1.5$','$\bar{\gamma}=2$');
% set(h1,'location','southwest')
% fig_name='./figure/coavr_collapse.eps';
% figurePostTreat

k=2;
figureParameter
f1=plot(plot_data{k,1,4},plot_data{k,1,5},'*r',plot_data{k,2,4},plot_data{k,2,5},'xb',plot_data{k,3,4},plot_data{k,3,5},'ok',plot_data{1,4,4},plot_data{1,4,5},'sg',plot_data{k,5,4},plot_data{k,5,5},'*cyan');
ylim([0 1]);
a1=xlabel('$\gamma \Delta_l/\sum_l \Delta_2^l$');
a2=ylabel('$\langle \alpha_l\rangle$');
h1=legend('$\bar{\gamma}=0$','$\bar{\gamma}=0.5$','$\bar{\gamma}=1$','$\bar{\gamma}=1.5$','$\bar{\gamma}=2$');
set(h1,'location','southeast')
fig_name='./figure/mean_collapse.eps';
figurePostTreat





figure, plot(plot_data{1,1,1},plot_data{1,1,2},'*r',plot_data{2,1,1},plot_data{2,1,2},'xb',plot_data{3,1,1},plot_data{3,1,2},'ok',plot_data{4,1,1},plot_data{4,1,2},'sg')

figure, plot(plot_data{1,1,1},plot_data{1,1,2},'*r',plot_data{2,1,1},plot_data{2,1,2},'xb',plot_data{3,1,1},plot_data{3,1,2},'ok',plot_data{4,1,1},plot_data{4,1,2},'sg')

%%
[factor,count,weight,mean_sequence,seq_corr]=sequence_mutation_Math(10,2);

select_data=seq_corr;

[NormVector_corr,orderEigValue_corr]=orderedEigSystem(select_data,0);

Delta=weight;
Delta=Delta./sqrt(sum(Delta.^2));


figureParameter
bar(1:count,mean_sequence,'r');
xlim([0 count]);
ylim([-1 1]);
a1=xlabel('$l$');
a2=ylabel('$\langle \alpha_l\rangle$');
fig_name='./figure/mean.eps';
figurePostTreat


figureParameter
bar(1:count,[abs(NormVector_corr(:,1)),Delta]);
xlim([0 count]);
%set(gca,'XTICKlabel',index2);
h1=legend('$\nu_1^l$','$\Delta_l$');%-\Delta_6^l\;$');
set(h1,'location','north')
%ylim([0 0.5]);
fig_name='./figure/corr_mode2.eps';
figurePostTreat

%%
I=sqrt(-1);
f=@(x) exp(I*x+100*log(cosh(I*x)));  
X=-1:0.02:1;
Y=0:0.02:3;

figure,
plot(X,real(f(X)),'-r',X,real(f(X+I*0.1)),'-b',X,real(f(X+I*0.5)),'-g',X,real(f(X+I*1)),'-k')

figure,
plot(X,real(f(X-I*10)),'-r')


Z0=X+I*Y;
Z1=Z0'*Z0;

figureParameter
f1=surf(X,Y,real(f(Z1)));
a1=xlabel('$x$');
a2=ylabel('$y$');
fig_name='./figure/distribution.eps';
figurePostTreat

figureParameter
f1=contour(X,Y,real(f(Z1)));
a1=xlabel('$x$');
a2=ylabel('$y$');
fig_name='./figure/distribution.eps';
figurePostTreat

%% check eigenvalues of the covariance matrix

% l controls the heterogeneity of Delta, or strength of the hotspots 
% m controls the shift gamma.  Here, m=1,  no shift. 

l=1;% only change l, which range from 1 to 5
m=1; 

weight1=data{l,m,2}; % Delta_1
weight2=data{l,m,3}; % Delta_2

mean_sequence=data{l,m,4};
sigma=sqrt(1-mean_sequence.^2);

seq_corr=data{l,m,5};

% weight1=Delta_lambda3';
% weight2=Delta_lambda4';

R=[[sum(weight1.*weight1), sum(weight1.*weight2)]; [sum(weight1.*weight2), sum(weight2.*weight2)]];
    
[NormVector_R,orderEigValue_R]=orderedEigSystem(R,0);

NormVector_R(:,1)=NormVector_R(:,1)*1/sqrt(orderEigValue_R(1,1));
NormVector_R(:,2)=NormVector_R(:,2)*1/sqrt(orderEigValue_R(2,1));

diagonal_form=NormVector_R'*R*NormVector_R;
D=NormVector_R';

new_weight1=D(1,1)*weight1+D(1,2)*weight2;
new_weight2=D(2,1)*weight1+D(2,2)*weight2;

%new_weight1=new_weight1./sqrt(sum(new_weight1.*new_weight1));
%new_weight2=new_weight2./sqrt(sum(new_weight2.*new_weight2));

%coupling=sum(weight1.*weight2)/(sqrt(sum(weight1.^2))*sqrt(sum(weight2.^2)));


select_data=seq_corr;

[NormVector_corr,orderEigValue_corr]=orderedEigSystem(select_data,0);


%Delta_lambda_plus=weight1./sqrt(sum(weight1.^2))+weight2./sqrt(sum(weight2.^2));
Delta_lambda_plus=weight1+weight2;
Delta_lambda_plus=Delta_lambda_plus./sqrt(sum(Delta_lambda_plus.*Delta_lambda_plus));

%Delta_lambda_mus=weight1./sqrt(sum(weight1.^2))-weight2./sqrt(sum(weight2.^2));
%Delta_lambda_mus=weight2;%weight1-weight2;
%Delta_lambda_mus=Delta_lambda_mus./sqrt(sum(Delta_lambda_mus.*Delta_lambda_mus));

clear f1 f2 a1 a2
figureParameter
bar(1:count,[NormVector_corr(:,1),-new_weight2]);
xlim([0 count]);
%set(gca,'XTICKlabel',index2);
h1=legend('$\nu_1^l$','$\sum_{j}D_{1j}\Delta_j^l\;$');%+\Delta_2^l$');%-\Delta_6^l\;$');
set(h1,'location','north')
%ylim([0 0.5]);
fig_name='./figure/corr_mode1.eps';
figurePostTreat

figureParameter
bar(1:count,[NormVector_corr(:,2),new_weight1]);
%set(gca,'XTICKlabel',index2);
%bar(1:count,[abs(NormVector_corr(:,1)),abs(Delta_lambda_plus)]);
xlim([0 count]);
h1=legend('$\nu_2^l$','$\sum_jD_{2j}\Delta_j^l\;$');%$\Delta_1^l-\Delta_2^l$');%-\Delta_6^l\;$');
set(h1,'location','north')
%ylim([0 0.5]);
fig_name='./figure/corr_mode2.eps';
figurePostTreat


%%
l=5;m=1;

weight=data{l,m,3};
mean_sequence=data{l,m,4};
seq_corr=data{l,m,5};


coupling=sum(weight1.*weight2)/(sqrt(sum(weight1.^2))*sqrt(sum(weight2.^2)));

select_data=seq_corr;

[NormVector_corr,orderEigValue_corr]=orderedEigSystem(select_data,0);

orderEigValue_corr(1,:)

%% finite size effect
close all
clear all;

count=100;
%M=floor(100*1.2.^(1:N));
kappa=[100 120 150 180 200 220 250 280 300 350 400 500 600 700 800 900 1000];%
%kappa=[300 500 800];
N=length(kappa);
ratio=zeros(N,count);
M_new=zeros(N,1);
for j=1:N
[factor,count,weight,mean_sequence,seq_corr,cover_ll,M_new(j)]=sequence_mutation_Math(1,0,100000,kappa(j));
select_data=seq_corr;
[NormVector_corr,orderEigValue_corr]=orderedEigSystem(select_data,0);
scale_K=sqrt(sum(weight.^2));
weight=weight/scale_K;
ratio(j,:)=(abs(weight)'*abs(NormVector_corr))';
end

clear h1;

figureParameter
f1=plot(scale_K./kappa,ratio(:,1),'*r');
%set(gca,'XTICK',[20 40 60]);
%set(gca,'XTICKlabel',index2);
a1=xlabel('$\kappa/\sqrt{\sum_l\Delta_l^2}$');
a2=ylabel('Similarlity');
%bar(1:count,[abs(NormVector_corr(:,1)),abs(Delta_lambda_plus)]);
%xlim([0 40]);
ylim([0 1]);
fig_name='./figure/corr_mode2.eps';
figurePostTreat


figureParameter
bar(1:count,ratio1,'r');
%set(gca,'XTICKlabel',index2);
a1=xlabel('Mode: $j$');
a2=ylabel('Similarlity');
%bar(1:count,[abs(NormVector_corr(:,1)),abs(Delta_lambda_plus)]);
xlim([0 100]);
ylim([0 1]);
fig_name='./figure/corr_mode2.eps';
figurePostTreat


%% dealing with real data
% mapping


figure,hist(efa1(:,7));
xlim([0 20]);


AA=randn(1500,1500);

[c,d]=orderedEigSystem(AA,0);

%% using our covariance matrix to test eigenvectors
count=100;
sigma=ones(count,1);
covar=zeros(count,count);
Delta=ones(count,1);
max_N=5;
Delta(1:max_N)=[2 3 4 5 20];
tot=sum(Delta.^2.*sigma.^2);
for i=1:count
    for j=1:count
       if i~=j
          covar(i,j)=-Delta(i)*Delta(j)*sigma(i)^2*sigma(j)^2/tot;
        else
            covar(i,j)=sigma(i)^2;
        end
    end
end

esti_eig=zeros(count,1);

for i=1:count
    esti_eig(i)=sigma(i)^2+(1-Delta(i)^2/(sum(Delta.^2)-Delta(i)^2))*sigma(i)^2*Delta(i)^2/tot;
end
esti_eig=sort(esti_eig);


select_data=covar;

[NormVector_corr,orderEigValue_corr]=orderedEigSystem(select_data,0);


figureParameter
f1=plot(1:count,orderEigValue_corr(:,1),'or',1:count,esti_eig,'*b');
a1=xlabel('Mode: $k$');
a2=ylabel('Eigenvalue: $\lambda_k^c$');
xlim([0 count]);
fig_name='./figure/corr_eig.eps';
figurePostTreat



Delta_weight=Delta./(sqrt(sum(Delta.*Delta)));

figureParameter
bar(1:count,[NormVector_corr(:,1),Delta_weight]);
%set(gca,'XTICKlabel',index2);
h1=legend('$\nu_{1}^l$','$\Delta_3^l$');%-\Delta_6^l\;$');
set(h1,'location','north')
xlim([0 count]);
fig_name='./figure/corr_mode2.eps';
figurePostTreat

Delta1=Delta;
Delta_weight=Delta1./(sqrt(sum(Delta1.*Delta1)));

big_v=NormVector_corr(:,1);
%big_v(5)=0;
big_v=big_v./sum(big_v.*big_v);

figureParameter
bar(1:count,[big_v./Delta_weight]);
%set(gca,'XTICKlabel',index2);
h1=legend('$\nu_{100}^l$','$\Delta_3^l$');%-\Delta_6^l\;$');
set(h1,'location','south')
xlim([0 count]);
fig_name='./figure/corr_mode2.eps';
figurePostTreat


%% inverse matrix 
N0=20;
Delta=ones(N0,1);
Delta(1:5)=3*[2 3 4 5 6];%10*(1:4);
Sigma=0.5*ones(N0,1)+0.1*randn(N0,1);
S=zeros(N0,N0);
d=Delta.*Sigma.^2/sqrt(sum(Delta.^2.*Sigma.^2));

for j=1:N0
    S(j,j)=Sigma(j)^2+d(j)^2;
end


C=S-d*d';


max_N=20;
figureParameter
f1=pcolor(1:max_N,1:max_N,C(1:max_N,1:max_N));
colorbar;
set(gca, 'clim', [-0.05 0.05]);
%set(gca, 'clim', [-1 1]);
%a1=xlabel('$x$');
%a2=ylabel('$y$');
fig_name='./figure/mode_correlation.eps';
figurePostTreat


% [NormVector_corr,orderEigValue_corr]=orderedEigSystem(C,0);
% Delta_1=Delta./sqrt(sum(Delta.^2));
% figure,bar(1:N0,[Delta_1,NormVector_corr(:,1)]);
% figure,plot(1:N0,orderEigValue_corr(:,1),'*b');

% analytical result
C_inv=inv(S)+inv(S)*d*d'*inv(S)/(1-d'*inv(S)*d);


epsilon=1-sum(d.^2./(d.^2+Sigma.^2));
v_bar=d./(Sigma.^2+d.^2);
Cinv2=zeros(N0,N0);

for j=1:N0
    for k=1:N0
        if j~=k
      Cinv2(j,k)=v_bar(j)*v_bar(k)/epsilon;
        else
            Cinv2(j,j)=1/(Sigma(j)^2+d(j)^2)+v_bar(j)^2/epsilon;
        
        end
    end
end

max_N=20;
figureParameter
f1=pcolor(1:max_N,1:max_N,Cinv2(1:max_N,1:max_N));
colorbar;
set(gca, 'clim', [-10 10]);
%set(gca, 'clim', [-1 1]);
%a1=xlabel('$x$');
%a2=ylabel('$y$');
fig_name='./figure/mode_correlation.eps';
figurePostTreat

max_N=20;
figureParameter
bar(1:max_N,[Cinv2(1:max_N,1)./Delta(1),Cinv2(1:max_N,2)./Delta(2),Cinv2(1:max_N,5)./Delta(5),Cinv2(1:max_N,10)./Delta(10)]);
xlim([0 20]);
ylim([0,0.5]);
h1=legend('$C^{-1}_{1n}/\Delta_1$','$C^{-1}_{2n}/\Delta_2$','$C^{-1}_{5n}/\Delta_5$','$C^{-1}_{10n}/\Delta_{10}$');
fig_name='./figure/Cij.eps';
figurePostTreat

max_N=20;
figureParameter
bar(1:max_N,[C(1:max_N,1)./Delta(1),C(1:max_N,2)./Delta(2),C(1:max_N,5)./Delta(5),C(1:max_N,10)./Delta(10)]);
ylim([-0.02,0.02]);
xlim([0 20]);
h1=legend('$C_{1n}/\Delta_1$','$C_{2n}/\Delta_2$','$C_{5n}/\Delta_5$','$C_{10n}/\Delta_{10}$');
fig_name='./figure/Cij1.eps';
figurePostTreat


%% WW domain

[Prob_orig,X]=hist(latent,20);
Prob_orig=Prob_orig/sum(Prob_orig);
[Prob_sp,X1]=hist(latent_sp,20);
Prob_sp=Prob_sp/sum(Prob_sp);

figureParameter
f1=plot(X,Prob_orig,'-*k',X1,Prob_sp,'-or');
%f1=plot(X,Prob,'-*b');
%b=bar(X,Prob);
%b.CData(1,:) = [0 0.8 0.8];
xlim([0 3  ]);
a1=xlabel('$\sum_l \beta_l\nu_{end}(l)$');
ylim([0 0.1]);
a2=ylabel('Distribution');
h1=legend('original','scrambled$\;\;$');
set(h1,'location','northeast');
fig_name='./figure/distribution.eps';
figurePostTreat


figure,hist(folds_cc);
hist(folds_n);
[ele_CC,cent_CC]=hist(ScoresCCCov,10);
ele_CC=ele_CC./sum(ele_CC);
[ele_IC,cent_IC]=hist(ScoresICCov,10);
ele_IC=ele_IC/sum(ele_IC);

[ele_N,cent_N]=hist(ScoresNCov,10);
ele_N=ele_N/sum(ele_N);

X_tick=
[ele_PC,cent_PC]=hist(ScoresPCov,10);
ele_PC=ele_PC/sum(ele_PC);

[ele_Scr,cent_Scr]=hist(ScoresScrCov,10);
ele_Scr=ele_Scr/sum(ele_Scr);



figure,plot(cent_CC,ele_CC,'*r',cent_IC,ele_IC,'*b',cent_N,ele_N,'*g',cent_PC,ele_PC,'*k',cent_Scr,ele_Scr,'or')



%% inverse matrix: eigenvalue
L=5;
amino_N=21;
count=L*amino_N;

Delta1=randn(1,amino_N*L);
Delta2=randn(1,amino_N*L);

Delta1(3)=3;
Delta1(amino_N+5)=6;
Delta1(2*amino_N+8)=3;
Delta1(3*amino_N+10)=-5;
Delta1(4*amino_N+6)=4;

%% second constraint
Delta2(8)=4;
Delta2(amino_N+10)=6;
Delta2(2*amino_N+6)=3;
Delta2(3*amino_N+20)=4;
Delta2(4*amino_N+16)=2;



% 
% count=10;
% Delta1=randn(1,count);
% Delta1(2)=4; Delta1(6)=3;Delta1(9)=5;
% 
% Delta2=randn(1,count);
% Delta2(1)=4; Delta2(4)=3;Delta2(7)=8;

G=[sum(Delta1.*Delta1) sum(Delta1.*Delta2); sum(Delta2.*Delta1) sum(Delta2.*Delta2)];

c0=1;
Cinv=Delta1'*Delta1+c0*Delta2'*Delta2;

for k=1:count
   Cinv(k,k)=0;
end

figureParameter
f1=image(1:count,1:count,5*abs(Cinv));
colorbar;
%set(gca, 'clim', [-0.25 0.25]);
%set(gca, 'clim', [-1 1]);
%a1=xlabel('$x$');
%a2=ylabel('$y$');
fig_name='./figure/mode_correlation.jpg';
set(gca,'YDir','normal')
print(fig_name,'-r500','-djpeg');


select_data=Cinv;

[NormVector_corr,orderEigValue_corr]=orderedEigSystem(select_data,0);


figureParameter
f1=plot(1:count,orderEigValue_corr(:,1),'or');
a1=xlabel('Mode: $k$');
a2=ylabel('Eigenvalue: $\lambda_k^c$');
xlim([1 count]);
%ylim([-1000 1000]);
fig_name='./figure/corr_eig.eps';
figurePostTreat




weight1=Delta1';
Delta_weight1=weight1./(sqrt(sum(weight1.*weight1)));

weight2=Delta2';
Delta_weight2=weight2./(sqrt(sum(weight2.*weight2)));


figureParameter
bar(1:count,[NormVector_corr(:,end-1),Delta_weight1]);
%set(gca,'XTICKlabel',index2);
h1=legend('$\nu_{end-1}^l$','$\tilde{\Delta}_l^{(1)}$');%-\Delta_6^l\;$');
set(h1,'location','north')
xlim([0 count+1]);
fig_name='./figure/corr_mode1.eps';
figurePostTreat

figureParameter
bar(1:count,[NormVector_corr(:,end),Delta_weight1]);
%set(gca,'XTICKlabel',index2);
h1=legend('$\nu_{end}^l$','$\tilde{\Delta}_l^{(1)}$');%-\Delta_6^l\;$');
set(h1,'location','north')
ylim([-0.6 1]);
xlim([0 count+2]);
fig_name='./figure/corr_mode2.eps';
figurePostTreat

%%
select_N=[ 150 160 180 190 200 210 220 230 240 250 260 270 280 290 300 310 320 330 340 350 360 370 380 390 400 410];
ratio_end=zeros(1,length(select_N));
M_new=zeros(1,length(select_N));
for j=1:length(select_N);
[M_new(j),ratio_end(j)]=generate_artificial_protein_sequence_inverse_matrix(select_N(j));
end

figureParameter
f1=plot(M_new,ratio_end,'*');
a1=xlabel('data number');
a2=ylabel('Recovery');
xlim([150 450]);
fig_name='./figure/corr_mode2.eps';
figurePostTreat

%%  DCA, pair identification
FW1=FW+FW';
count=79;
figureParameter
f1=image(1:count,1:count,8*FW1);
colorbar;
fig_name='./figure/corr1.jpg';
set(gca,'YDir','normal')
print(fig_name,'-r500','-djpeg');

%% matrix compression
inv_compressed=zeros(L,L);
for j=1:length(index_k)
    temp1=index_k(j)+1:index_k(j)+conserved_num(j);
   for k=1:length(index_k)
    temp2=index_k(k)+1:index_k(k)+conserved_num(k);   
    tot_N=conserved_num(j)*conserved_num(k);
    data=inv_corr_off_diag(temp1,temp2);
    inv_compressed(j,k)=sqrt(sum(sum(data.^2))/tot_N);
       
   end
end

figureParameter
f1=image(1:L,1:L,0.5*inv_compressed);
colorbar;
fig_name='./figure/corr1.jpg';
set(gca,'YDir','normal')
print(fig_name,'-r500','-djpeg');



%% analyze data from Anne-Flo
count=79*21;
L=79;
data=zeros(count,count);
for j=1:L
      temp1=(j-1)*21+1:j*21;
    for k=1:L
      temp2=(k-1)*21+1:k*21;
      data(temp1,temp2)=W(j,k,:,:);  
    end
end

data1=data+data';
for j=1:count
    data1(j,j)=0;
end

figureParameter
f1=image(1:count,1:count,100*data1);
colorbar;
fig_name='./figure/corr1.jpg';
set(gca,'YDir','normal')
print(fig_name,'-r500','-djpeg');

%% compress 

L=79;
data=zeros(L,L);
for j=1:L
    for k=1:L
        tot_N=21*21;
        data0=W(j,k,:,:);
        data(j,k)=sqrt(sum(sum(data0.^2))/tot_N);
    end
end

data1=data+data';
for j=1:L
    data1(j,j)=0;
end

figureParameter
f1=image(1:L,1:L,200*data1);
colorbar;
fig_name='./figure/corr1.jpg';
set(gca,'YDir','normal')
print(fig_name,'-r500','-djpeg');

%%
close all;
bias=0;
[predict_mean_beta_0,mean_sequence_0,predict_covar_0,covar_1d_0]=sequence_mutation_Math(bias);

bias=2;
[predict_mean_beta_2,mean_sequence_2,predict_covar_2,covar_1d_2]=sequence_mutation_Math(bias);

bias=5;
[predict_mean_beta_5,mean_sequence_5,predict_covar_5,covar_1d_5]=sequence_mutation_Math(bias);

bias=10;
[predict_mean_beta_10,mean_sequence_10,predict_covar_10,covar_1d_10]=sequence_mutation_Math(bias);

bias=20;
[predict_mean_beta_20,mean_sequence_20,predict_covar_20,covar_1d_20]=sequence_mutation_Math(bias);


figureParameter
f1=plot(predict_mean_beta_0,mean_sequence_0,'*r',predict_mean_beta_2,mean_sequence_2,'xb',predict_mean_beta_5,mean_sequence_5,'ok',predict_mean_beta_10,mean_sequence_10,'sg',predict_mean_beta_20,mean_sequence_20,'*cyan');
ylim([-0.5 0.5]);
xlim([-0.5 0.5]);
a1=xlabel('$\gamma \Delta_l/\sum_l \Delta_2^l$');
a2=ylabel('$\langle \alpha_l\rangle$');
h1=legend('$\bar{\gamma}=0$','$\bar{\gamma}=0.5$','$\bar{\gamma}=1$','$\bar{\gamma}=1.5$','$\bar{\gamma}=2$');
set(h1,'location','southeast')
fig_name='./figure/mean_collapse.eps';
figurePostTreat




figureParameter
f1=plot(predict_covar_0,covar_1d_0,'*r',predict_covar_2,covar_1d_2,'xb',predict_covar_5,covar_1d_5,'ok');%predict_covar_10,covar_1d_10,'sg',predict_covar_20,covar_1d_20,'*cyan');
ylim([-0.1 0]);
a1=xlabel('$\Delta_l\Delta_{l''}\sigma_l^2\sigma_{l''}^2/\sum_l \Delta_l^2\sigma_l^2$');
%a1=xlabel('$\Delta_l\Delta_{l''}/\sum_l \Delta_2^l$');
a2=ylabel('$C_{ll''}$');
h1=legend('$\bar{\gamma}=0$','$\bar{\gamma}=0.5$','$\bar{\gamma}=1$','$\bar{\gamma}=1.5$','$\bar{\gamma}=2$');
set(h1,'location','southwest')
fig_name='./figure/coavr_collapse.eps';
figurePostTreat

%%
x=0:0.01:1;

q0=0.5;
entropy=@(x) x.*log(x./q0)+(1-x).*log((1-x)./(1-q0));
phi_conserv=@(x) -log((1-x)./x)+log((1-q0)/q0);
var=x.*(1-x); 

figureParameter
f1=plot(x,entropy(x),'b');
a1=xlabel('$\langle S^*\rangle$');
%a2=ylabel('KL divergence');
fig_name='./figure/KLdivergence.eps';
figurePostTreat

figureParameter
left_color = [0 0 255]/255;
right_color = [255 0 0 ]/255;
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
f1=plot(x,phi_conserv(x));

yyaxis right
f2=plot(x,var);
a1=xlabel('$\langle S^*\rangle$');
fig_name='./figure/phi_var.eps';
figurePostTreat


f1=plot(x,phi_conserv(x),'b');
a1=xlabel('$\langle S^*\rangle$');
%a2=ylabel('$\phi$');
h1=legend('$\phi_l$','$\sigma_l^2$');
fig_name='./figure/phi.eps';
figurePostTreat

%%
figureParameter
%fig = figure;
red = [255 0 0]/255;
dark_green = [0 100 0]/255;
set(fig,'defaultAxesColorOrder',[red; dark_green]);

x=1:10;
y = x.^2;
yyaxis left
bar(x,y);

z = 10./x;
yyaxis right
f1=plot(x,z,'*');
fig_name='./figure/test.eps';
figurePostTreat


%% paper figure
clear all;close all;
% sample1=0:0.5:2;
% L1=length(sample1);
% bias(1:length(sample1))=sample1;
% sample2=3:10:500;
% L2=length(sample2);
% bias(L1:L1+L2-1)=sample2;
sample1=0:0.02:0.1;
L1=length(sample1);
bias(1:length(sample1))=sample1;
sample2=0.2:0.1:4;
L2=length(sample2);
bias(L1:L1+L2-1)=sample2;

result=zeros(length(bias),5);
for j=1:length(bias)
%[ratio_us,ratio_rama,ratio_inv]=sequence_mutation_Math(bias(j)) ;
 [ratio_us,ratio_rama,ratio_inv,ratio_PCA,ratio_rama_1,ratio_Cocco]=sequence_mutation_Math(bias(j));%,Delta_0) %(scale,shift,M,kappa)
result(j,1)=ratio_us(end);
result(j,2)=ratio_rama(1);
result(j,3)=ratio_inv(1);
result(j,4)=ratio_Cocco(end);
result(j,5)=ratio_PCA(end);
end

%save synthetic_sector_size_100_data_size_50000

% save paper_protein_conformation_Delta
%save paper_figure_Comparison_Delta_first_100;

clear all
load synthetic_sector_size_100_data_size_500000

figureParameter
f1=plot(bias,result(:,3),'-r',bias,result(:,5),'-.k',bias,result(:,2),'-.g',bias,result(:,4),'-b');
a1=xlabel('$\gamma$');
h1=legend('ICOD','PCA','SCA','Cocco');
set(h1,'location','southwest');
ylim([0 1]);
%xlim([0 40]);
fig_name='./figure/performance_comp.eps';
figurePostTreat


figureParameter
f1=plot(bias,result(:,3),'+r',bias,result(:,4),'sb');
a1=xlabel('$\delta E^*/Std(\delta E)$');
h1=legend('ICOD','Cocco');
set(h1,'location','southwest');
ylim([0 1]);
%xlim([0 40]);
fig_name='./figure/performance_comp_2.eps';
figurePostTreat


%% figure for papers: combine several different cases

result_1=result;

result_2=result;

result_100=result;

result_aver=(result_1+result_2+result_3+result_4+result_5+result_6+result_7+result_8+result_9+result_10+result_100)/11;

bias_2=bias;
result_1=result;  % Delta_first_1_1
bias_1=bias;
result_3=result;
bias_3=bias;
result_4=result;
bias_4=bias;

result_protein=result;

figureParameter
f1=plot(bias,result_1(:,2),'--b',bias,result_protein(:,2),'--r',bias,result_9(:,2),'--k',bias,result_1(:,3),'-b',bias(1:2:end),result_protein(1:2:end,3),'-r',bias,result_9(:,3),'-k',bias,result_aver,'-cyan');
a1=xlabel('$\delta E^*/Std(\delta E)$');
%h1=legend('Model 1','Model 1','Model 1','Model 1','Model 1','Model 1');
%set(h1,'location','south');
%h1=legend('Rama','Inv');
%set(h1,'location','south');
xlim([0 4]);
ylim([0 1]);
fig_name='./figure/performance_comp.eps';
figurePostTreat
