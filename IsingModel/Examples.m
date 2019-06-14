%% use the existing mutation vector
load ./Data/Delta_for_conformational_change

bias=0.1; %relative_selection_bias
Delta_0=Delta_lambda; % single-mutation effect
M=50000; % initial_random_sequence_number


% method, two-entry variable [a,b], how to perform selection: 
% .  [a=1,b=1], reweight each sequence using the quadratic free energy,x^2
% .  [a=1,b=0], reweight each sequence using the quartic free energy, x^4
%    [a=0,b=1], Window-selection to get some sequences, then each sequence has the same weight
%    [a=0,b=0], Threshold-selection to get some sequences, then each one has the same weight      
method=[1,1]; %reweight each sequence with the quadratic free energy
plot_fig=1;


sequence_mutation_Math(bias,Delta_0,M,method,plot_fig);

% or, add an option
selection_strength=1; 
sequence_mutation_Math(bias,Delta_0,M,method,plot_fig,selection_strength);

%% inside the function "sequence_mutation_Math", there is an section on initialization, where you choose different parameters

%%%%%% scenario=1;  % 0, for -1,1;  1 for 0, 1;  others
%%%%%% remove_conservation=0; % no remove by default

%% you can also synthesize your own Delta
%% Here is an example of the average performance for a given sector size

clear all;close all;

sample1=0:0.02:0.1;
L1=length(sample1);
bias(1:length(sample1))=sample1;
sample2=0.2:0.1:4;
L2=length(sample2);
bias(L1:L1+L2-1)=sample2;


method=[1,1]; %reweight each sequence with the quadratic free energy
plot_fig=0;

%bias=1:0.5:4;

Average_run=1; % you can change it too 100, which takes longer!
sector_size=10;

result_inv=zeros(length(bias),Average_run);
result_us=zeros(length(bias),Average_run);
result_rama=zeros(length(bias),Average_run);
result_rama_1=zeros(length(bias),Average_run);
result_Cocco=zeros(length(bias),Average_run);
result_PCA=zeros(length(bias),Average_run);
M=50000;

for j=1:length(bias)
    for k=1:Average_run
    
    Delta_0=1*randn(100,1);
    %factor=scale*floor(20*(randn(95,1)+0.5))+1;
    %factor=(4:1:10);
    
    %sector_size=k;
    Delta_0(1:sector_size)=20*randn(sector_size,1);

    
    [ratio_us,ratio_rama,ratio_inv,ratio_PCA,ratio_rama_1,ratio_Cocco]=sequence_mutation_Math(bias(j),Delta_0,M,method,plot_fig) ;
    result_inv(j,k)=ratio_inv(1);
    result_us(j,k)=ratio_us(end);
    result_rama(j,k)=ratio_rama(1);
    result_rama_1(j,k)=ratio_rama_1(1);
    result_PCA(j,k)=ratio_PCA(end);
    result_Cocco(j,k)=ratio_Cocco(end);
    end

end
save ./Data/20190228_paper_figure_Comparison_Single_Delta_Average_over_100_SectorSizes_M500000;



%% plot figures

clear all
load ./Data/20190228_paper_figure_Comparison_Single_Delta_Average_over_100_SectorSizes_M500000
% clear all;
% %load paper_figure_Comparison_Single_Delta_Average_sector_size_50_sequence_500000;
% load paper_figure_Comparison_Single_Delta_Average_elastic_network_Delta;

result_rama_aver=mean(result_rama(:,1:Average_run),2);
result_rama_1_aver=mean(result_rama_1(:,1:Average_run),2);
result_inv_aver=mean(result_inv(:,1:Average_run),2);
result_us_aver=mean(result_us(:,1:Average_run),2);
result_Cocco_aver=mean(result_Cocco(:,1:Average_run),2);
result_PCA_aver=mean(result_PCA(:,1:Average_run),2);


figureParameter
f1=plot(bias,result_inv_aver,'-r',bias,result_PCA_aver,'-.k',bias,result_rama_aver,'-.g',bias,result_Cocco_aver,'-b');
a1=xlabel('Relative bias');
%h1=legend('ICOD','PCA','SCA','Cocco');
%set(h1,'location','south');
%h1=legend('Rama','Inv');
%set(h1,'location','southwest');
xlim([0 4]);
ylim([0 1]);
fig_name='./figure/performance_comp.eps';
figurePostTreat

%% for two-sector case. This is less well organized.  Please look into the file
%% The two constraints Delta1 and Delta2 are generated inside the function
%% we use ICA to sepaerate the two constraints

bias1=1;
bias2=-1;
[ratio_us_1,ratio_us_2,ratio_rama_1,ratio_rama_2,ratio_inv_1,ratio_inv_2]=sequence_mutation_Math_double_select_new(bias1,bias2);
