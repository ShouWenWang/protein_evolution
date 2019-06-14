function [SCA_Exp_Corr,Conserv_Exp_Corr,Conserv_SCA_Corr,correct_hits_N20]=dealing_Protein_data_SCA(pseudo_ratio,cutoff_gap,using_PDZ_data,plot_fig)


%%%%%%%%%%%%%%%%%
% last updated, 2019-2-26
% This program is used to generate Fig.S18 in our paper
% Method of SCA is employed. The key difference from ICOD: we compute
% seq_corr without using a reference state. 

% Here, we reweight each state, and then apply Forbenius norm to compress the matrix

%%%%%%%%%%%%%%%%%%%


% clear all; 
% close all;

%% Initialization
%using_PDZ_data=1;
%cutoff_gap=0.3; % throw away the gap state (default, 0.3)
%pseudo_ratio=0.018; % maximum: 1
using_square_root=1; %use the square root of the SCA vector for prediction, consistent with the rest of the paper

%plot_fig=0;
seq_reweighting=0;   % 1, apply reweighting to each sequence; 0, do not apply. In this real data, we do not apply reweighting
similarity_threshold=0.9; % add more weight to sequences that are small than 90% similar to other sequences
save_all_results=0;

%%
addpath('./functions/');

if ~exist('figure','dir')
    mkdir figure;
end

if ~exist('Data','dir')
    mkdir hh
end

%% laod data
if using_PDZ_data
    %rand_sample=cumsum(floor(10*rand(2000,1))+1);
    %input_seq_data=efa1(rand_sample,:);
    load ./Data/PDZal;
    input_seq_data=efa1;


    load ./Data/McLaughlin_Ranganathan_data.mat;

    starting_index=312;  % in reference to the real protein index
    amino_N=21;
    referenced_mutation_effect=single_mutagenesis_cript.data;
    referenced_mutation_effect=referenced_mutation_effect(:,3:81); %20*79 is the dimension
else
    disp("Running with syntehtic protein data");
    load ./Data/artificial_protein_sequences_and_Delta
    referenced_mutation_effect=Delta_new_matrix;
    input_seq_data=sequence_new;
    starting_index=0;
end
%% remove sites dominated by gaps [threshold is cutoff_gap]

% first convert to Ising representation
[binary_sequence_old,inverse_index_old,conserved_num_old,index_aminoacid_old]=naive_binary_representation(input_seq_data,amino_N);

mean_sequence_0=mean(binary_sequence_old,1);


L0=size(input_seq_data,2);
index=(0:L0-1)*amino_N+1;
if plot_fig
    figureParameter, bar(1:L0,mean_sequence_0(index),'r');
    ylabel('Gap state mean');
    xlabel('Site index');
    fig_name='./figure/beta_mean_gap_state_0.eps';
    figurePostTreat
end


low_gap_index=mean_sequence_0(index)<cutoff_gap;

% update the mutation vector
referenced_mutation_effect=referenced_mutation_effect(:,low_gap_index);
[exp_s1,exp_s2]=size(referenced_mutation_effect);
compressed_mutation_vector=zeros(exp_s2,1);
for j=1:exp_s2  %remove the first 2 and last 2
    temp=referenced_mutation_effect(:,j);
    compressed_mutation_vector(j)=norm(temp);
end

% update input sequences
input_seq_data2=input_seq_data(:,low_gap_index);
L=size(input_seq_data2,2);

gap_p=mean(mean_sequence_0(low_gap_index));

reverse_index=reverse_index_from_logic_input(low_gap_index);

%% using this updated data, and again, convert to Ising representation: a vector with size 21*L
[binary_sequence,inverse_index,conserved_num,index_aminoacid]=naive_binary_representation(input_seq_data2,amino_N);
%[binary_sequence,inverse_index,conserved_num,index_aminoacid]=reduced_representation(input_seq_data2,amino_N);

%% covariance analysis
M_new=size(binary_sequence,1);

if seq_reweighting==1
    Weight=calc_weights(input_seq_data,similarity_threshold); %put more weights on different sequences, suppress the contribution of similar ones. 
    %reweighting at the sequence level, default is reweighting=0 for our paper
    Re_Weight=(1-seq_reweighting)*ones(M_new,1)+seq_reweighting*Weight; 
else
    Re_Weight=ones(M_new,1);
end

site_size=21;
[seq_corr,mean_sequence]=covariance_matrix(binary_sequence,Re_Weight,pseudo_ratio,site_size);


%% Compute the conservation vector phi
if using_PDZ_data
    bg_prob=[0.073, 0.025, 0.050, 0.061, 0.042, 0.072, 0.023, 0.053, 0.064, 0.089, 0.023, 0.043, 0.052, 0.040, 0.052, 0.073, 0.056, 0.063, 0.013, 0.033];
    p=zeros(1,amino_N);
    p(1)=gap_p;
    p(2:amino_N)=(1-gap_p)*bg_prob;
else
    p=zeros(1,amino_N)+1/amino_N; % the background frequency is the same for synthetic data
end

[phi_conserv,compress_phi]=compute_conservation_phi(binary_sequence,pseudo_ratio,p,amino_N,plot_fig);

L=length(conserved_num);
index_k=zeros(L,1);
for k=2:L
    index_k(k)=index_k(k-1)+conserved_num(k-1);
end

%% SCA, or Ranganathan's method

SCA_seq_corr=abs(seq_corr).*abs(phi_conserv*phi_conserv');

if plot_fig
    figureParameter
    f1=image(1:length(phi_conserv),1:length(phi_conserv),1000*SCA_seq_corr);
    colorbar;
    title('SCA C');
    fig_name='./figure/SCA_C.jpg';
    set(gca,'YDir','normal')
    print(fig_name,'-r400','-djpeg');
end

compres_SCA_seq_corr=matrix_compression(index_k,conserved_num,SCA_seq_corr,0);

if plot_fig
    figureParameter
    f1=image(1:L,1:L,1000*abs(compres_SCA_seq_corr));
    colorbar;
    %set(gca, 'clim', [-0.25 0.25]);
    %set(gca, 'clim', [-1 1]);
    %a1=xlabel('$x$');
    %a2=ylabel('$y$');
    title('SCA C');
    fig_name='./figure/SCA_C.jpg';
    set(gca,'YDir','normal')
    print(fig_name,'-r400','-djpeg');
end
%%  eigenvalue analysis of Corr: less accurate scheme
select_data=compres_SCA_seq_corr;%SCA_seq_corr;%
% if use the covariance matrix, change the eigenvector used for prediction.
count=size(select_data,2);
[NormVector_corr,orderEigValue_corr]=orderedEigSystem(select_data,0);

if plot_fig
    clear h1
    figureParameter
    f1=plot(1:count,orderEigValue_corr(1:end,1),'or');
    a1=xlabel('Mode index: $k$');
    a2=ylabel('Eigenvalue: $\lambda_k$');
    xlim([0 count]);
    %set(gca,'YTICK',[10^(-6) 10^(-4) 10^(-2) 1]);
    fig_name='./figure/corr_eig.eps';
    figurePostTreat

    figureParameter
    bar(1:count,[NormVector_corr(:,1)],'r');
    xlabel('Site: l');
    xlim([0 count]);
    fig_name='./figure/corr_mode.eps';
    figurePostTreat
end

%%  comparison between SCA, conservation and experimentally measured mutation effect


 % take the Frobenius norm of a vector
%compress_phi=vector_compression(index_k,conserved_num,phi_conserv);
conserv_prediction=compress_phi;
EXP_mutation_vector=compressed_mutation_vector;

if using_square_root==1
    SCA_prediction=sqrt(abs(NormVector_corr(:,1)));
    disp("Using square root of the eigenvector for prediction");
    
else
    disp("Using original eigenvector for prediction");
    SCA_prediction=abs(NormVector_corr(:,1)); 
end



[SCA_Exp_Corr,top20_index_SCA,top20_index_exp]=comparing_two_predictions(SCA_prediction,"SCA",EXP_mutation_vector,"Experiment",plot_fig);
disp("SCA-Exp correlation is: "+num2str(SCA_Exp_Corr));

[Conserv_Exp_Corr,top20_index_exp,top20_index_conserv]=comparing_two_predictions(EXP_mutation_vector,"Experiment",conserv_prediction,"Conservation",plot_fig);
disp("Exp-Conserv correlation is: "+num2str(Conserv_Exp_Corr));

[Conserv_SCA_Corr,top20_index_SCA,top20_index_conserv]=comparing_two_predictions(SCA_prediction,"SCA",conserv_prediction,"Conservation",plot_fig);
disp("Conserv-SCA correlation is: "+num2str(Conserv_SCA_Corr));



%% Comparison with experimental data

%%%% old data, from mannual estimation.  The ranking is not entirely
%%%% correct. Also, it contains 367, which actually has a smaller mutation
%%%% effect than site 357
ranked_true_hot=reverse_index(top20_index_exp)+starting_index;
%ranked_true_hot=[327   329   359   372   323   324   325   330   341   347   375   379   328   336  338   376   388   353   362   367];

my_vector=abs(NormVector_corr(:,1));
correct_hits_N20=prediction_analysis(my_vector,reverse_index,ranked_true_hot,starting_index,plot_fig);



if save_all_results
%     tt=today;
%     aa=string(datestr(tt,'yy'))+string(datestr(tt,'mm'))+string(datestr(tt,'dd'));  
    output_name="./Data/SCA_data_"+"pseudoCount_"+num2str(pseudo_ratio)+"_cutoffGap_"+num2str(cutoff_gap);
    save(output_name);
end

%% save relevant predictions, combine the prediction from "dealing_protein_data",
%% the actual comparison between different prediction is done in figure_generator, via loading data from "20190228_prediction_wildtype_ref"


% SCA_prediction=abs(NormVector_corr(:,1));
% SCA_prediction_top20sites=reverse_index(index_SCA(1:20))+starting_index;
% conserv_prediction_top20sites=reverse_index(index_conserv(1:20))+starting_index;
% true_sector_top20sites=ranked_true_hot;

% save prediction ICOD_prediction SCA_prediction SCA_prediction_top20sites ICOD_prediction_top20sites conserv_predicton conserv_prediction_top20sites true_sector_top20sites referenced_mutation_vector inverse_index


% load prediction
% figureParameter
% f1=plot(true_sector,1*ones(20),'+r',ICOD_prediction,2*ones(20),'ok',SCA_prediction,3*ones(20),'*g',conserv_prediction,4*ones(20),'sb');
% %xlim([0.5,100.5]);
% ylim([0,10]);
% h1=legend("True sector sites","ICOD","SCA","Conservation");
% set(gca, 'YTick' , []);
% fig_name='./figure/prediction_comparison.eps';
% figurePostTreat



