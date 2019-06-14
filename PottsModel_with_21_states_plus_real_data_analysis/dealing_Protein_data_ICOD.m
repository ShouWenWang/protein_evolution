function [ICOD_Exp_Corr,ICOD_Conserv_Corr,Conserv_Exp_Corr,correct_hits_N20]=dealing_Protein_data_ICOD(pseudo_ratio,cutoff_gap,using_PDZ_data,plot_fig)
%%%%%%%%%%%%%%%%%%%
% last updated, 2019-2-26
% This program is used to generate Fig.S18 in our paper
% Method of ICOD is employed 


% clear all; 
% close all;

%% Initialization 
%cutoff_gap=0.3; % throw away the gap state (default, 0.3)
%using_PDZ_data=1; % 1, use the PDZ data; 0, use the syntehtic data (IN this case, q is predefined and will be changed accordingly when loading data)


%%% -1 for removing the least conserved sites, the default choice for the paper; 
   % positive number, removing the corresponding residue; 
   % -2, this will be convert to the wildtype index later, which becomes a vector
   % a vector, remove corresponding residues as each site
   % 
q=-2;

%plot_fig=1; %do not plot figures
%pseudo_ratio=0.010; %  default: 0.005 for the PDZ data in the paper; maximum is 1;
zero_sum_gauge=0;  % 1: apply zero sum gauge,  0: do not apply, the new default after 20190226
compress=1; %0, do not apply Frobenius norm; 1, apply Frobenius compression. 1 is the default setting for our paper

save_all_results=1;
predict_inverse_entries=0;
clean_contact=0; %0, do not clean contacts [the default for our paper]
seq_reweighting=0;   % 1, apply reweighting to each sequence; 0, do not apply. In this real data, we do not apply reweighting
similarity_threshold=0.9;

%% 
addpath('./functions/');

if ~exist('figure','dir')
    mkdir figure;
end

if ~exist('Data','dir')
    mkdir hh
end

%% print out a summary of the condition used!

if using_PDZ_data
    disp("Running with PDZ experimental data");
    %% loading the experimentally measured mutation effect from McLaughlin_Ranganathan_data, and construct the wildtype references
    % wildtype residues for PDZ, starting from position 313 to 391 
    % wildtype=char("RIVIHRGSTGLGFNIVGGEDGEGIFISFILAGGPADLSGELRKGDQILSVNGVDLRNASHEQAAIALKNAGQTVTIIAQ");

    % in our representation, at each site, the residue ordering is [the same as
    % the McLaughlin12 Nature paper, map is given in "number2letter_map" 

    load ./Data/McLaughlin_Ranganathan_data.mat;
    
    starting_index=312;  % in reference to the real protein index
    amino_N=21;
    % the mapping to Aminoacid is provided in "number2letter_map" from this
    % data

    wildtype_index=zeros(79,1);
    for j=1:79
        wildtype_index(j)=find(number2letter_map==wildtype(j),1);
    end

    if q==-2
        q=wildtype_index;
        disp("Using wildtype reference");
    else
        disp("Reference state: "+num2str(q));
    end

    referenced_mutation_effect=single_mutagenesis_cript.data;
    referenced_mutation_effect=referenced_mutation_effect(:,3:81); %20*79 is the dimension
    
    %% Loading the available sequence ensemble for PDZ from PDZal.
    load ./Data/PDZal;
    %rand_sample=cumsum(floor(10*rand(2000,1))+1);
    %input_seq_data=efa1(rand_sample,:);
    input_seq_data=efa1;
    
else
    disp("Running with syntehtic protein data");
    load ./Data/artificial_protein_sequences_and_Delta
    referenced_mutation_effect=Delta_new_matrix;
    input_seq_data=sequence_new;
    starting_index=0;
    disp("Reference state: "+num2str(q)); % q is also defined here
    
end

%% remove sites dominated by gaps [threshold is cutoff_gap]

% first convert to Ising representation
[binary_sequence_old,inverse_index_old,conserved_num_old,index_aminoacid_old]=naive_binary_representation(input_seq_data,amino_N);

% then compute the average presence for each residue
mean_sequence_0=mean(binary_sequence_old,1);


L0=size(input_seq_data,2);
index=(0:L0-1)*amino_N+1;
    
if plot_fig
    figureParameter, bar(1:L0,mean_sequence_0(index),'r');
    ylabel('Gap state mean');
    xlabel('Site: l');
    fig_name='./figure/beta_mean_gap_state_0.eps';
    figurePostTreat
end

% L0=size(input_seq_data,2);
% index=(0:L0-1)*amino_N+1;
% figureParameter, bar(1:(L0*amino_N),mean_sequence_0,'r');
% ylabel('Residue-level state');
% xlabel('site: l');
% %fig_name='./figure/beta_mean_residue.eps';
% figurePostTreat



low_gap_index=mean_sequence_0(index)<cutoff_gap;

referenced_mutation_effect=referenced_mutation_effect(:,low_gap_index);
[exp_s1,exp_s2]=size(referenced_mutation_effect);
referenced_mutation_vector=zeros(exp_s1*exp_s2,1);
compressed_mutation_vector=zeros(exp_s2,1);
for j=1:exp_s2  %remove the first 2 and last 2
    temp=referenced_mutation_effect(:,j);
    referenced_mutation_vector((j-1)*(amino_N-1)+1:j*(amino_N-1))=temp;
    compressed_mutation_vector(j)=norm(temp);
end


if length(q)>1 % q is a vector
    q=q(low_gap_index);
end

% update data
input_seq_data2=input_seq_data(:,low_gap_index);
L=size(input_seq_data2,2); % the new sequence length
gap_p=mean(mean_sequence_0(low_gap_index));

% reverse_index take the index from the new representation and convert it
% to the original representation
reverse_index=reverse_index_from_logic_input(low_gap_index);

%% using this updated data, and again, convert to Ising representation: a vector with size 21*L
[binary_sequence_new,inverse_index_new,conserved_num_new,index_aminoacid_new]=naive_binary_representation(input_seq_data2,amino_N);
index_aminoacid_new=index_aminoacid_new';
mean_sequence_old_2=mean(binary_sequence_new,1);
%[binary_sequence_new,inverse_index_new,conserved_num_new,index_aminoacid_new]=reduced_representation(input_seq_data2,amino_N);

%% compute conservation
% the background probability for each residue, used for computing the
% conservation vector phi_l, needed if we want to get rid of degeneracy
% according to conservation
if using_PDZ_data
    bg_prob=[0.073, 0.025, 0.050, 0.061, 0.042, 0.072, 0.023, 0.053, 0.064, 0.089, 0.023, 0.043, 0.052, 0.040, 0.052, 0.073, 0.056, 0.063, 0.013, 0.033];
    p=zeros(1,amino_N);
    p(1)=gap_p;
    p(2:amino_N)=(1-gap_p)*bg_prob;
else
    p=zeros(1,amino_N)+1/amino_N; % the background frequency is the same for synthetic data
end

[phi_conserv_0,compress_phi]=compute_conservation_phi(binary_sequence_new,pseudo_ratio,p,amino_N,plot_fig);
%% remove degeneracy by getting rid of one residue at each site, determined by q

% q determines which residue to remove, -1 for most abundant; 
[conserved_num,index_k,index,remove_q]=remove_residue_q(binary_sequence_new,index_aminoacid_new,conserved_num_new,q,pseudo_ratio);
% actually get rid of the least conserved site, and update corresponding
% variables
binary_sequence=binary_sequence_new(:,index);
index_aminoacid=index_aminoacid_new(index);
inverse_index=inverse_index_new(index);
phi_conserv=phi_conserv_0(index); 



%% covariance analysis
M_new=size(binary_sequence,1);

if seq_reweighting==1
    Weight=calc_weights(input_seq_data,similarity_threshold); %put more weights on different sequences, suppress the contribution of similar ones. 
    %reweighting at the sequence level, default is reweighting=0 for our paper
    Re_Weight=(1-seq_reweighting)*ones(M_new,1)+seq_reweighting*Weight; 
else
    Re_Weight=ones(M_new,1);
end
    

site_size=amino_N-1;
[seq_corr,mean_sequence]=covariance_matrix(binary_sequence,Re_Weight,pseudo_ratio,site_size);

figure,bar(1:length(mean_sequence),mean_sequence);

if plot_fig
    count=size(seq_corr,2);
    figureParameter
    f1=image(1:count,1:count,10000*seq_corr);
    colorbar;
    title('seq corr')
    fig_name='./figure/seq_corr.jpg';
    set(gca,'YDir','normal')
    print(fig_name,'-r400','-djpeg');
end

%% reweighting at the covariance matrix level, disabled here
%compres_seq_corr=matrix_compression(index_k,conserved_num,seq_corr,0);

% count_1=size(seq_corr,1);
% reweight_seq_corr=zeros(count_1,count_1);
% for j=1:count_1
%     for k=j:count_1
%        if j==k
%        reweight_seq_corr(j,j)=0; % this works equally well compared with more elaborated choice  
%        else
%        reweight_seq_corr(j,k)=seq_corr(j,k)/(seq_corr(j,j)*seq_corr(k,k));
%        %reweight_seq_corr(j,k)=seq_corr(j,k);
%        reweight_seq_corr(k,j)=reweight_seq_corr(j,k);
%        end
%     end
% end


%% inverse matrix 
inv_corr=inv(seq_corr);
inv_corr_off_diag=inv_corr;
for k=1:L
    temp_index=index_k(k)+1:index_k(k)+conserved_num(k);
    inv_corr_off_diag(temp_index,temp_index)=0;
end

count=size(seq_corr,2);

if plot_fig
    figureParameter
    f1=image(1:count,1:count,2*abs(inv_corr_off_diag));
    colorbar;
    %set(gca, 'clim', [-0.25 0.25]);
    %set(gca, 'clim', [-1 1]);
    title('inv corr zero off diag')
    fig_name='./figure/inv_corr_off_diag.jpg';
    set(gca,'YDir','normal')
    print(fig_name,'-r400','-djpeg');
end

%% predicting inverse matrix entries
if predict_inverse_entries
    predicting_inverse_matrix_entries(mean_sequence,mean_sequence_old_2,remove_q,inv_corr,amino_N,referenced_mutation_vector);
end

%% compression: taking Frobenius norm
if compress
    compres_inv_C=matrix_compression(index_k,conserved_num,inv_corr_off_diag,zero_sum_gauge);
    
    if zero_sum_gauge    
        disp("Compress the inverse correlation matrix with zero sum gauge"); 
    else
        disp("Compress the inverse correlation matrix without zero sum gauge")
    end
    
    if plot_fig
        figureParameter
        f1=image(1:L,1:L,0.5*abs(compres_inv_C));
        colorbar;
        fig_name='./figure/compres_inv_C.jpg';
        title('compress inv C');
        set(gca,'YDir','normal')
        print(fig_name,'-r400','-djpeg');

    end
    
    new_inv_C=compres_inv_C;
   
else
    disp("Do not compress the inverse correlation matrix"); 
    new_inv_C=inv_corr_off_diag;
    
end
%% clean close contacts

if clean_contact==1
    disp("Clean contacts");

    C_vector=new_inv_C(:);

    figure,hist(C_vector,100)
    
    for j=1:L
      for k=1:L
         if abs(new_inv_C(j,k))>500 %&& abs(j-k)<3
             new_inv_C(j,k)=0;
         end
      end
    end


    if plot_fig
        figureParameter
        f1=image(1:L,1:L,0.5*abs(new_inv_C));
        colorbar;
        title('inv C: remove contact');
        fig_name='./figure/new_inv_C.jpg';
        set(gca,'YDir','normal')
        print(fig_name,'-r400','-djpeg');
    end

end

%%  eigenvalue analysis of Corr: less accurate scheme [SVD is more accurate]. But for our problem, it does not matter
if compress==1
  select_data=new_inv_C;  % contacts are not cleaned for default option. But can be changed 
  
  % take the Frobenius norm of a vector
  %conserv_prediction=vector_compression(index_k,conserved_num,phi_conserv);
  conserv_prediction=compress_phi;
  EXP_mutation_vector=compressed_mutation_vector; 
  %EXP_mutation_vector=EXP_mutation_vector';

else
  
  select_data=inv_corr_off_diag;  % contacts are not cleaned
  conserv_prediction=phi_conserv;
  EXP_mutation_vector=referenced_mutation_vector;
end

count=size(select_data,2);
[NormVector_corr,orderEigValue_corr]=orderedEigSystem(select_data,0);



clear h1
if plot_fig
    figureParameter
    f1=plot(1:count,orderEigValue_corr(1:end,1),'or');
    xlabel('Mode index');
    ylabel('Eigenvalue');
    xlim([0 count]);
    %set(gca,'YTICK',[10^(-6) 10^(-4) 10^(-2) 1]);
    fig_name='./figure/corr_eigenvalue.eps';
    figurePostTreat


    figureParameter
    bar(1:count,NormVector_corr(:,1),'r');
    %plot(1:count,50*[NormVector_corr(:,1)],'.r',1:length(phi_conserv),phi_conserv,'.b');%,1:count,inverse_index(1:count),'.g');
    xlabel('Site index');
    xlim([0 count]);
    fig_name='./figure/corr_eigenmode.eps';
    figurePostTreat
end

%% Comparion with conservation, and experimental data

ICOD_prediction=abs(NormVector_corr(:,1));
%conserv_prediction=abs(phi_conserv);

[ICOD_Exp_Corr,top20_index_ICOD,top20_index_exp]=comparing_two_predictions(ICOD_prediction,"ICOD",EXP_mutation_vector,"Experiment",plot_fig);
disp("Exp-ICOD correlation is: "+num2str(ICOD_Exp_Corr));

[Conserv_Exp_Corr,top20_index_exp,top20_index_conserv]=comparing_two_predictions(EXP_mutation_vector,"Experiment",conserv_prediction,"Conservation",plot_fig);
disp("Exp-Conserv correlation is: "+num2str(Conserv_Exp_Corr));

[ICOD_Conserv_Corr,top20_index_ICOD,top20_index_conserv]=comparing_two_predictions(ICOD_prediction,"ICOD",conserv_prediction,"Conservation",plot_fig);
disp("ICOD-Conservation correlation: "+num2str(ICOD_Conserv_Corr));

%% comparison with experimental data
correct_hits_N20=[];
if compress==1
    
%%%% old data, from mannual estimation.  The ranking is not entirely
%%%% correct. Also, it contains 367, which actually has a smaller mutation
%%%% effect than site 357
% ranked_true_hot=[327   329   359   372   323   324   325   330   341   347   375   379   328   336  338   376   388   353   362   367];

%%% new data, based on the top 20 sites from experiment: referenced_mutation_vector,
%%% ranked in the descending order, this is the default option after
%%% 2019-02-28

      ranked_true_hot=reverse_index(top20_index_exp)+starting_index;
    % This is the resulting hotspot for the PDZ model
    % ranked_true_hot=[338   327   359   325   323   329   379   324   388   375   328   341   347   353   372   330   336   362 357   376];


     correct_hits_N20=prediction_analysis(NormVector_corr(:,1),reverse_index,ranked_true_hot,starting_index,plot_fig);

%else
%    [ICOD_Exp_Corr,top20_index_ICOD,top20_index_exp]=comparing_two_predictions(ICOD_prediction,"ICOD-noCompre",referenced_mutation_vector,"Exp-noCompre");

%     % compress first
%     compres_ICOD=vector_compression(index_k,conserved_num,NormVector_corr(:,1));
%     compres_ICOD=compres_ICOD';
%     compres_exp=vector_compression(index_k,conserved_num,referenced_mutation_vector);
%     compres_exp=compres_exp';
%     corr_x1x2=comparing_two_predictions(compres_ICOD,"ICOD-Comp",compres_exp,"Experiment-Comp");
% 
%     prediction_analysis(compres_ICOD,reverse_index,ranked_true_hot,starting_index);

end


if save_all_results
%     tt=today;
%     aa=string(datestr(tt,'yy'))+string(datestr(tt,'mm'))+string(datestr(tt,'dd'));
%     output_name="./Data/"+aa+"_ICOD_data_"+"pseudoCount_"+num2str(pseudo_ratio)+"_cutoffGap_"+num2str(cutoff_gap);
    output_name="./Data/ICOD_data_"+"pseudoCount_"+num2str(pseudo_ratio)+"_cutoffGap_"+num2str(cutoff_gap);
    save(output_name);
end


%% %% save relevant predictions, combine the prediction from "dealing_protein_data_SCA", the actual comparison between different prediction is done in figure_generator


% ICOD_prediction_top20sites=reverse_index(top20_index_ICOD)+starting_index;
% conserv_prediction_top20sites=reverse_index(top20_index_conserv)+starting_index;
% true_sector_top20sites=reverse_index(top20_index_exp)+starting_index;
% reverse_index_ICOD=reverse_index;
% 
% save prediction ICOD_prediction ICOD_prediction_top20sites conserv_prediction conserv_prediction_top20sites true_sector_top20sites referenced_mutation_vector reverse_index_ICOD

%% check how the wildtype reference differs from the current reference
% figureParameter
% plot(remove_q,wildtype_index,'*b');
% xlabel("Most abundant residue");
% ylabel("Wilde type")
% % xlim([0.00001,0.5]);
% set(gca,'XTICK',[0 5 10 15 20]);
% % ylim([0.0001,5])
% %legend("Eigenvector","conservation");
% %xlim([0 count]);
% fig_name='./figure/abundant_residue_wild_type.eps';
% figurePostTreat

