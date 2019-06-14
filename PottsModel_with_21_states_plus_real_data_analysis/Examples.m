%% if using PDZ data provided in the Data folder
clear all;
close all;
using_PDZ_data=1; % use syntehtic data
plot_fig=1; %plot figure
pseudo_ratio=0.02;
cutoff_gap=0.15;

[ICOD_Exp_Corr,ICOD_Conserv_Corr,Conserv_Exp_Corr,correct_hits_N20]=dealing_Protein_data_ICOD(pseudo_ratio,cutoff_gap,using_PDZ_data,plot_fig);
[SCA_Exp_Corr,Conserv_Exp_Corr,Conserv_SCA_Corr,correct_hits_N20]=dealing_Protein_data_SCA(pseudo_ratio,cutoff_gap,using_PDZ_data,plot_fig);

%% Other options for ICOD
% if you want to use your own protein family data, you need to go into the file 
% "dealing_Protein_data_ICOD" and "dealing_Protein_data_SCA", and modify the data directory and name  
% there are many other options for the program, 

% For the ICOD program, you can go into the file, and in the initialization
% section, you can modify: 

%%%%%%% the reference sequence "q"
% -1 for removing the least conserved sites, the default choice for the paper; 
% positive number, removing the corresponding residue; 
% -2, this will be convert to the wildtype index later, which becomes a vector
% a vector, remove corresponding residues as each site

%%%%%%%  plot_fig: 1 plot figure; 0, do not plot figures, 

%%%%%%% zero_sum_gauge: 1,apply zero sum gauge;  0, do not apply, the 
% option for real data analysis is not to apply zero sum gauge

%%%%%%% compress:  %0, do not apply Frobenius norm; 1, apply Frobenius compression. 1 is the default setting for our paper

%%%%%%% save_all_results: 1, save results in the ./Data for later analysis;
% 0, not to save data

%%%%%%% predict_inverse_entries: 0, predict and plot; 1, not to. This is
% used to check our theory about the inverse entries.

%%%%%%% clean_contact: 0, do not clean contacts [the default for our
% paper]; 1, clear contacts. The threshold need to be tuned at the step of
% execution

%%%%%%% seq_reweighting: % 1, apply reweighting to each sequence (say, similar sequences have a less weight than fresh new sequences); 
% 0, do not apply. In this real data, we do not apply reweighting

%%%%%%% similarity_threshold: This is associated to the sequence reweighting. If reweighting is
%true, sequence with a similarlity larger than 90% is less weighted.

%% Similarly, you can decide how to run SCA in the initialization section
%%%%%% using_square_root: 1 use the square root of the SCA vector for prediction, consistent with the rest of the paper
% 0, use the original vector

%%%%%% plot_fig: 0 or 1

%%%%%% seq_reweighting: 0 or 1

%%%%%% similarity_threshold: [0,1]. The default is 0.9. Associated with
% seq_reweighting

%%%%%% save_all_results: 0 or 1;

%% if using artificial protein sequences
clear all;
close all;

generate_artificial_protein_sequence; % the data will be saved to /Data folder
using_PDZ_data=0; % use syntehtic data
plot_fig=1; %plot figure
pseudo_ratio=0.02;
cutoff_gap=0.15;

[ICOD_Exp_Corr,ICOD_Conserv_Corr,Conserv_Exp_Corr,correct_hits_N20]=dealing_Protein_data_ICOD(pseudo_ratio,cutoff_gap,using_PDZ_data,plot_fig);
[SCA_Exp_Corr,Conserv_Exp_Corr,Conserv_SCA_Corr,correct_hits_N20]=dealing_Protein_data_SCA(pseudo_ratio,cutoff_gap,using_PDZ_data,plot_fig);



