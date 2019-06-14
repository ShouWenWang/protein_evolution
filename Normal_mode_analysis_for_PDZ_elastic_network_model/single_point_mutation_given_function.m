function Delta_lambda=single_point_mutation_given_function


%% important parameters
%PDB_name_1='1be9.pdb'; reference_index_1=301;
PDB_name_1='1BFE.pdb'; reference_index_1=306;
cut_index_low=309;cut_index_high=398;
cutoff_distance=7.5;
factor=0.8; % perturbation factor, close to 1 is a small perturbation

[Norm_Involvement,displace_vector_rescale,displace_vector_orig]=conformational_change_Cbeta;
displace_vector=displace_vector_orig;
%% start

[Atom_type,coords]=read_PDB_data_C_beta(PDB_name_1,reference_index_1,cut_index_low,cut_index_high);

% get the positions of the C beta atoms
N=length(Atom_type);
beta_atom_index=zeros(2,1); % dynamically expand
count=0;

reverse_index_beta=zeros(1); % reverse index of beta along the chain
count_alpha=0;
for i=1:N
   if Atom_type(i)==-3; % C alpha
       count_alpha=count_alpha+1;
   end
   
   if Atom_type(i)==0
       count=count+1;
       beta_atom_index(count)=i;
       reverse_index_beta(count)=count_alpha;
   end
end

reverse_index_beta=reverse_index_beta+308;
%% original eigenvalues
[NormVector,orderEigValue,coord_normVector,spring,Hessen_2d]=normal_mode_computation_Cbeta_mutation(coords,Atom_type,cutoff_distance,factor);
eigen_ref=displace_vector'*Hessen_2d*displace_vector;
%% generate a random sequence

Atom_type_mutant=Atom_type;
%Atom_type_mutant(2)=-1;

M=count;
eigen_array=zeros(M,1);
sequence_array=zeros(count,M);
for j=1:M
sequence=zeros(count,1);
sequence(j)=1;
Atom_type_mutant(beta_atom_index)=sequence; 
[NormVector_mutant,orderEigValue_mutant,coord_normVector_mutant,spring,Hessen_2d_mutant]=normal_mode_computation_Cbeta_mutation(coords,Atom_type_mutant,cutoff_distance,factor);
eigen_array(j)=displace_vector'*Hessen_2d_mutant*displace_vector;
sequence_array(:,j)=sequence;
end

Delta_lambda=eigen_array-eigen_ref;


figureParameter
bar(1:count,Delta_lambda,'r');
a1=xlabel('Residue index: $l$');
fig_name='./figure/Delta.eps';
a2=ylabel('$\Delta_l$');
%set(h1,'location','south')
%ylim([ 0 0.0025]);
%set(gca,'YTICKLABEL',[-0.1 0 0.1 0.2 0.3 0.4 0.5]);
figurePostTreat


% 
%     temp_c=0;
%     index_Delta=0;
%     value_Delta=0;
%     for i=1:length(Delta_lambda)
%         if abs(Delta_lambda(i))>0.0025
%         temp_c=temp_c+1;
%         index_Delta(temp_c)=i; 
%         value_Delta(temp_c)=Delta_lambda(i);
% 
%         end
%     end
%     index_Delta=index_Delta';
%     value_Delta=value_Delta';
% 
%     plot_distribution(Delta_lambda,[-0.04 0],20);
% 
%     Amino_acid_loc=visualizing_hotspot(coords,beta_atom_index,index_Delta,spring,cut_index_low,Atom_type);

