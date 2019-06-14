function single_point_mutation


%% important parameters
%PDB_name_1='1be9.pdb'; reference_index_1=301;
PDB_name_1='1BFE.pdb'; reference_index_1=306;
cut_index_low=309;cut_index_high=398;
cutoff_distance=7.5;
factor=0.8; % perturbation factor, close to 1 is a small perturbation

%eig_index=3;
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
[NormVector,orderEigValue,coord_normVector,spring]=normal_mode_computation_Cbeta_mutation(coords,Atom_type,cutoff_distance,factor);

%% generate a random sequence

Atom_type_mutant=Atom_type;
%Atom_type_mutant(2)=-1;

M=count;

eigen_array_many=zeros(6,M);
vector_array_many=zeros(6,M);

sequence_array_many=zeros(count,M);
for j=1:M
sequence=zeros(count,1);
sequence(j)=1;
Atom_type_mutant(beta_atom_index)=sequence; 
[NormVector_mutant,orderEigValue_mutant,coord_normVector_mutant,spring,Hessen_2d]=normal_mode_computation_Cbeta_mutation(coords,Atom_type_mutant,cutoff_distance,factor);
eigen_array_many(:,j)=orderEigValue_mutant(7:12,1)-orderEigValue(7:12,1);
vector_array_many(:,j)=sum(NormVector_mutant(:,7:12).*NormVector(:,7:12));
sequence_array(:,j)=sequence;
end



%%
Delta_lambda2=eigen_array_many(2,:);
Delta_lambda3=eigen_array_many(3,:);
Delta_lambda4=eigen_array_many(4,:);
Delta_lambda6=eigen_array_many(6,:);
%Delta_lambda=eigen_array-orderEigValue(6+eig_index,1);

plot_distribution(Delta_lambda,[-0.003 0],40,'$\Delta_6^l$');

figureParameter
bar(1:count,-Delta_lambda3','r');
fig_name='./figure/corr_mode.eps';
a1=ylabel('$-\Delta_4^l$');
%set(h1,'location','south')
%ylim([ 0 0.0025]);
%set(gca,'YTICKLABEL',[-0.1 0 0.1 0.2 0.3 0.4 0.5]);
figurePostTreat


%%
temp_c=0;
index_Delta=0;
value_Delta=0;
for i=1:length(Delta_lambda)
    if abs(Delta_lambda(i))>0.0005
    temp_c=temp_c+1;
    index_Delta(temp_c)=i; 
    value_Delta(temp_c)=Delta_lambda(i);
    
    end
end
index_Delta=index_Delta';
value_Delta=value_Delta';



Amino_acid_loc=visualizing_hotspot(coords,beta_atom_index,index_Delta,spring,reference_index_1,Atom_type);

