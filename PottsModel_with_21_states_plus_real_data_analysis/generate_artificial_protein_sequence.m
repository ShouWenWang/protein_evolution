%%
% Delta is generated directly for the reduced representation
% each site can be a value from 1 to 21

close all; clear all;
addpath('./functions');

%% initialization
M=500000; % number of samples
q=2; % the residue to remove at each site
bias_L=0; % Window selection, the left threshold
bias_R=1; %window selection, the right threshold
amino_N=21; % this is the default. Otherwise, there will be error in the down stream analysis
L=10; % length of the sequence
 
%% loading the mutation vector Delta
% lower modes seems to be relevant in our system
Delta_0=randn(1,amino_N*L);

s=1.2;
sector_region=1:5*amino_N;
Delta_0(sector_region)=10*randn(1,5*amino_N);

Delta_0(8)=s*4;
% Delta_0(10)=5;
% Delta_0(20)=2;
Delta_0(amino_N+10)=-s*5;
Delta_0(2*amino_N+6)=s*3;
Delta_0(3*amino_N+20)=-s*4;
Delta_0(4*amino_N+16)=s*2;
figure,plot(Delta_0);

%% generate sequence ensemble
sequence=zeros(M,L);

for k=1:M
    sequence(k,:)=ceil(amino_N*rand(1,L));
end



%% naive binary representation
[binary_sequence_0,inverse_index,conserved_num,index_aminoacid]=naive_binary_representation(sequence,amino_N);
count_0=L*amino_N;


%% selection. Only window selection is provded here
sum_K=zeros(M,1);
for k=1:M
   sum_K(k)=sum(Delta_0.*binary_sequence_0(k,:));
end

plot_distribution(sum_K,[-40 40],40,'$\sum_l\alpha_l\Delta_0_l$');

figureParameter
bar(1:count_0,Delta_0,'b')
xlim([0 count_0]);
ylim([-30 30]);
a1=xlabel('Site index');
a2=ylabel('Old Delta');
fig_name='./figure/old-Delta.eps';
figurePostTreat


index_center=sum_K>mean(sum_K)+bias_L*std(sum_K) & sum_K<mean(sum_K)+bias_R*std(sum_K);

% the selected sequence
sequence_new=sequence(index_center,:);

%% use the reference q to generate the corresponding referenced mutation effect
Delta_new_matrix=zeros(amino_N-1,L); % a matrix
Delta_new_vector=zeros(1,(amino_N-1)*L);
for j=1:L
    old_index_1=(j-1)*amino_N+1:(j-1)*amino_N+q-1;
    old_index_2=(j-1)*amino_N+q+1:j*amino_N;
    old_index=cat(2,old_index_1,old_index_2);
    
    
    Delta_new_matrix(:,j)=Delta_0(old_index)-Delta_0((j-1)*amino_N+q); 
    new_index=(j-1)*(amino_N-1)+1:j*(amino_N-1);
    Delta_new_vector(new_index)=Delta_0(old_index)-Delta_0((j-1)*amino_N+q);
end
    
figureParameter
bar(1:length(Delta_new_vector),Delta_new_vector,'b')
xlim([0 count_0]);
ylim([-30 30]);
a1=xlabel('Site index');
a2=ylabel('Referenced Delta');
fig_name='./figure/referenced-Delta.eps';
figurePostTreat



save ./Data/artificial_protein_sequences_and_Delta  sequence_new Delta_new_matrix q Delta_new_vector amino_N
