function Delta_new=Delta_conversion(Delta_0,conserved_num_old,index_aminoacid_old,amino_N)
%%
% after using the reduced representation,  we need also suppress certain
% sites in Delta 

%%
L=length(conserved_num_old);
index_k=zeros(L,1); % the starting index at kth amino acid in the new representation
for k=2:L
    index_k(k)=index_k(k-1)+conserved_num_old(k-1);
end

Delta_new=zeros(1,sum(conserved_num_old));

for k=1:L

 index_temp_new=index_k(k)+1:index_k(k)+conserved_num_old(k);
 index_temp_old=(k-1)*amino_N+index_aminoacid_old(index_temp_new);
 Delta_new(index_temp_new)=Delta_0(index_temp_old);
end