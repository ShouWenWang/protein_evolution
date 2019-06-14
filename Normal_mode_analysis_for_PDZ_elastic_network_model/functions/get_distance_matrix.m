function Dist_ij=get_distance_matrix

PDB_name_1='1be9.pdb'; reference_index_1=301;
PDB_name_2='1BFE.pdb'; reference_index_2=306;
%cut_index_low=309;cut_index_high=398;
cut_index_low=312;cut_index_high=390;

coords_1=read_PDB_data_C_alpha(PDB_name_1,reference_index_1,cut_index_low,cut_index_high);

%% Get the distance matrix
data=coords_1;
[L1,L2]=size(data);
Dist_ij=zeros(L1,L1)+0.0001;

for j=1:L1
    for k=1:L1
        temp=data(j,:)-data(k,:);
        Dist_ij(j,k)=sqrt(sum(temp.^2));
        
    end
end

inv_Dist_ij=1./Dist_ij;
for j=1:L1
    inv_Dist_ij(j,j)=0;     
end

% 
% figureParameter
% f1=image(1:L1,1:L1,200*inv_Dist_ij);
% colorbar;
% fig_name='./figure/corr1.jpg';
% set(gca,'YDir','normal')
% print(fig_name,'-r500','-djpeg');