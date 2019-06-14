function Norm_Involvement=conformational_change_Calpha
% conformational chanage
PDB_name_1='1be9.pdb'; reference_index_1=301;
PDB_name_2='1BFE.pdb'; reference_index_2=306;
%cut_index_low=309;cut_index_high=398;
cut_index_low=312;cut_index_high=390;

coords_1=read_PDB_data_C_alpha(PDB_name_1,reference_index_1,cut_index_low,cut_index_high);
coords_2=read_PDB_data_C_alpha(PDB_name_2,reference_index_2,cut_index_low,cut_index_high);


figureParameter
f1=plot3(coords_1(:,1),coords_1(:,2),coords_1(:,3),'-*b');
hold on;
plot3(coords_2(:,1),coords_2(:,2),coords_2(:,3),'-*r')
hold off;
fig_name='./conform_correlation.eps';
figurePostTreat

target_atom_index=20;% for rotation
rotation_on=0; % disable rotaton
new_coords_1=translation_rotation(coords_1,target_atom_index,rotation_on);
new_coords_2=translation_rotation(coords_2,target_atom_index,rotation_on);




figureParameter
f1=plot3(new_coords_1(:,1),new_coords_1(:,2),new_coords_1(:,3),'-*b');
hold on;
plot3(new_coords_2(:,1),new_coords_2(:,2),new_coords_2(:,3),'-*r')
hold off;
fig_name='./conform_correlation.eps';
figurePostTreat


%% remove the rotational effect and translational effect



%%

coords_1=new_coords_1;
coords_2=new_coords_2;

displacement=coords_1-coords_2;
info=size(coords_1);
N=info(1);
distance=0;
for i=1:N
   distance=distance+sum(displacement(i,:).^2);       
end

displacement_norm=displacement./sqrt(distance);

   
correlation=zeros(N,N);


for i=1:N
    for j=1:N
       correlation(i,j)=sum(displacement_norm(i,:).*displacement_norm(j,:));    
    end
end


meshsize=4;
coarse_graining_corr(correlation,meshsize);



%%  involvement coefficient


projection=zeros(3*N-6,1);
cutoff_distance=7.5;
[NormVector,orderEigValue,coord_normVector,spring]=normal_mode_computation_Calpha(coords_1,cutoff_distance);

% relies on the case that the zero eigenvalue are in the region: 3*N-5:3*N
for i=1:3*N-6
    projection(i)=0;
    for j=1:N
        for s=1:3
        projection(i)=projection(i)+displacement_norm(j,s).*coord_normVector(j,6+i,s);
        end
    end
end

Norm_Involvement=zeros(3*N-6,1);
Norm_Involvement=abs(projection./orderEigValue(7:3*N,1));
Norm_Involvement=Norm_Involvement./sqrt(sum(Norm_Involvement.^2));

figureParameter
f1=plot(1:length(Norm_Involvement),Norm_Involvement,'-*r');
xlim([0 20]);
a1=xlabel('$n$');
a2=ylabel('$T_n$');
fig_name='./thermal_involvement.eps';
figurePostTreat


figureParameter
f1=plot(1:length(projection),abs(projection),'-*r');
xlim([0 20]);
a1=xlabel('$n$');
a2=ylabel('$|I_n|$');
fig_name='./thermal_projection.eps';
figurePostTreat



