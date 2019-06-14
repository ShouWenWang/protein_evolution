function search_distant_hotspot_inverse_C(inv_C,true_hot_index,starting_index)
%% find hotspot
% in_C, the inverse matrix
% starting_index=312 for our data 
%%
comp_index=true_hot_index-starting_index;
L1=length(comp_index);


[hot_index_pairs,hotspot_data]=search_hotspot_pairs(inv_C,0.05);
Dist_ij=get_distance_matrix;
L0=length(hotspot_data);
distance_data=zeros(L0,1);
temp_1=0;
temp_2=0;
count_1=0;
count_2=0;
contact_ratio=zeros(L0,1);
functional_ratio=zeros(L0,1);

cutoff_distance=11;

for j=1:L0
    distance_data(j)=Dist_ij(hot_index_pairs(j,1),hot_index_pairs(j,2));
    if abs(hot_index_pairs(j,1)-hot_index_pairs(j,2))>2
        count_1=count_1+1;
        if distance_data(j)<cutoff_distance 
        temp_1=temp_1+1;
        else
          % matching the hotspot: at least one of the pair is the hotspot 
           count_2=count_2+1;
           if find(hot_index_pairs(j,1)==comp_index) & find(hot_index_pairs(j,2)==comp_index)
               temp_2=temp_2+1;  
           end
           functional_ratio(count_2)=temp_2/count_2;
        end
        contact_ratio(count_1)=temp_1/count_1;
         
    end
    
end



figureParameter
f1=plot(1:L0,distance_data(1:L0),'*r');
a1=xlabel('count');
a2=ylabel('pair distance');
xlim([1,150]);
fig_name='./figure/DCA_pair_distance.eps';
figurePostTreat

figureParameter
f1=plot(1:L0,contact_ratio(1:L0),'*r');
a1=xlabel('count');
a2=ylabel('Contact ratio');
xlim([1,150]);
ylim([0,1]);
fig_name='./figure/DCA_pair_distance.eps';
figurePostTreat

figureParameter
f1=plot(1:L0,functional_ratio(1:L0),'*r');
a1=xlabel('count');
a2=ylabel('Functional ratio');
xlim([1,100]);
ylim([0,1]);
fig_name='./figure/DCA_pair_distance.eps';
figurePostTreat