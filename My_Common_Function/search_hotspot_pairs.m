function [hot_index_pairs,hotspot_data]=search_hotspot_pairs(data_org,show_ratio)

%%
  % select the data that lies outside the range of the distance from the
  % mean

  % for 1-d vector, hot_index_x or hot_index_y could be 1, which does not
  % mean that index 1 is involved.

%%

data=abs(data_org);


[L1,L2]=size(data);

if L1~=L2
    error('2-d square matrix required');
end


pairs=zeros(L1,2);
value=zeros(L1,1);
count=0;
for i=1:L1
    for j=i+3:L1 %skip the next two neighbors
      count=count+1;
      pairs(count,1)=i; 
      pairs(count,2)=j; 
      value(count)=data(i,j);
    end
end
    




   [value_new,index]=sort(value,'descend');
    pairs_new=pairs(index,:);
    
     %show_ratio=0.1;
     upper_index=floor(show_ratio*count);
     hot_index_pairs=pairs_new(1:upper_index,:);
     hotspot_data=value_new(1:upper_index);
   
     
    figureParameter
    f1=semilogx(1:upper_index,hotspot_data,'*');
    a1=xlabel('count');
    a2=ylabel('$|e_{ij}|$');
    xlim([0 100]);
    fig_name='./figure/hotspot_pairs.eps';
    figurePostTreat




        
        
    
