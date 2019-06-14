function [hot_index_x,hotspot_data]=search_hotspot(data,expect_value,range)

%%
  % select the data that lies outside the range of the distance from the
  % mean

  % for 1-d vector, hot_index_x or hot_index_y could be 1, which does not
  % mean that index 1 is involved.

%%

  show_hotspot=1;
  rank_yes=0;
  

info=size(data);

if info(1)==1
    data=data'; % vertical vector
    info=size(data); % now info(2)=1;
end


hot_index_x=0;
x_count=0;

for i=1:info(1)
    for j=1:info(2)
     if abs(data(i,j)-expect_value)>range(1) &&  abs(data(i,j)-expect_value)<range(2)    
         if  x_count==0
             x_count=x_count+1;
             hot_index_x(x_count)=i;

             
         else 
             flag_x=sum(hot_index_x==i);
             if flag_x==0
                 x_count=x_count+1;
                 hot_index_x(x_count)=i;
             end
             
             
         end
             
     end
    end
        
end

hot_index_x=hot_index_x';



if x_count==0
   error('No hotspot found, please change range')
end





%% rank the hotspot index by its relative importance

%% 2-d data set


%% 1-d data set
if info(2)==1
    
    
    if rank_yes==1
     temp_data=data(hot_index_x);
    [value,index]=sort(temp_data,'descend');
    hot_index_x=hot_index_x(index);
    end
    
    
    hotspot_data=data(hot_index_x);
    
else
    %% 2-d data
    if rank_yes==1
    temp_data=data(hot_index_x,hot_index_x);
    temp_data=max(abs(temp_data));
    [value,index]=sort(temp_data,'descend');
    hot_index_x=hot_index_x(index);
    end
    
    hotspot_data=data(hot_index_x,hot_index_x);

    % get hot spot entries


  

    if show_hotspot==1

        for i=1:length(hot_index_x)
            for j=1:length(hot_index_x)

               if ~(abs(hotspot_data(i,j)-expect_value)>range(1) &&  abs(hotspot_data(i,j)-expect_value)<range(2) )   
                     hotspot_data(i,j)=0;
               end

            end
        end

     end

    
    
end

 
    
    
end




        
        
    
