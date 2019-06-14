function reverse_index=reverse_index_from_logic_input(forward_index)
%% 
% forward index is a string of zeros and ones, obtained from logic
% operation. This generates a truncation of the original data that satisfy this
% condition.
% Here, we want to get the original index in the new representation
%%

L=length(forward_index);
L1=sum(forward_index);
reverse_index=zeros(1,L1);
count=0;
for j=1:L
    if forward_index(j)==1
       count=count+1; 
       reverse_index(count)=j;
    end
end


