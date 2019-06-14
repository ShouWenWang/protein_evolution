function correct_hits_N20=prediction_analysis(data,reverse_index,ranked_true_hot_index,starting_index,plot_fig)
%% data should be one-dimensional. We only care about the magnitude of data
%  hot_spot_index should be in terms of the index for this truncated
%  sequences, i.e., minus 312, which is the starting_index
% reverse_index refers to the index for data,  which is truncated from the
% original sequences. 

%% 
L=length(data);

if L>40
    range_cutoff=40;
else
    range_cutoff=L;
end



comp_index=ranked_true_hot_index-starting_index;
L1=length(comp_index);
%index_org=1:L;
[value_new,index]=sort(abs(data),'descend');
rank_index=reverse_index(index);

count=0;
correct_ratio=zeros(1,L);
correct_N=zeros(1,L);
temp=zeros(L1,1);
correct_matrix=zeros(L1,L);
for j=1:L
     temp_index=find(rank_index(j)==comp_index);
       if temp_index  % true if there is a match; otherwise false
           count=count+1;
           temp(temp_index)=1; 
       end
    correct_matrix(:,j)=temp;
    correct_N(j)=count;
    correct_ratio(j)=count/j; % the correct_ratio
end

if L>=20
    correct_hits_N20=correct_N(20);
else
    correct_hits_N20=correct_N(L);
end

disp("Correct hits for top20 predictions: "+num2str(correct_hits_N20));

if plot_fig
    
    figureParameter
    f1=plot(1:L,correct_ratio,'*r');
    xlabel('Count');
    ylabel('Correct ratio');
    xlim([1,range_cutoff]);
    ylim([0,1]);
    fig_name='./figure/correct_ratio.eps';
    figurePostTreat

    figureParameter
    f1=plot(1:L,correct_N,'*r');
    xlabel('Count');
    ylabel('Correct #');
    xlim([1,range_cutoff]);
    %ylim([0,1]);
    fig_name='./figure/correct.eps';
    figurePostTreat


    figureParameter
    pcolor(1:range_cutoff,1:L1,abs(correct_matrix(:,1:range_cutoff)));
    colorbar;
    fig_name='./figure/hot_matrix.eps';
    figurePostTreat

end
