function [corr_x1x2, to20_prediction_1, to20_prediction_2]=comparing_two_predictions(predict_1,name1,predict_2,name2,plot_fig)
% highlight the top 20 values for X1 or X2.

x1=abs(predict_1);
x2=abs(predict_2);

% red = [255 0 0]/255;
% dark_green = [0 180 0]/255;
% light_green=[0 255 0]/255;
% blue=[0 0 255]/255;
% grey=[17 17 17];

size_1=size(x1);
size_2=size(x2);
if sum(size_1==size_2)~=2
    error("The two vectors does not have the same length or one vector need to be transposed");
end

count=length(x1);

if count>20  
    cutoff=20; % the experimentally defined sector size for PDZ
else
    cutoff=count; % for artificial protein data, use its own length
    
end

[y1,index_1]=sort(x1,'descend');
[y2,index_2]=sort(x2,'descend');
index_sort_2=x2>=y2(cutoff);
index_sort_1=x1>=y1(cutoff);




if plot_fig
    %the vector has been normalized by their standard deviation
    figureParameter
    plot(1:count,x1/std(x1),'*r',1:count,x2/std(x2),'*b');
    xlabel("Residue index");
    legend(name1,name2);
    %ylim([0 3]);
    fig_name="./figure/Comparing_"+name1+"_"+name2+"_1.eps";
    figurePostTreat


    figureParameter
    plot(x1,x2,'*k',x1(index_sort_1),x2(index_sort_1),'*g',x1(index_sort_2),x2(index_sort_2),'*cyan',x1(index_sort_2 & index_sort_1),x2(index_sort_2 & index_sort_1),'*r');
    xlabel(name1);
    ylabel(name2);
    fig_name="./figure/Comparing_"+name1+"_"+name2+"_2.eps";
    %xlim([0 count]);
    figurePostTreat
end

corr_x1x2=corr(x1,x2);

to20_prediction_1=index_1(1:cutoff);
to20_prediction_2=index_2(1:cutoff);
