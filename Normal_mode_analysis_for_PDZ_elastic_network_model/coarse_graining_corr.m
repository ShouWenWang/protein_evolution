function coarse_graining_corr(correlation,mesh_size)

info=size(correlation);
N=info(1);
%mesh_size=4;

coarse_corr_size=ceil(N/mesh_size);
coarse_grained_corr=zeros(coarse_corr_size,coarse_corr_size);
for i=1:coarse_corr_size
    for j=1:coarse_corr_size
        x_low=(i-1)*mesh_size+1;
        y_low=(j-1)*mesh_size+1;

        x_high=i*mesh_size;
        y_high=j*mesh_size;
        if x_high>N
           x_high=N;
        end
        
        if y_high>N
           y_high=N;
        end
        
        area=(x_high-x_low)*(y_high-y_low);
        
        coarse_grained_corr(i,j)=sum(sum(correlation(x_low:x_high,y_low:y_high)))/area;
                
    end
end


figureParameter
f1=pcolor(1:coarse_corr_size,1:coarse_corr_size,coarse_grained_corr);
colorbar;
set(gca, 'clim', [-0.05 0.05]);
a1=xlabel('$x$');
a2=ylabel('$y$');
fig_name='./figure//mode_correlation.eps';
figurePostTreat