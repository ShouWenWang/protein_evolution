function plot_eigenmode(index,NormVector,coords,spring)
coord_info=size(coords);
N=coord_info(1);

figureParameter, f1=plot3(coords(:,1),coords(:,2),coords(:,3),'.r');
hold on;
 plot3(coords(:,1),coords(:,2),coords(:,3),'-k');
hold on;

Plot_connection=0;
if Plot_connection==1
for i=1:N
    for j=1:N
        if spring(i,j)==1
            vector(1,:)=coords(i,:); 
            vector(2,:)=coords(j,:);
            plot3(vector(:,1),vector(:,2),vector(:,3),'-r');
            hold on;
        end
    end
end
end

scaling=10;

for i=1:N
    temp=scaling*NormVector(3*(i-1)+1:3*i,index)';
    if sqrt(sum(temp.^2))>0.5
        vector(1,:)=coords(i,:); 
        vector(2,:)=coords(i,:)+temp;
        plot3(vector(:,1),vector(:,2),vector(:,3),'-b');
        hold on;
    end
end

axis off;
hold off;
%view(-37.5, 40);
%view(-160, -60);
%view(45,80);
%view(-137.5, 40);
view(-130, 10);
   set(gcf,'color','white');
fig_name='./figure/eigen_mode.eps';
figurePostTreat

