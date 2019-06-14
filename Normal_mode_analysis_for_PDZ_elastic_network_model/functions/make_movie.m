function make_movie(index,NormVector,coords,spring)
% spring_flag=0;
% cutoff_distance=4; 
index=index+6;
coord_info=size(coords);
N=coord_info(1);
temp_coord=zeros(coord_info(1),coord_info(2));
    
nframe=20;


vedio_name=cat(2,'mode_',num2str(index),'.avi');
delete vedio_name;
v = VideoWriter(vedio_name);
set(v,'FrameRate',10);
open(v);

% mov(1:nframe)= struct('cdata',[],'colormap',[]);
set(gca,'nextplot','replacechildren')


scaling=0.2;

temp=NormVector(:,index)';
    
for k=1:nframe     

    for i=1:N
    temp_coord(i,1)=coords(i,1)+k*scaling*temp(3*i-2);
    temp_coord(i,2)=coords(i,2)+k*scaling*temp(3*i-1);
    temp_coord(i,3)=coords(i,3)+k*scaling*temp(3*i);
    end
    f1=plot3(temp_coord(:,1),temp_coord(:,2),temp_coord(:,3),'.b');
       axis off;
   set(gcf,'color','white');
    set(f1,'linewidth',2,'MarkerSize',20)
          hold on;
        f1=plot3(temp_coord(:,1),temp_coord(:,2),temp_coord(:,3),'-r');
    set(f1,'linewidth',1)
          hold on;

   
      
    Plot_connection=0;
if Plot_connection==1
for i=1:N
    for j=1:N
        if spring(i,j)==1
            vector(1,:)=temp_coord(i,:); 
            vector(2,:)=temp_coord(j,:);
            plot3(vector(:,1),vector(:,2),vector(:,3),'-r');

        end
    end
end
end


    hold on;
    
    %%%%%%%%%

%     
%     for i=1:N
%     for j=1:N
%             
%             if spring_flag~=1 && sum((coords(i,:)-coords(j,:)).^2)<cutoff_distance^2 && i~=j
%                 vector(1,:)=temp_coord(i,:); 
%                 vector(2,:)=temp_coord(j,:);
%                 plot3(vector(:,1),vector(:,2),vector(:,3),'-.r');
%                 hold on;
%             end
%             
%             if spring_flag==1 && spring(i,j)==1
%                 vector(1,:)=temp_coord(i,:); 
%                 vector(2,:)=temp_coord(j,:);
%                 plot3(vector(:,1),vector(:,2),vector(:,3),'-.r');
%                 hold on;
%         
%         
%             
%             end
%     end
%     
%     end
    %%%%%%%%%%%
%     xlim([20 60]);
%     %xlabel('x');
%     ylim([40 80]);
%     zlim([20 50]);

     xlim([15 45]);
     xlabel('x');
     ylim([50 90]);
     zlim([20 60]);

    
    hold off;
    %mov(k)=getframe(gcf);  
    view(-137.5, 40);
    frame=getframe(gcf);
    %pause(1);
    writeVideo(v,frame);
end

close(v);
%movie2avi(mov, '1moviename.avi', 'compression', 'None');

