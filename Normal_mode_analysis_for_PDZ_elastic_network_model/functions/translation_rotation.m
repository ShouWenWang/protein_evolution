function coords=translation_rotation(coords,target_Calpha,rotation_on)

%% translational correction 
for i=1:3
coords(:,i)=coords(:,i)-mean(coords(:,i));
end

%rotation_on=0;


if rotation_on==1
    %% rotational operation,  make the first point has y=0 and z=0
    old_coords=coords;

    info=size(coords);
    N=info(1);


    %% rotate along z axis, making y=0
    vec_1=coords(target_Calpha,1)+sqrt(-1)*coords(target_Calpha,2);
    angle_z=angle(vec_1); % the rotation angle along z axis

    for i=1:N
          r0=sqrt(coords(i,1)^2+coords(i,2)^2);
          angle_0=angle(coords(i,1)+sqrt(-1)*coords(i,2));
          angle_1=angle_0-angle_z;
          coords(i,1)=r0*cos(angle_1);
          coords(i,2)=r0*sin(angle_1);
    end

    %% rotate along y axis,  making y=0 and z=0

    vec_2=coords(target_Calpha,1)+sqrt(-1)*coords(target_Calpha,3);
    angle_y=angle(vec_2); % the rotation angle along z axis

    for i=1:N
          r0=sqrt(coords(i,1)^2+coords(i,3)^2);
          angle_0=angle(coords(i,1)+sqrt(-1)*coords(i,3));
          angle_1=angle_0-angle_y;
          coords(i,1)=r0*cos(angle_1);
          coords(i,3)=r0*sin(angle_1);
    end

    %% plot

%     figureParameter
%     f1=plot3(old_coords(:,1),old_coords(:,2),old_coords(:,3),'-*b');
%     a1=xlabel('x');
%     a2=ylabel('y');
%     a3=zlabel('z');
%     hold on;
%     plot3(coords(:,1),coords(:,2),coords(:,3),'-*r')
%     hold off;
%     fig_name='./conform_rotation.eps';
%     figurePostTreat

end
