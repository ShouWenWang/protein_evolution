
%PDB_name_1='1be9.pdb'; reference_index_1=301;
PDB_name_1='1BFE.pdb'; reference_index_1=306;
cut_index_low=309;cut_index_high=398;

index=7; % second softest mode that is not zero


[Atom_type,coords]=read_PDB_data_C_beta(PDB_name_1,reference_index_1,cut_index_low,cut_index_high);




%index=7; % first softest mode that is not zero



plot_eigenmode(index,NormVector,coords,spring);
plot_mode_correlation(index,NormVector);



%figure,  plot3(coords(:,1),coords(:,2),coords(:,3),'*b'); hold on;

            
% 
% 
% Prob=NormVector(:,1)';
% 
% norm_Prob=zeros(Nx,Ny);
% for j=1:Ny
%    tp=Nx*(j-1);
%    norm_Prob(1:Nx,j)=Prob(tp+1:tp+Nx); 
% end
% 
% % figure, surf(X_index,Y_index,norm_Prob');
% % xlabel('x');
% % ylabel('y');
% % zlabel('Prob');
% 
% figureParameter
% f1=pcolor(X_index,Y_index,norm_Prob');
% a1=xlabel('$x$');
% a2=ylabel('$y$');
% fig_name='./figure/distribution.eps';
% figurePostTreat
% 
% 
% distr_x=zeros(1,Nx);
% for i=1:Nx
%     distr_x(i)=sum(norm_Prob(i,:));
% end
% distr_x=distr_x./sum(distr_x)*1/dx;
% 
% distr_y=zeros(1,Ny);
% for i=1:Ny
%     distr_y(i)=sum(norm_Prob(:,i));
% end
% distr_y=distr_y./sum(distr_y)*1/dy;
% 
% figure, plot(X_index,distr_x,'-r');
% xlabel('x');
% ylabel('Prob-x');
% 

%%


% vector(1,:)=[1 3 4];
% vector(2,:)=[4 5 6];
% vector(3,:)=[3 4 3];
% vector(4,:)=[4 2 1];
% 
% figure, plot3(vector(1:2,1),vector(1:2,2),vector(1:2,3),'-r');
%  hold on;
%  plot3(vector(3:4,1),vector(3:4,2),vector(3:4,3),'-r');

% 
% 
% figureParameter
% f1=pcolor(1:N,1:N,spring);
% a1=xlabel('index');
% a2=ylabel('index');
% fig_name='./springnetwork.eps';
% figurePostTreat
% 
% 
% 
% figureParameter
% f1=pcolor(1:N,1:N,spring);
% a1=xlabel('index');
% a2=ylabel('index');
% fig_name='./springnetwork.eps';
% figurePostTreat
% 
% 
% 
% hist(distance(1,:))
% 
% 
% 
% Z=[[1 1  -1 -1];[1 1 -1 -1];[-1 -1 1 1];[-1 -1 1 1]];
% [NormVector,orderEigValue]=orderedEigSystem(Z);
% 
% 
