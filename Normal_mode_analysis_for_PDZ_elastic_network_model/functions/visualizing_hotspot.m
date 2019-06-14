function Amino_acid_loc=visualizing_hotspot(coords,beta_atom_index,hotspot_index,spring,reference_index,Atom_type)
%%
% hotspot_index is the index within the set of Cbeta atoms,
% we need to use the index information for the C beta atom 
% to put the hotspot in the context of the original structure
%
% beta_atom_index is the index of beta atom in the original structure
%%

figureParameter, plot3(coords(:,1),coords(:,2),coords(:,3),'-*b')
xlim([25 50]);
ylim([40 80]);
zlim([20 50]);

hold on

index=beta_atom_index(hotspot_index);
N=length(index);

Amino_acid_loc=zeros(N,1);
for i=1:N
    Amino_acid_loc(i)=sum(Atom_type(1:index(i)))/(-3)+reference_index-1;
end



plot3(coords(index,1),coords(index,2),coords(index,3),'or')
hold on;

view(-137.5, 40);

fig_name='./figure/coordinates.eps';
figurePostTreat

  plot_edge=1;
  %% plot the edge
  if plot_edge==1
      vector=zeros(2,3);

      for i=1:N
         for j=1:N
             if spring(index(i),index(j))==1
                vector(1,:)=coords(index(i),:);
                vector(2,:)=coords(index(j),:);

                plot3(vector(:,1),vector(:,2),vector(:,3),'-r')
             end
         end
      end
  
  end
  
hold off

