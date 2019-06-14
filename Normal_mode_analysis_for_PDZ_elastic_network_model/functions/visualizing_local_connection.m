function visualizing_local_connection(coords,beta_atom_index,target_index,spring,Atom_type)
%%
% hotspot_index is the index within the set of Cbeta atoms,
% we need to use the index information for the C beta atom 
% to put the hotspot in the context of the original structure
%
% beta_atom_index is the index of beta atom in the original structure
%%


figure, plot3(coords(:,1),coords(:,2),coords(:,3),'-*k')
hold on

info=size(coords);



index=beta_atom_index(target_index);
%target_coord=coords(index,:);

plot3(coords(index,1),coords(index,2),coords(index,3),'sb')
hold on

for i=1:info(1)

    if spring(i,index)==1 && Atom_type(i)==1 % connection with C alpha
     plot3(coords(i,1),coords(i,2),coords(i,3),'or')
     hold on
    end
    
     if spring(i,index)==1 && Atom_type(i)==0 % connection with C beta
     plot3(coords(i,1),coords(i,2),coords(i,3),'og')
     hold on
    end
end
hold off
    
