function coords=read_PDB_data_C_alpha(PDB_name,reference_index,cut_index_low,cut_index_high)


%load the pdb file of PDZ
pdbstruct=pdbread(PDB_name);
%should use chain A - this is the actual protein. There is
%another chain, chain B, but it's just a small peptide binding to the
%protein.

% Get the first model in pdbstruct - CONTAINS ATOM COORDS
modelstruct = pdbstruct.Model(1);

% get the rows of modelstruct.Atom that contain alpha carbon atoms AND that are in chain A (binary array with 1 for those lines and 0 for others)
relevantAtom=zeros(size(modelstruct.Atom(:),1),1);


for i=1:size(modelstruct.Atom(:),1) %loop over atoms
    atom_name_temp=char(modelstruct.Atom(i).AtomName);
    if length(atom_name_temp)==2  && char(modelstruct.Atom(i).chainID)==char('A') &&  atom_name_temp(1)==char('C') &&  atom_name_temp(2)==char('A') 
        relevantAtom(i,1)=1;
            
    end
end


%% sequence selection based on Calpha chain
% erase the higher branch
count=0;

%cut_index_1=309;

for i=1:size(modelstruct.Atom(:),1)
    if relevantAtom(i,1)==1
        count=count+1;
         
        % this should be computed within the (larger) if condition,
        % otherwise, i will continue change with a fixed count, resulting
        % in a shift of low_index
        if count==cut_index_low-reference_index+1
             low_index=i;
        end
    
        if count==cut_index_high-reference_index+1
            high_index=i;
        end
    end
    
    
end
    
relevantAtom(1:low_index-1,1)=0;
relevantAtom(high_index+1:end,1)=0;




% Create Nx3 matrix of 3D coords of the alpha carbons of PDZ chain A

%% coordinate selection

% with both C alpha and C beta coordinate
coords = [modelstruct.Atom(logical(relevantAtom(:,1))).X;
modelstruct.Atom(logical(relevantAtom(:,1))).Y;
modelstruct.Atom(logical(relevantAtom(:,1))).Z; ]';


%just to visualize







