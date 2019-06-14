function [NormVector,orderEigValue,coord_normVector,spring,Hessen_2d]=normal_mode_computation_Cbeta_mutation(coords,Atom_type,cutoff_distance,factor)

chain_info=size(coords);
N=chain_info(1);
distance=zeros(N,N);
spring=zeros(N,N);
Hessen=zeros(N,N,3,3);
%Hessen_1=zeros(N,N);
Hessen_2d=zeros(3*N,3*N);

%factor=0.8;

%% off diagonal terms
for i=1:N
    for j=1:N
        distance(i,j)=sqrt(sum((coords(i,:)-coords(j,:)).^2)); 
        if distance(i,j)<cutoff_distance && i~=j
            
         
            
            %%  determination of the spring_strength
            if Atom_type(i)==0 && Atom_type(j)==0  % Cbeta-Cbeta
             spring_strength=0.5; 
            end
            
            if (Atom_type(i)==1 && Atom_type(j)==0) || (Atom_type(i)==0 && Atom_type(j)==1)
             spring_strength=0.5*factor;   % single mutation of Cbeta
            end
            
            if Atom_type(i)==1 && Atom_type(j)==1
             spring_strength=0.5*factor^2;   %double mutation
            end
            
            if Atom_type(i)==-3 && Atom_type(j)==-3 && abs(i-j)==1 % close by C alpha atoms
                   spring_strength=2; %C alpha,  contact along the chain
            end
            
            if Atom_type(i)==-3 && Atom_type(j)==-3 && abs(i-j)~=1
                    spring_strength=1;  % C alpha, distant contact
               
            end
            
            
            if (Atom_type(i)==-3 && Atom_type(j)==0) ||  (Atom_type(i)==0 && Atom_type(j)==-3)
                    spring_strength=1; %C alpha-C beta
               
            end
            
            if (Atom_type(i)==-3 && Atom_type(j)==1) ||  (Atom_type(i)==1 && Atom_type(j)==-3)
                    spring_strength=1*factor; %C alpha- mutation of C beta
               
            end
            
            
            %%
            
            for l=1:3
                for s=1:3
                  Hessen(i,j,l,s)=-spring_strength*(coords(j,l)-coords(i,l))*(coords(j,s)-coords(i,s))/(distance(i,j)^2);
                end
            end
            
           
            % plot
            spring(i,j)=1;

        end
    end
end
 
%% diagonal term 

for i=1:N
    Hessen(i,i,l,s)=0;
    for j=1:N
         if distance(i,j)<cutoff_distance && i~=j
            
             %% determination of the spring_strength
            if Atom_type(i)==0 && Atom_type(j)==0
             spring_strength=0.5; 
            end
            
            if (Atom_type(i)==1 && Atom_type(j)==0) || (Atom_type(i)==0 && Atom_type(j)==1)
             spring_strength=0.5*factor; 
            end
            
            if Atom_type(i)==1 && Atom_type(j)==1
             spring_strength=0.5*factor^2; 
            end
            
            if Atom_type(i)==-3 && Atom_type(j)==-3 && abs(i-j)==1 % close by C alpha atoms
                   spring_strength=2;
            end
            
            if Atom_type(i)==-3 && Atom_type(j)==-3 && abs(i-j)~=1
                    spring_strength=1;
               
            end
            
            
            if (Atom_type(i)==-3 && Atom_type(j)==0) ||  (Atom_type(i)==0 && Atom_type(j)==-3)
                    spring_strength=1;
               
            end
            
            if (Atom_type(i)==-3 && Atom_type(j)==1) ||  (Atom_type(i)==1 && Atom_type(j)==-3)
                    spring_strength=1*factor;
               
            end
            
            %%
            for l=1:3
                for s=1:3
                 Hessen(i,i,l,s)=Hessen(i,i,l,s)+spring_strength*(coords(j,l)-coords(i,l))*(coords(j,s)-coords(i,s))/(distance(i,j)^2);
                end
            end
        end
    end
end

%% matrix transform to 2-d
for i=1:N
    for j=1:N
        for l=1:3
            for s=1:3
             Hessen_2d(3*(i-1)+l,3*(j-1)+s)=Hessen(i,j,l,s);
            end
        end
    end
end
    

%% eigenvalues,  eigenmodes
[NormVector,orderEigValue]=orderedEigSystem(Hessen_2d,0);

coord_normVector=zeros(N,3*N,3);

for i=1:3*N
    for j=1:N
        coord_normVector(j,i,1)=NormVector(3*j-2,i);
        coord_normVector(j,i,2)=NormVector(3*j-1,i);
        coord_normVector(j,i,3)=NormVector(3*j,i);
    end
end
