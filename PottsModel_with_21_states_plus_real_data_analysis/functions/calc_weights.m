function Weight=calc_weights(seq,similarity_value_threshold)
%% Weighting each sequence according to how different they are from other sequences. Different sequences get more weight
% for the PDZ data (24000*79), it takes 600s, which is 10 mins
%%


[M,L]=size(seq);
W=zeros(M,1); % the inverse weight
for j=1:M
    W(j)=W(j)+1;
    for k=j+1:M
       id=0; % number of identical residues for two sequences
        for n=1:L
            if seq(j,n)==seq(k,n)
                id=id+1;
            end
        end
        
        if id<=similarity_value_threshold*L
            W(j)=W(j)+1;
            W(k)=W(k)+1;
        end
        
    end

end

Weight=1./W;
