function compres_vector=vector_compression(index_k,conserved_num,vector)
%%
% Implement zero-sum gauge directly
%  Below, we take Frobenius norm of the matrix within each block of 21*21
%%


L=length(index_k);
compres_vector=zeros(1,L);


    for j=1:L
        temp_index_1=index_k(j)+1:index_k(j)+conserved_num(j);
        
        S1=length(temp_index_1);
        data0=zeros(S1,1);
        data0(1:S1)=vector(temp_index_1);
        
%         mean_1=mean(data0,1);
%             for n1=1:S1+1
%                 data0(n1,1)=data0(n1,1)-mean_1; 
%             end
     
        
%       tot_N=(S1+1);
        compres_vector(j)=sqrt(sum(data0.^2)/S1);

            %compres_vector(j,k)=max(max(abs(data0)));

    end

 
