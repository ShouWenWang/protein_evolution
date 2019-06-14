function compres_inv_C=matrix_compression(index_k,conserved_num,inv_matrix,zero_sum_gauge)
%%
% inv_matrix, the matrix to be compressed
% index_k: when plus 1, it is the starting index of the site k
% conserved_num: the number of available kinds of amino acid, within each site 
%  Below, we take Frobenius norm of the matrix within each block of 21*21
%%


L=length(index_k);
compres_inv_C=zeros(L,L);

%zero_sum_gauge=1;
if zero_sum_gauge==1
    for j=1:L
        temp_index_1=index_k(j)+1:index_k(j)+conserved_num(j);
        S1=length(temp_index_1);
      for k=1:L
        temp_index_2=index_k(k)+1:index_k(k)+conserved_num(k);
        S2=length(temp_index_2);
        data0=zeros(S1+1,S2+1);
        data0(1:S1,1:S2)=inv_matrix(temp_index_1,temp_index_2);
        mean_1=mean(data0,1);
        mean_2=mean(data0,2);
        mean_tot=mean(mean_1);
        for n1=1:S1+1
            for n2=1:S2+1
                data0(n1,n2)=data0(n1,n2)-mean_1(n2)-mean_2(n1)+mean_tot; 
            end
        end
        
        tot_N=(S1+1)*(S2+1);
        compres_inv_C(j,k)=sqrt(sum(sum(data0.^2))/tot_N);

            %compres_inv_C(j,k)=max(max(abs(data0)));
      end
    end
    
else
    
    for j=1:L
        temp_index_1=index_k(j)+1:index_k(j)+conserved_num(j);
      for k=1:L
        temp_index_2=index_k(k)+1:index_k(k)+conserved_num(k);
            data0=inv_matrix(temp_index_1,temp_index_2);
            tot_N=length(temp_index_1)*length(temp_index_2);
            compres_inv_C(j,k)=sqrt(sum(sum(data0.^2))/tot_N);

            %compres_inv_C(j,k)=max(max(abs(data0)));
      end
    end
end
