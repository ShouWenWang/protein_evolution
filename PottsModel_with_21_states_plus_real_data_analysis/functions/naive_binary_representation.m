function [binary_sequence,inverse_index,conserved_num,index_aminoacid]=naive_binary_representation(efa1,amino_N)

%amino_N=21;
clear_zero_state=0; % by default, we do not do so. This step is being taken care by adding pseudo count now. 
new_efa1=efa1;
[M_new, L]=size(new_efa1);
L1=L*amino_N;


binary_sequence=zeros(M_new,L1);
for j=1:M_new
    for k=1:L
        binary_sequence(j,(k-1)*amino_N+new_efa1(j,k))=1;
    end
end

%%
index_aminoacid=zeros(L1,1);
inverse_index=zeros(L1,1);

for k=1:L
   inverse_index((k-1)*amino_N+1:k*amino_N)=k;
   index_aminoacid((k-1)*amino_N+1:k*amino_N)=1:amino_N;
end

mean_sequence=mean(binary_sequence,1);
%figure, bar(1:L1,mean_sequence,'r'); % this is not a very good measure of importance

%%
conserved_num=zeros(L,1);
if clear_zero_state
    for k=1:L
        site_mean=mean_sequence((k-1)*amino_N+1:k*amino_N);
        conserved_num(k)=sum(site_mean>0); 
        
        
            for j=1:amino_N
                if site_mean(j)==0
                    index_aminoacid((k-1)*amino_N+j)=0; % if never exists, reset to zero
                end
            end
    end
    
else
    conserved_num=zeros(L,1)+amino_N;
    
end

%% conservation
% Prob_x=zeros(1,amino_N); % x is the label for a particular amino acid   
% 
% for x=1:amino_N
%     
%    index=x:amino_N:L1;
%    Prob_x(x)=sum(mean_sequence(index));
%    
% end
% Prob_x=Prob_x/sum(Prob_x);
% 
% pseudocount=0.0001; % a sufficiently small number to avoid divergence
% KL_diverg_1=zeros(L,1);
% KL_diverg_2=zeros(L,1);
% for k=1:L
% Prob_x_k=mean_sequence((k-1)*amino_N+1:k*amino_N)+pseudocount; % the frequency of x in a given site k of the sequence 
% Prob_x_k=Prob_x_k/sum(Prob_x_k);
% 
% KL_diverg_1(k)=-sum((Prob_x).*log(Prob_x_k./Prob_x));
% KL_diverg_2(k)=sum((Prob_x_k).*log(Prob_x_k./Prob_x));
% end
% 
% 
% sector=[323 324 325 327 328 329 330 336 338 341 347 353 359 362 367 372 375 376 379 388]-312;
% 
% 
% figureParameter,
% bar(312+(1:L),[KL_diverg_1,KL_diverg_2])
% a2=ylabel('Conservation');
% fig_name='./figure/conserv.eps';
% figurePostTreat
% 
% 
% hot_index_x=find(KL_diverg_2>1.3);
% (312+hot_index_x)'



