function [new_binary_sequence,inverse_index,conserved_num,index_aminoacid]=reduced_representation(efa1,amino_N)

% this is no longer used, updated 2019-03-02

% generate a reduced binary representation that eliminates sites that do
% not occur at all 
% critical_N: critical number of counts needed to be represented
%%
% number2letter_map(1) = '-';
% number2letter_map(2) = 'A';
% number2letter_map(3) = 'C';
% number2letter_map(4) = 'D';
% number2letter_map(5) = 'E';
% number2letter_map(6) = 'F';
% number2letter_map(7) = 'G';
% number2letter_map(8) = 'H';
% number2letter_map(9) = 'I';
% number2letter_map(10) = 'K';
% number2letter_map(11) = 'L';
% number2letter_map(12) = 'M';
% number2letter_map(13) = 'N';
% number2letter_map(14) = 'P';
% number2letter_map(15) = 'Q';
% number2letter_map(16) = 'R';
% number2letter_map(17) = 'S';
% number2letter_map(18) = 'T';
% number2letter_map(19) = 'V';
% number2letter_map(20) = 'W';
% number2letter_map(21) = 'Y';

%%

critical_N=0;
new_efa1=efa1;
s=size(new_efa1);
L=s(2); % number of sites in the original representation
L1=L*amino_N;
M_new=s(1);

%% naive binary representation
binary_sequence=zeros(M_new,L1);
for j=1:s(1)
    for k=1:L
        binary_sequence(j,(k-1)*amino_N+new_efa1(j,k))=1;
    end
end

mean_sequence=mean(binary_sequence,1);
sum_sequence=sum(binary_sequence,1);

%% critical number
%critical_N=1; % only sites which have counts larger this critical number is included
%%

    
%% rank the  amino acid type at each site  according to their appearance frequence
% 
% % most conserved m amino acids

conserved_info=zeros(L,amino_N);
conserved_num=zeros(L,1);
sort_by_conservation=0;

    for k=1:L
        site_sum=sum_sequence((k-1)*amino_N+1:k*amino_N);
        
        temp=1:amino_N; 
        temp=temp(site_sum>=critical_N); 
        
        if sort_by_conservation==1
        temp1=[site_sum(site_sum>=critical_N)',temp'];
        temp1=sortrows(temp1);
        temp2=temp1(:,2);
        else 
            temp2=temp;
        end
        
        conserved_num(k)=length(temp2);
        conserved_info(k,1:conserved_num(k))=temp2; 

        %temp=[site_mean',(1:amino_N)'];
%         temp2=sortrows(temp);
% %        appearance_info(k,:)=temp2(:,1);
%         
%         s=0;
%         while temp2(s+1)==0
%             s=s+1;
%         end
        
%         for j=1:amino_N-s
%         conserved_info(k,j)=temp2(amino_N+1-j,2); % ranked amnio acid preference at each site
%         end
                
        % rank according to the alphabeta order
%        conserved_info(k,1:amino_N-s)=sort(conserved_info(k,1:amino_N-s));
%        conserved_num(k)=amino_N-s;

       
        
    end
    
    %% the initial index of site k in the new representation
    index_k=zeros(L,1); % the starting index at kth amino acid in the new representation

    for k=2:L
        index_k(k)=index_k(k-1)+conserved_num(k-1);
    end
    
    %% the total sequence site for the new representation
  count=sum(conserved_num);
    
  %% generating reduced description
  new_binary_sequence=zeros(M_new,count);
  
    for j=1:M_new
        for k=1:L
            for l=1:conserved_num(k)
                if new_efa1(j,k)==conserved_info(k,l)
                     index=index_k(k)+l;
                     new_binary_sequence(j,index)=1;
                end
            end
        end
    end

    %% inverse indexing: know index in the new representation, 
    % generate the index for the original representation
    
    inverse_index=zeros(count,1);
    temp=0;
    for i=1:L
        inverse_index(temp+1:temp+conserved_num(i))=i;
        temp=temp+conserved_num(i);
    end
    
    
%     for i=1:count
%         for k=1:L-1
%             if i>=index_k(k) && i<index_k(k+1)
%                inverse_index(i)=k; 
%             end
%             
%         end
%            
%            if i>=index_k(L)
%                  inverse_index(i)=L; 
%            end
%                
%     end
    
    
    
    
    index_aminoacid=[];

    for k=1:L

    index_aminoacid=cat(2,index_aminoacid,conserved_info(k,1:conserved_num(k)));

    end

    
    
    