function [conserved_num,index_k,index,remove_q]=remove_residue_q(binary_sequence_old,index_aminoacid_old,conserved_num_old,q,pseudo_ratio)
%%
% here,  we provide two method: 
% 1, a constant q [if q>0]
% 2, the least conserved aminoacide at a given site [if q=-1, which is the case for our paper]
% 3, if q is a vector, remove the residue at the corresponding site
%%
%q_constant=1;
%amino_N=21;
% p is the background rate of a given amino acid 

if ~exist('pseudo_ratio')
   pseudo_ratio=0;
end


L=length(conserved_num_old);
index_k_old=zeros(L,1);
for k=2:L
    index_k_old(k)=index_k_old(k-1)+conserved_num_old(k-1);
end

conserved_num=conserved_num_old-1;
index_k=zeros(L,1);
for k=2:L
    index_k(k)=index_k(k-1)+conserved_num(k-1);
end

index=[];
    
%% average mutation, plus pseudo count
[M_new,seq_L]=size(binary_sequence_old);
mean_sequence_old=mean(binary_sequence_old,1);

amino_N=21; 

mean_sequence_new=pseudo_ratio/amino_N+(1-pseudo_ratio)*mean_sequence_old;

%    conservation=zeros(1,amino_N*L);
phi_conserv=zeros(1,sum(conserved_num_old));


%% if a constant q, requires information from conserved_num_old

if length(q)==1

    if q>0

        for j=1:L

             %q=2;     %    throw away 
             index_relative1=1:q-1;
             index_relative2=q+1:conserved_num(j)+1;
             index_relative=cat(2,index_relative1,index_relative2);
             index=cat(2,index,index_k_old(j)+index_relative);
             remove_q=q;
             
            % to get phi_conserv 
%             temp=index_k_old(j)+1:index_k_old(j)+conserved_num_old(j);
%             beta_l=mean_sequence_new(temp);
%             temp_p=p(index_aminoacid_old(temp));
%             conservation=(1-beta_l).*log((1-beta_l)./(1-temp_p))+beta_l.*log(beta_l./temp_p);
%             phi_conserv(temp)=-log((1-beta_l)./(1-temp_p))+log(beta_l./temp_p);

        end

    else
        %% if selecting the least conserved sites
        % requires information from mean_sequence,index_aminoacid_old, back ground probability: p

        %% conservation measure
        %gap_p=0.1;
        remove_q=zeros(1,L);


        for j=1:L
%             temp=index_k_old(j)+1:index_k_old(j)+conserved_num_old(j);
%             beta_l=mean_sequence_new(temp);
%             temp_p=p(index_aminoacid_old(temp));
%             conservation=(1-beta_l).*log((1-beta_l)./(1-temp_p))+beta_l.*log(beta_l./temp_p);
%             phi_conserv(temp)=-log((1-beta_l)./(1-temp_p))+log(beta_l./temp_p);
            
            %% suppress the state that is most abundantly present
            %[C,I]=min(conservation);
            [C,I]=min(beta_l); % remove the least abundant residue, not recommended, used in the paper for the version before 2019-02
            %[C,I]=max(beta_l); % remove the most abundant residue  used after 2019-02, recommended
            index=cat(2,index,index_k_old(j)+cat(2,1:(I-1),(I+1):conserved_num_old(j)));
            remove_q(j)=index_aminoacid_old(index_k_old(j)+I);
        end

    end
    
else % if q is a vector, i.e., it indicates which residue to serve as the reference at each site
    
        for j=1:L

             %q=2;     %    throw away 
             index_relative1=1:q(j)-1;
             index_relative2=q(j)+1:conserved_num(j)+1;
             index_relative=cat(2,index_relative1,index_relative2);
             index=cat(2,index,index_k_old(j)+index_relative);
             
%             % to get phi_conserv 
%             temp=index_k_old(j)+1:index_k_old(j)+conserved_num_old(j);
%             beta_l=mean_sequence_new(temp);
%             temp_p=p(index_aminoacid_old(temp));
%             conservation=(1-beta_l).*log((1-beta_l)./(1-temp_p))+beta_l.*log(beta_l./temp_p);
%             phi_conserv(temp)=-log((1-beta_l)./(1-temp_p))+log(beta_l./temp_p);

        end    
        
        remove_q=q;
end


%phi_conserv=phi_conserv(index);
