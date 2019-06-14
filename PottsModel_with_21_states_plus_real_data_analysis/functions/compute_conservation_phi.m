function [phi_conserv,compress_phi]=compute_conservation_phi(binary_sequence,pseudo_ratio,p,amino_N,plot_fig)

%%%% add pseudo count to the average residue abundance to compute the phi_conserv
[M_new,count]=size(binary_sequence);
mean_sequence_old=mean(binary_sequence,1);
mean_sequence_new=pseudo_ratio/amino_N+(1-pseudo_ratio)*mean_sequence_old;

L=count/amino_N;
% L=length(conserved_num);
% index_k=zeros(L,1);
% for k=2:L
%     index_k(k)=index_k(k-1)+conserved_num(k-1);
% end

phi_conserv=zeros(amino_N*L,1);
compress_phi=zeros(L,1);
for j=1:L
    temp=(j-1)*amino_N+1:j*amino_N;
    beta_l=mean_sequence_new(temp);
    temp_p=p;
    %temp_p=p(index_aminoacid(temp));
%      conservation=(1-beta_l).*log((1-beta_l)./(1-temp_p))+beta_l.*log(beta_l./temp_p);
    phi_conserv(temp)=-log((1-beta_l)./(1-temp_p))+log(beta_l./temp_p);
    compress_phi(j)=norm(phi_conserv(temp));
end


if plot_fig
    figureParameter
    f1=plot(1:length(phi_conserv),phi_conserv,'*k');
    a1=xlabel('$l$');
    a2=ylabel('$\phi_l$');
    %a2=ylabel('$1/\sigma_l^2$');
    %h1=legend('$\phi_l$');
    axis tight
    %ylim([0 6]);
    fig_name='./figure/weight.eps';
    figurePostTreat
end