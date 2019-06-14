function [seq_corr,mean_sequence]=dealing_Protein_data_reduced_description(efa1)

amino_N=21;
s=size(efa1);
count=s(2);
L1=L*amino_N;
M_new=s(1);

binary_sequence=zeros(M_new,L1);
for j=1:s(1)
    for k=1:L
        binary_sequence(j,(k-1)*amino_N+efa1(j,k))=1;
    end
end


%figure,hist(sum_rule,20)

mean_sequence0=mean(binary_sequence,1);

%% find the most frequent K amino acid type at each site

% most conserved m amino acids
cons_N=2;
conserved_info=zeros(count,cons_N);

    for k=1:count
        site_mean=mean_sequence0((k-1)*amino_N+1:k*amino_N);
        temp=[site_mean',(1:length(site_mean))'];
        temp2=sortrows(temp);
        for j=1:cons_N
        conserved_info(k,j)=temp2(amino_N+1-j,2);
        end
    end
    
  %% generating reduced description
  new_sequence=zeros(M_new,count);
  
for j=1:M_new
    for k=1:L
        ty=1; % default mutant
        for l=1:cons_N
            if efa1(j,k)==conserved_info(k,l)
                ty=(l-1)/cons_N;
            end
        end
        new_sequence(j,k)=ty;
    end
end


mean_sequence=mean(new_sequence,1);
%%


%figure,hist(mean_sequence);
plot_distribution(mean_sequence,[0 1],2,'$\langle \beta_l\rangle_*$');

figureParameter
bar(1:count,mean_sequence,'r');
fig_name='./figure/corr_mode.eps';
a1=ylabel('$\langle \beta_l\rangle*$');
%set(h1,'location','south')
xlim([0 80]);
%set(gca,'YTICKLABEL',[-0.1 0 0.1 0.2 0.3 0.4 0.5]);
figurePostTreat


%[hot_index_x,hot_index_y,hotspot_data]=search_hotspot(mean_sequence,0.5,[0.12 1]);


%% covariance analysis
seq_corr=zeros(count,count);


for k=1:M_new
    
%     for i=1:count
%         for j=1:count
% 
%            seq_corr(i,j)=seq_corr(i,j)+(binary_sequence(k,i)-mean_sequence(i))*(binary_sequence(k,j)-mean_sequence(j));
% 
%         end
%     end

      %% this is faster     
       seq_corr=seq_corr+(new_sequence(k,:)-mean_sequence)'*(new_sequence(k,:)-mean_sequence);

end



seq_corr=seq_corr./M_new;



%%


%% analyze the correlation matrix
% histogram without the self-correlation
seq_corr_off_diag=seq_corr;
seq_diag=zeros(count,2);
for i=1:count
    seq_corr_off_diag(i,i)=0;
    seq_diag(i,1)=i;
    seq_diag(i,2)=seq_corr(i,i);
end

corr_vector=zeros(count*count,1);
for i=1:count
  corr_vector((i-1)*count+1:i*count)=seq_corr_off_diag(:,i);
end

new_corr_vector=corr_vector(abs(corr_vector)>=0 );%& abs(corr_vector)<0.24);
    
%figure, hist(new_corr_vector,100);
plot_distribution(new_corr_vector,[-0.25 0.25],40,'$C_{ll''}^-$');


figureParameter
bar(1:count,seq_diag(:,2),'r');
fig_name='./figure/corr_mode_var.eps';
%set(gca,'XTICKlabel',index2);
a1=ylabel('$\sigma_l^*$');
%set(h1,'location','south')
xlim([0 count]);
%set(gca,'YTICKLABEL',[-0.1 0 0.1 0.2 0.3 0.4 0.5]);
figurePostTreat


%%  eigenvalue analysis of Corr
%select_data=seq_corr_off_diag;
select_data=seq_corr;

[NormVector_corr,orderEigValue_corr]=orderedEigSystem(select_data,0);


figureParameter
f1=plot(1:count,orderEigValue_corr(:,1),'*r');
a1=xlabel('mode');
xlim([1 count])
a2=ylabel('$\lambda$');
fig_name='./figure/corr_eig.eps';
figurePostTreat


%% absolute version
figureParameter
bar(1:count,[abs(NormVector_corr(:,1)),abs(NormVector_corr(:,2))]);
%set(gca,'XTICKlabel',index2);
h1=legend('$|\nu_{1}^l|$','$|\nu_{2}^l|$');%-\Delta_6^l\;$');
set(h1,'location','northwest')
xlim([0 count]);
ylim([0 0.3])
fig_name='./figure/corr_mode2.eps';
figurePostTreat

figureParameter
bar(1:count,[abs(NormVector_corr(:,end-1)),abs(NormVector_corr(:,end))]);
%set(gca,'XTICKlabel',index2);
h1=legend('$|\nu_{end-1}^l|$','$|\nu_{end}^l|$');%-\Delta_6^l\;$');
set(h1,'location','north')
xlim([0 count]);
ylim([0 0.3])
fig_name='./figure/corr_mode3.eps';
figurePostTreat


%% original version
figureParameter
bar(1:count,[NormVector_corr(:,1),NormVector_corr(:,3),NormVector_corr(:,4),NormVector_corr(:,5)]);
%set(gca,'XTICKlabel',index2);
h1=legend('$\nu_1^l$','$\nu_3^l$','$\nu_{4}^l$','$\nu_{5}^l$');%-\Delta_6^l\;$');
set(h1,'location','south')
xlim([0 count]);
%ylim([-0.2 0.2])
fig_name='./figure/corr_mode4.eps';
figurePostTreat

figureParameter
bar(1:count,[real(NormVector_corr(:,45)),imag(NormVector_corr(:,45))]);
%set(gca,'XTICKlabel',index2);
h1=legend('$Re(\nu_{45}^l)$','$Imag(\nu_{45}^l)$');%-\Delta_6^l\;$');
set(h1,'location','southwest')
xlim([0 count]);
ylim([-0.3 0.3])
fig_name='./figure/corr_mode5.eps';
figurePostTreat

