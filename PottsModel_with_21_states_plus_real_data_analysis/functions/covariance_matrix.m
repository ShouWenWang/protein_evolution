function [seq_corr,mean_sequence]=covariance_matrix(binary_sequence,Weight,pseudo_ratio,site_size)

%reweighting_threshold=0.1/21;
if pseudo_ratio==0
    
    
[M_new,count]=size(binary_sequence);

%reweighting_threshold=0.1  for the artificial binary sequences
%reweighting_threshold=0.1/21  for the protein sequences (the summation
%within each site is guaranteed to be 1, which generates a lot of spurious
%similarlity


%Weight=calc_weights(binary_sequence,reweighting_threshold);

%count=length(binary_sequence(1,:));
%mean_sequence=mean(binary_sequence,1);

temp_1=zeros(count,1);
temp_2=zeros(count,count);

for j=1:count
temp_1(j)=Weight'*binary_sequence(:,j);
end
mean_sequence=temp_1/sum(Weight);

% compute the upper triangle: off diagonal terms
for j=1:count
    for k=j+1:count
        temp_2(j,k)=(Weight.*binary_sequence(:,j))'*binary_sequence(:,k);
    end
end
% the lower triangle terms
temp_2=temp_2+temp_2';
% the diagonal terms
for j=1:count
    temp_2(j,j)=(Weight.*binary_sequence(:,j))'*binary_sequence(:,j);
end

joint_prob=temp_2/sum(Weight);



seq_corr=joint_prob-mean_sequence*mean_sequence';

else
    %% only for protein data
    
%pseudo_count=0.5; % pseudo count
%% assuming that every kind of amino acid are observed in each site
%amino_N=20; % after removing one residue 
amino_N=21;
pseudo_count=sum(Weight)*pseudo_ratio/(1-pseudo_ratio);

[M_new,count]=size(binary_sequence);

%reweighting_threshold=0.1  for the artificial binary sequences
%reweighting_threshold=0.1/21  for the protein sequences (the summation
%within each site is guaranteed to be 1, which generates a lot of spurious
%similarlity


%Weight=calc_weights(binary_sequence,reweighting_threshold);

%count=length(binary_sequence(1,:));
%mean_sequence=mean(binary_sequence,1);

temp_1=zeros(count,1);
temp_2=zeros(count,count);

for j=1:count
temp_1(j)=Weight'*binary_sequence(:,j);
end
mean_sequence=(pseudo_count/amino_N+temp_1)/(pseudo_count+sum(Weight));

% compute the upper triangle: off diagonal terms
for j=1:count
    for k=j+1:count
        temp_2(j,k)=(Weight.*binary_sequence(:,j))'*binary_sequence(:,k);
    end
end
% the lower triangle terms
temp_2=temp_2+temp_2';
% the diagonal terms
for j=1:count
    temp_2(j,j)=(Weight.*binary_sequence(:,j))'*binary_sequence(:,j);
end

joint_prob=(temp_2+pseudo_count/(amino_N^2))/(pseudo_count+sum(Weight));

%% change the diagonal block
L=round(count/site_size);
for k=1:L
   temp=(k-1)*(site_size)+1:k*(site_size); % because one state is removed
   joint_prob(temp,temp)=zeros(site_size,site_size);
   for j=1:site_size
   joint_prob(temp(j),temp(j))=mean_sequence(temp(j));
   end
end


seq_corr=joint_prob-mean_sequence*mean_sequence';

%%   this is faster 
% mean_sequence=mean(binary_sequence,1);
% seq_corr=zeros(count,count);
% for k=1:M_new
% 
%          
%        seq_corr=seq_corr+(binary_sequence(k,:)-mean_sequence)'*(binary_sequence(k,:)-mean_sequence);
% 
% end
% 
% seq_corr=seq_corr./M_new;

%%

end


% max_N=count;
% figureParameter
% f1=image(1:max_N,1:max_N,1000000*seq_corr(1:max_N,1:max_N));
% colorbar;
% %set(gca, 'clim', [-0.25 0.25]);
% %set(gca, 'clim', [-1 1]);
% %a1=xlabel('$x$');
% %a2=ylabel('$y$');
% fig_name='./figure/mode_correlation.jpg';
% set(gca,'YDir','normal')
% print(fig_name,'-r400','-djpeg');
