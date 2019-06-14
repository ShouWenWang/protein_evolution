function predicting_inverse_matrix_entries(mean_sequence,mean_sequence_old,remove_q,inv_corr,amino_N,Delta)
% mean_sequence_old, the mean of each residue before removing the
% referenced state
%% the effective Delta
%Delta=Delta_new(index);
new_residue_N=amino_N-1;
L=round(length(mean_sequence)/new_residue_N);
temp_1=(0:L-1)*amino_N+remove_q;
mean_removed=zeros(1,L);
fitting_parameter=0.15; 
row_to_check=50;

for j=1:L
   mean_removed(j)=mean_sequence_old(temp_1(j)); 
end


%%%%% compute theoretical predictions
[count,count]=size(inv_corr);
inv_diag=zeros(1,count);
predict_inv_diag=zeros(1,count);
for j=1:L
    for k=1:new_residue_N
     temp_index=new_residue_N*(j-1)+k;
     inv_diag(1,temp_index)=inv_corr(temp_index,temp_index);
     predict_inv_diag(temp_index)=1/mean_removed(j)+1/mean_sequence(temp_index)+fitting_parameter*abs(Delta(temp_index))^2;
    end
end

predict_inv=zeros(count,count);
conservation_diag_block=zeros(count,count); % the effect of conservation in diagonal block
for j=1:L
    for j1=1:L
      for k=1:new_residue_N
          for k1=1:new_residue_N
               temp_index1=new_residue_N*(j-1)+k;
               temp_index2=new_residue_N*(j1-1)+k1;
               %%%%%%%%%%%%%%%
               if j==j1 && k==k1
               predict_inv(temp_index1,temp_index2)=fitting_parameter*abs(Delta(temp_index1))^2+1/mean_removed(j)+1/mean_sequence(temp_index1);
               conservation_diag_block(temp_index1,temp_index2)=1/mean_removed(j)+1/mean_sequence(temp_index1);
               else
                   if j==j1 && k~=k1
                     predict_inv(temp_index1,temp_index2)=fitting_parameter*Delta(temp_index1)*Delta(temp_index2)+1/mean_removed(j);
                     conservation_diag_block(temp_index1,temp_index2)=1/mean_removed(j);
            


                   else
                     predict_inv(temp_index1,temp_index2)=fitting_parameter*Delta(temp_index1)*Delta(temp_index2);


                   end
               end
             %%%%%%%%%%%%%%
          end
      end
    end
end



figureParameter
f1=plot(1:count,inv_corr(row_to_check,:),'-r',1:count,predict_inv(row_to_check,:),':k')
xlim([0 count+1]);
xlabel('Site index');
ylabel('Entry value');
title("Row of "+num2str(row_to_check));
h=legend('Data','Prediction');
%a2=ylabel('Distribution');
fig_name='./figure/predicted_inv_1.eps';
figurePostTreat




figureParameter
f1=plot(1:count,inv_diag,'-r',1:count,predict_inv_diag,':k')
xlim([0 count+1]);
xlabel('Site index');
ylabel('Entry value');
title('Diagonal terms');
h=legend('Data','Prediction');
%a2=ylabel('Distribution');
fig_name='./figure/predicted_inv_diagonal.eps';
figurePostTreat


% figureParameter
% f1=image(1:count,1:count,2*abs(inv_corr_off_diag));
% colorbar;
% %set(gca, 'clim', [-0.25 0.25]);
% %set(gca, 'clim', [-1 1]);
% %a1=xlabel('$x$');
% %a2=ylabel('$y$');
% fig_name='./figure/corr1.jpg';
% set(gca,'YDir','normal')
% print(fig_name,'-r400','-djpeg');
