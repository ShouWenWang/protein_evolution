function [seq_corr,mean_sequence]=selection_covariance(sequence_array,sum_K,d,gamma,method)
%%  The select the sequences, then compute the covariance matrix

% sum_K, the summation of mutation effects (i.e., the total energetic cost for our PDZ model) for a given sequence
% d, selection range;  
% gamma, the selection bias
% method, two-entry variable [a,b], how to perform selection: 
% .  [a=1,b=1], reweight each sequence using the quadratic free energy
% .  [a=1,b=0], reweight each sequence using the quartic free energy
%    [a=0,b=1], Window-selection to get some sequences, then each sequence has the same weight
%    [a=0,b=0], Threshold-selection to get some sequences, then each one has the same weight      
%
%



%% sequence_array(M,count)
if method(1)  % reweight
    
     if method(2) % quartic selection    
        kappa=1/d^2;
        weighting_fun=@(x) sqrt(kappa/(2*pi))*exp(-0.5*kappa*(x-gamma).^2);

     else %% quartic
        kappa=1/d^4;
        weighting_fun=@(x) power(4*kappa,0.25)/3.62*exp(-0.25*kappa*(x-gamma).^4);
     end


    %%

    Weight=weighting_fun(sum_K);
    % covariance matrix 
    pseudo_ratio=0;
    [seq_corr,mean_sequence]=covariance_matrix(sequence_array,Weight,pseudo_ratio);
    
else % not reweight, but select
    
    if method(2) 
        index_center=sum_K>-d+gamma & sum_K<d+gamma; % window selection
    else
        index_center=sum_K>gamma; % threshold selection
    end

    M_new=sum(index_center);
    sequence_array_new=sequence_array(index_center,:);


    %%
    mean_sequence=mean(sequence_array_new,1);
    count=length(mean_sequence);
    %% covariance analysis
    seq_corr=zeros(count,count);
    for k=1:M_new

           seq_corr=seq_corr+(sequence_array_new(k,:)-mean_sequence)'*(sequence_array_new(k,:)-mean_sequence);

    end

    seq_corr=seq_corr./M_new;

    
end