function [NormVector,orderEigValue]=orderedEigSystem(A,Is_Markov_Process)
%% Purpose
% input: evolution Matrix A for a Markov Process
% Output: eigValue: all set of eigenvalues for A in a descending order
%         eigVector: a complete set of eigenvectors corresponding to the order of eigValue   

%%  
[vector,eigvalue]=eig(A,'balance');
D=size(A);

orderEigValue=zeros(D(1),2); % the first input is eigenvalue, 
%the second is the column index in the original problem

for i=1:D(1)
    orderEigValue(i,1)=eigvalue(i,i);
    orderEigValue(i,2)=i;
end


%Is_Markov_Process=0;

if Is_Markov_Process==1
%% eigenvalue ordering in a descending way, for Markov process
for i=1:D(1)
    temp=orderEigValue(i,:);
    for j=i+1:D(1)
        if temp(1)<orderEigValue(j,1)
           orderEigValue(i,:)=orderEigValue(j,:);
           orderEigValue(j,:)=temp;
            temp=orderEigValue(i,:);
         
        end
    end
end
        
% Eigenvector ordering and Normalization
NormVector=zeros(D(1),D(1));

%Markov process

for i=1:D(1)
    norm=sqrt(dot(vector(:,orderEigValue(i,2)),vector(:,orderEigValue(i,2))));
    NormVector(:,i)=vector(:,orderEigValue(i,2))/norm;
    
     %make sure that the majority is positive
    sign_flag=sign(sum(NormVector(:,i)));
    NormVector(:,i)=sign_flag*NormVector(:,i);
end

%probability normalizaton
   NormVector(:,1)=vector(:,orderEigValue(1,2))/sum(vector(:,orderEigValue(1,2)));
   
   orderEigValue(:,2)=(1:D(1))';
   
   
   %% order the eigenvalue in the descending way
else 

        % eigenvalue ordering
    for i=1:D(1)
        temp=orderEigValue(i,:);
        for j=i+1:D(1)
            if temp(1)<orderEigValue(j,1)
               orderEigValue(i,:)=orderEigValue(j,:);
               orderEigValue(j,:)=temp;
                temp=orderEigValue(i,:);

            end
        end
    end

    % Eigenvector ordering and Normalization
    NormVector=zeros(D(1),D(1));

    %Markov process

    for i=1:D(1)
        norm=sqrt(dot(vector(:,orderEigValue(i,2)),vector(:,orderEigValue(i,2))));
        NormVector(:,i)=vector(:,orderEigValue(i,2))/norm;
        %make sure that the majority is positive
        sign_flag=sign(sum(NormVector(:,i)));
        NormVector(:,i)=sign_flag*NormVector(:,i);
    end

       orderEigValue(:,2)=(1:D(1))';

end



