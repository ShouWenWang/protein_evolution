function output_cell=correlation_response(alpha_0,beta_0,phi_0,lambda_0)
 % Note that lambda is defined to be  positive now.

%% re-ordering based on amplitude of alpha0, beta0 and phi0 
    N0=length(alpha_0);
    index=zeros(N0,2);
    
    index(:,1)=abs(alpha_0)+abs(beta_0)+abs(phi_0);
    index(:,2)=1:N0;


% index ordering
for i=2:N0
    temp=index(i,:);
    for j=i+1:N0
        if temp(1)<index(j,1)
           index(i,:)=index(j,:);
           index(j,:)=temp;
            temp=index(i,:);
         
        end
    end
end
        
alpha_1=alpha_0(index(:,2));
beta_1=beta_0(index(:,2));
phi_1=phi_0(index(:,2));
lambda_1=lambda_0(index(:,2));

% Use less mode to reduce the computation 
C0=sum(index(2:end,1)>0.001*index(2,1))+1; % threshold cutoff
Max_N=490;
if N0<Max_N 
    N1=N0;
% else if C0<Max_N
%         N1=Max_N;
%     else
%     N1=min([N0,C0]);
%     end
% end
else 
    N1=Max_N;
end
 

%% consider correlation function first


        C_time=@(t) 0;
        for j=2:N1
        C_time=@(t) C_time(t)+alpha_1(j).*beta_1(j).*exp(-lambda_1(j).*t);
        end

        C_fre=@(omega) 0;
        for j=2:N1
        C_fre=@(omega) C_fre(omega)+2.*alpha_1(j).*beta_1(j).*lambda_1(j)./(lambda_1(j)^2+omega.^2);
        end

        C_time_velo=@(t) 0;
        for j=2:N1
        C_time_velo=@(t) C_time_velo(t)-lambda_1(j).^2.*alpha_1(j).*beta_1(j).*exp(-lambda_1(j).*t);
        end

        C_fre_velo=@(omega) 0;
        for j=2:N1
        C_fre_velo=@(omega) C_fre_velo(omega)-2.*alpha_1(j).*beta_1(j).*lambda_1(j).*(1./(1+(omega./lambda_1(j)).^2)-1);
        end
    

%% consider response function

        R_time=@(t) 0;
        for j=2:N1
           R_time=@(t) R_time(t)+alpha_1(j).*phi_1(j).*exp(-lambda_1(j).*t);
        end
        
        
        
        R_fre=@(omega) 0;
        for j=2:N1
            R_fre=@(omega)  R_fre(omega)+alpha_1(j).*phi_1(j)./(lambda_1(j)-sqrt(-1)*omega);
        end
        
        R_time_velo=@(t) 0;
        for j=2:N1
            R_time_velo=@(t) R_time_velo(t)-alpha_1(j).*phi_1(j).*lambda_1(j).*exp(-lambda_1(j).*t);
        end
        
        R_fre_velo=@(omega) 0;
        for j=2:N1
            R_fre_velo=@(omega) R_fre_velo(omega)+alpha_1(j).*phi_1(j).*(1-1./(1-sqrt(-1).*(omega./lambda_1(j))));
        end
    
%% FRR violation for velocity observable
% Scheme I
% FRR_vio_velo_fre=@(omega) real(C_fre_velo(omega)-2*R_fre_velo(omega));
          
 
%Scheme II
%this naive scheme is wrong if the eigenvalue is complex ( this means
%alpha, beta, and phi are also complex)
        FRR_vio_velo_fre=@(omega) 0;
        for j=2:N1
           FRR_vio_velo_fre=@(omega) FRR_vio_velo_fre(omega)+2.*alpha_1(j).*(phi_1(j)-lambda_1(j).*beta_1(j))./(1+(omega./lambda_1(j)).^2);
            
        end
%           check consistency (there are consistent, although lambda is a complex number)
%            Omega=0.00001.*2.^(1:30);
%            figure, semilogx(Omega,FRR_vio_velo_fre_1(Omega),'-r',Omega,FRR_vio_velo_fre(Omega),'+k')
%            legend('analytical','direct');
%% susceptibility: integral of R_x(t)

        Chi_time=@(t) 0;
        for j=2:N1
           Chi_time=@(t) Chi_time(t)-alpha_1(j)*phi_1(j)/lambda_1(j)*(exp(-lambda_1(j).*t)-1);
        end
        
%% FRR violation integral
        %% relative contribution
        FRR_mode=zeros(1,N1-1);
        T_mode=zeros(1,N1-1);
        for j=2:N1
        FRR_mode(j-1)=lambda_1(j)*alpha_1(j)*(phi_1(j)-lambda_1(j)*beta_1(j));
        T_mode(j-1)=(lambda_1(j)*beta_1(j))/phi_1(j);
        end
        
%          index=(abs(phi_1)>10^(-5));
%         
%         figureParameter
%         semilogx(lambda_1(index),FRR_mode(index),'+r');
%         xlabel('Eigenvalue','FontSize',16);
%         ylabel('Dissipation','FontSize',16);
%         %fig_name='./figure/Adaptation-mode-a-FRR.eps';
%         fig_name='./figure/LaserSwitching-mode-FRR.eps';
%         figurePostTreat;
%         
%         figureParameter, semilogx(lambda_1(index),T_mode(index),'+r');
%         xlabel('Eigenvalue','FontSize',16);
%         ylabel('Temperature','FontSize',16);
%         %ylim([0 10]);
%         fig_name='./figure/LaserSwitching-mode-T.eps';
%         figurePostTreat;
%         %figure, plot(lambda_1(2:N1),FRR_mode,'-*r');
        
         FRR_integral=sum(FRR_mode);
%         FRR_integral=0;
%         for j=2:N1        %0.5 is due to transformation factor 1/2pi.  minus sign due to |\lambda|
%             FRR_integral=FRR_integral+0.5*2*lambda_1(j)*alpha_1(j)*(phi_1(j)-lambda_1(j)*beta_1(j)); % lambda is negative here.  So, different from that in the paper
%         end
        %FRR_integral_1=quadgk(FRR_vio_velo_fre,0,Inf)/pi;
        
         
        
        
%%        
        lambda_N=lambda_0(end);
        effective_gamma=1/R_fre_velo(100*abs(lambda_N));
        FRR_integral_diss=effective_gamma*FRR_integral;
        
        %% output
        output_cell(1,:)={R_time};
        output_cell(2,:)={R_fre};
        output_cell(3,:)={R_time_velo};
        output_cell(4,:)={R_fre_velo};
        output_cell(5,:)={C_time};
        output_cell(6,:)={C_fre};
        output_cell(7,:)={C_time_velo};
        output_cell(8,:)={C_fre_velo};
        output_cell(9,:)={FRR_vio_velo_fre};
        output_cell(10,:)={Chi_time};
        output_cell(11,:)={FRR_integral};
        output_cell(12,:)={FRR_integral_diss};
        output_cell(13,:)={effective_gamma};
        %%
        
end
