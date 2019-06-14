function [flux,Diss,A]=Solve_periodic_Master_equation(birthRate1,birthRate2,deathRate1,deathRate2,transR1to2,transR2to1)
%%
%% also accept two inputs (for 1-d periodic )
%% birthRate1 is the forward rate, birthRate2 is the backward rate

%%
switch nargin
    case 6
        N=length(birthRate1);
        A=zeros(2*N,2*N); %operate on the vector: (P1(1),P1(2),P1(3),...,P1(N),P2(1),P2(2),...,P2(N))
        A(1,N)=birthRate1(N); A(1,2)=deathRate1(1); A(1,1)=-(birthRate1(1)+deathRate1(N)+transR1to2(1)); A(1,N+1)=transR2to1(1);

        for i=2:N-1
        A(i,i-1)=birthRate1(i-1); A(i,i+1)=deathRate1(i); A(i,i)=-(birthRate1(i)+deathRate1(i-1)+transR1to2(i)); A(i,i+N)=transR2to1(i);
        A(i+N,i+N-1)=birthRate2(i-1);A(i+N,i+N+1)=deathRate2(i); A(i+N,i+N)=-(birthRate2(i)+deathRate2(i-1)+transR2to1(i)); A(i+N,i)=transR1to2(i);
        end

        A(N,N-1)=birthRate1(N-1);A(N,1)=deathRate1(N); A(N,N)=-(birthRate1(N)+deathRate1(N-1)+transR1to2(N)); A(N,2*N)=transR2to1(N);
        A(N+1,2*N)=birthRate2(N); A(N+1,N+2)=deathRate2(1); A(N+1,N+1)=-(birthRate2(1)+deathRate2(N)+transR2to1(1)); A(N+1,1)=transR1to2(1);
        A(2*N,2*N-1)=birthRate2(N-1); A(2*N,N+1)=deathRate2(N); A(2*N,2*N)=-(birthRate2(N)+deathRate2(N-1)+transR2to1(N)); A(2*N,N)=transR1to2(N);

        [NormVector,orderEigValue]=orderedEigSystem(A,1);

        Prob=NormVector(:,1);

        % figure, plot(1:N,Prob(1:N),'-r',1:N,Prob(N+1:2*N),'-+b',1:N,Prob(1:N)+Prob(N+1:2*N),'*k')
        %  xlabel('position index');
        %  ylabel('Distribution');
        %  legend('state: 0','state: 1','state:1+0');



        % flux in the active state
        NetFlux_1=zeros(1,N);
        NetFlux_2=zeros(1,N); %flux in the inactive state

           NetFlux_1(1:N-1)=Prob(1:N-1).*birthRate1(1:N-1)'-Prob(2:N).*deathRate1(1:N-1)';
           NetFlux_1(N)=Prob(N).*birthRate1(N)-Prob(1)*deathRate1(N);

           NetFlux_2(1:N-1)=Prob(N+1:2*N-1).*birthRate2(1:N-1)'-Prob(N+2:2*N).*deathRate2(1:N-1)';
           NetFlux_2(N)=Prob(2*N)*birthRate2(N)-Prob(N+1)*deathRate2(N);

        totFlux=abs(NetFlux_1(1)+NetFlux_2(1));
        flux=zeros(1,2*N);
        flux(1:N)=NetFlux_1; flux(N+1:2*N)=NetFlux_2;


        Diss=sum(NetFlux_1.*log(birthRate1./deathRate1)+NetFlux_2.*log(birthRate2./deathRate2));

    case 2
     
        % birthRate1 is the forward rate, birthRate2 is the backward rate
        b=birthRate1;
        d=birthRate2;
        
        N=length(b);
        
        
        % transition rate matrix
            A=zeros(N,N);
            A(1,N)= b(N); A(1,1)=-(b(1)+d(N)); A(1,2)=d(1); 

        for i=2:N-1
            A(i,i-1)=b(i-1); A(i,i)=-(b(i)+d(i-1)); A(i,i+1)=d(i); 
        end

            A(N,N-1)=b(N-1); A(N,N)=-(b(N)+d(N-1));  A(N,1)=d(N);

        %% eigenvalues,  distribution,  dissipation
        [NormVector,orderEigValue]=orderedEigSystem(A,1);

        Prob=NormVector(:,1);

        flux=zeros(1,N);
        for i=1:N-1
           flux(i)=Prob(i)*b(i)-Prob(i+1)*d(i);
        end
        flux(N)=Prob(N)*b(N)-Prob(1)*d(N);

        Diss=sum(flux.*log(b./d));
        
        
        
        %v1=-sum(Prob.*log(b./d));
        %v2=mean(totFlux);


%         %figure, plot(X,-log(Prob./Prob(1)),'-+r',X,U(X)-U(X(1))-net_force*X,'*b')
%         figure, plot(X,U(X)-U(X(1))-net_force*X,'*b')
%         xlabel('x','FontSize',18);
%         xlim([0 1]);
%         legend('U(x)','FontSize',18);

        
        
        
        
    otherwise
        error('no input');
        
end
