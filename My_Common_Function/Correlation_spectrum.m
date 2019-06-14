  
function [Cx_fre,Cv_fre,fre]=Correlation_spectrum(x,dt,cutoff_time,plot_on)
%%
% obtain correlation spectrum by truncating the long time noise for
% Cx(time)

% cutoff_time,  choose to be around 10*tau_s (tau_s:  the slow timescale of the system)

%%


L2=length(x);

%velo=x(end)/(L2*dt);
%x=x-velo*(0:L2-1)'*dt;
Cx_time_0=ifft(abs(fft(x,L2)).^2)/L2; % zero padding is irrelevant since the data is essentially uncorrelated for L

Cx_time=Cx_time_0(1:L2/2+1);
time=dt*(1:L2/2+1);
% figureParameter
% f1=plot(time(2:end),Cx_time(2:end),'-r');
% xlim([0.001 1000]);
% a1=xlabel('Time');
% a2=ylabel('$\tilde{C}_v$');
% fig_name='./figure/response_velo_time.eps';
% figurePostTreat

%% Frequency spectrum from Cx_time: truncation of long time  noise
% this truncation introduce slow modulation power spectrum in high
% frequency domain. 


if round(cutoff_time/dt)<length(Cx_time)
    L3=round(cutoff_time/dt);
else
    L3=length(Cx_time);
end
newL=2^nextpow2(L3);


newTime=dt*(1:2*newL);
newCx=zeros(2*newL,1);

cutoff=1/dt;
newCx(1:newL)=Cx_time_0(1:newL);
%alpha=1; %obtained by fitting to log Cx(t)
%newCx(cutoff+1:newL)=Cx_time_0(cutoff+1).*exp(-alpha*dt*(0:newL-cutoff-1));

%newCx(1:newL)=Cx_time_0(1:newL);
newCx(2*newL:-1:newL+1)=newCx(1:newL);



%% final result
fre=2*pi*(0:newL-1)'/(2*newL*dt);
Cx_fre=abs(fft(newCx))*dt;
Cx_fre=Cx_fre(1:newL);


Cv_fre=fre.^2.*Cx_fre;

if plot_on
    figureParameter
    plot(fre,Cx_fre);
    xlim([0.01 1]);
    fig_name='./figure/spectrum.eps';
    figurePostTreat;

end
