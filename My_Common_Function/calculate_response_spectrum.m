function [Rx_fre,Rv_fre]=calculate_response_spectrum(x,h,Hz,dt,plot_on)

L=length(x);
X_fre=(fft(x,L)/L);
h_fre=fft(h,L)/L;
fre=2*pi*(0:L/2)'/(L*dt); 
Omega=2*pi*Hz;  %perturbed frequency

index=round(L*dt*Hz+1); % index of given perturbed frequency
Y=X_fre(index);     
Z=h_fre(index);
Rx_fre=Y'./Z;
Rv_fre=Omega.*Rx_fre;

if plot_on
    figure,plot(fre,real(X_fre(1:length(fre))),'-r',fre,imag(X_fre(1:length(fre))),'-.b');
    xlim([1,10]);
    legend('real','imag');
end
