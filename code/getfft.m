function getfft(N,L)
clf
x = 0:L/N:L - L/N;
f = sin(x);
df = cos(x);
fhat = fft(f);
%xhat = fft(x)/N;
xhat = fftshift((1:N));

%xhat = fftshift(2*pi/L*(0:N-1));
dfhat = (pi*2i*xhat)/L.*fhat;
ddfhat = (pi*2i*xhat)/L.*dfhat;
dfFFT = real(ifft(dfhat));
ddfFFT = real(ifft(ddfhat));
plot(x,df)
hold on 
plot(x,dfFFT)
plot(x,ddfFFT)
end