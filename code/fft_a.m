
acf = fft(val.ac,500,2);

for i  = 1:20:length(acf)
clf
plot(abs(acf(i,:)/500))
title(val.t(i))
xlim([0 10])
pause(0.01)
end