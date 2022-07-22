test_value = shape_solveh;
test_value.n = 256;
test_value.ps = 1;
test_value.R = 1;
test_value.equation = 1;
del = [0.1 0.5 1 1.5 2 ];
L = pi/4:pi/20:6*pi;
theta = atan(3*L.^3./(2*pi*(4*pi^2-L.^2))) +pi;
figure(15),clf, hold on 
xlabel('$L$')
ylabel('$A$')
plot(L,-cos(theta))

figure(16),clf, hold on 
xlabel('$L$')
ylabel('$\theta$')
plot(L,theta/(2*pi))

for d = del
    A = [];
    Theta = [];
    
    test_value.del  = d;
for l = L
    test_value.L = l;
    test_value = test_value.get_h;
    if test_value.eflag<=0
        amp = 0;
        theta = 0;
        
    else
        [mh0,loc] = max(test_value.h0);
        amp =(mh0-min(test_value.h0))/2/test_value.del;
    end
    A = [A , amp];
    Theta = [Theta,loc/test_value.n];
    
end
figure(15)
plot(L,A)
figure(16)
plot(L,Theta)

end
