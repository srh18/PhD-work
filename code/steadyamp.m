test_value = shape_solveh;
test_value.n = 256;
test_value.ps = 1;
test_value.del = .1;
test_value.equation = 1;
R = [0.6 0.75 0.9 1 1.1 1.25 1.5 2 ];
L = pi/4:pi/20:6*pi;
figure(10),clf, hold on 
xlabel('$L$')
ylabel('$A$')
figure(11),clf, hold on 
xlabel('$L$')
ylabel('$\theta$')
figure(12),clf, hold on 
xlabel('$L$')
ylabel('$q$')
for r = R
    A = [];
    Theta = [];
    Q = [];
    test_value.R = r;
for l = L
    test_value.L = l;
    test_value = test_value.get_h;
    if test_value.eflag<=0
        amp = 0;
        theta = 0;
        q = 0;
    else
        [mh0,loc] = max(test_value.h0);
        amp =(mh0-min(test_value.h0))/2;
    end
    A = [A , amp];
    Theta = [Theta,loc/test_value.n];
    Q = [Q, test_value.q];
end
figure(10)
plot(L,A)
figure(11)
plot(L,Theta)
figure(12)
plot(L,Q)
end
