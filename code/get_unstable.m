function A = get_unstable()
clf,hold on
del = [0 0.1 0.5 1 1.5 2];
L = (2:0.25:10)*pi;
j = length(del);
k = length(L);
A = zeros(j,k);

for i = 1:j
    for l = 1:k

        addlegend = 0;
        
            value = loadbigpR1(del(i),L(l));
            A(i,l) = value.most_unstable;
            

    end
end