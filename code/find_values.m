function [vals]=find_values()
syms k S
syms A B C D s
syms x r gam a Bo
p = A*I_0(x)+B*K_0(x);
u = -(1i*C+A/k)*I_1(x)+(1i*D+B/k)*K_1(x)+A*x/(2*k)*I_0(x)+B*x/(2*k)*K_0(x);
w = C*I_0(x)+D*K_0(x)+1i*A*x/(2*k)*I_1(x)-1i*B*x/(2*k)*K_1(x);

tang = 1i*k*u + k*diff(w);

eqn1 = subs(u,x,k) == 0;
eqn2 = subs(w,x,k) == 0;
eqn3 = subs(k*diff(u),x,k) == 0;
eqn4 = subs(tang,x,k) == 0;
eqn5 = p -a == (1-k^2)*s;

[M,y] = equationsToMatrix([eqn1, eqn2, eqn3, eqn4, eqn5], [A, B, C, D, s]);
vals = simplify(linsolve(M,y));
end