function finite_diff_stal(F,delta,B, a, b,dt,nz,T,L)
dz = L/nz;
z = 0:dz:L;

eta = sin(2*pi*z/L);
etan =0*eta;

for t= 0:dt:T
    etan(1) = eta(1) + dt*F*(1+delta^2*(-a*z(1)+b+1/3*(((-3/2*eta(1)+2*eta(2)-1/2*eta(3))/dz).^2 - 2*eta(1) - 2*B/dz^3*(-5/2*eta(1)+9*eta(2)-12*eta(3)+7*eta(4)-3/2*eta(5))))); 
    etan(2) = eta(2) + dt*F*(1+delta^2*(-a*z(2)+b+1/3*(((eta(3)-eta(1))/(2*dz)).^2 - 2*eta(2) - 2*B/dz^3*(-5/2*eta(2)+9*eta(3)-12*eta(4)+7*eta(5)-3/2*eta(6))))); 
    etan(nz) =  eta(nz) + dt*F*(1+delta^2*(-a*z(nz)+b)); 
    etan(nz+1) =  eta(nz+1) + dt*F*(1+delta^2*(-a*z(nz+1)+b));
%     for j = [1 2 nz nz+1]
%         etan(j) = eta(j) + dt*F*(1+delta^2*(-a*z(j) +b-2/3*eta(j)));
%     end
    for i =3:nz-1
        etan(i) = eta(i) + dt*F*(1+delta^2*(-a*z(i)+b+1/3*(((eta(i+1)-eta(i-1))/(2*dz)).^2 - 2*eta(i) - 2*B*((1/2*eta(i+2) -eta(i+1)+eta(i-1)-1/2*eta(i-2))/dz^3))));
        dddeta(i)  = (1/2*eta(i+2) -eta(i+1)+eta(i-1)-1/2*eta(i-2))/dz^3;
    end
    eta = etan;
    plot(eta,z)

    pause(1)
    
    
    
end
plot(eta(1:nz-1),z(1:nz-1))
end