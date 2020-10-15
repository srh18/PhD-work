function F = acalc(a,n,L, Bo, Re, ep,del)
z = 0:L/n:L-L/n;
F = 0*z;
F(n+1) = 0;
for m = 1:n
    [az,azz,azzz,azzzz] = getdiv(a,m,n,L);
    
    
    F(m) = a(m)^2*az + ep*(4/3*a(m)^3*az+3*a(m)^2*az*del*cos(2*pi/L*z(m))-2*pi/L*del*sin(2*pi/L*z(m))*a(m)^3+a(m)^2*az/Bo*(az+azzz+(8*pi^3/L^3-2*pi/L)*del*sin(2*pi/L*z(m)))+a(m)^3/(3*Bo)*(azz+azzzz+(16*pi^4/L^4 - 4*pi^2/L^2)*del*cos(2*pi/L*z(m)))-3/40*Re*(a(m)^6*azz+6*a(m)^5*az^2));
    F(n+1) = F(n+1) + a(m) + ep*(a(m)^2/2 + a(m)*del* cos(2*pi/L*z(m)));
end
F(n+1) = F(n+1) - (1+ep/2);
end

function [az,azz,azzz,azzzz] = getdiv(a,m,n,L)
if m ==1
    az = n/L*(-1/2*a(n)+1/2*a(m+1));
    azz =(n/L)^2*(a(n)-2*a(m)+ a(m+1));
    azzz= (n/L)^3*(-1/2*a(n-1)+a(n)-a(m+1)+1/2*a(m+2));
    azzzz =(n/L)^4 *(a(n-1)-4*a(n)+6*a(m)-4*a(m+1)+a(m+2));
elseif m==2
    az = n/L*(-1/2*a(m-1)+1/2*a(m+1));
    azz =(n/L)^2*(a(m-1)-2*a(m)+ a(m+1));
    azzz= (n/L)^3*(-1/2*a(n)+a(m-1)-a(m+1)+1/2*a(m+2));
    azzzz =(n/L)^4 *(a(n)-4*a(m-1)+6*a(m)-4*a(m+1)+a(m+2));
elseif m ==n
    az = n/L*(-1/2*a(m-1)+1/2*a(1));
    azz =(n/L)^2*(a(m-1)-2*a(m)+ a(1));
    azzz= (n/L)^3*(-1/2*a(m-2)+a(m-1)-a(1)+1/2*a(2));
    azzzz =(n/L)^4 *(a(m-2)-4*a(m-1)+6*a(m)-4*a(1)+a(2));
elseif m==n-1
    az = n/L*(-1/2*a(m-1)+1/2*a(m+1));
    azz =(n/L)^2*(a(m-1)-2*a(m)+ a(m+1));
    azzz= (n/L)^3*(-1/2*a(m-2)+a(m-1)-a(m+1)+1/2*a(1));
    azzzz =(n/L)^4 *(a(m-2)-4*a(m-1)+6*a(m)-4*a(m+1)+a(1));
else
    az = n/L*(-1/2*a(m-1)+1/2*a(m+1));
    azz =(n/L)^2*(a(m-1)-2*a(m)+ a(m+1));
    azzz= (n/L)^3*(-1/2*a(m-2)+a(m-1)-a(m+1)+1/2*a(m+2));
    azzzz =(n/L)^4 *(a(m-2)-4*a(m-1)+6*a(m)-4*a(m+1)+a(m+2));
end
end
