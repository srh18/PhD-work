function [k1p,k1m,k2p,k2m,ka,Kw,K2,KH] = rate_constants(T,Ca)
if nargin ==1
    Ca = 1e-4;
    end
I = 3*Ca;
gOH = Debye(T,I,3,-1);
gH = Debye(T,I,9,1);
gHCO3 = Debye(T,I,4,-1);
gCO3 = Debye(T,I,4.5,-2);
gCa = Debye(T,I,6,2);
gamma = gH*gHCO3;
k1p = 10^(329.850-110.54*log10(T) - 17265.4/T);
k1m = 10^(13.558-3617.1/T);
k2p = 10^(13.635-2895/T);
k2m = 10^(14.09-5308/T);
K5 = 1.707e-4;
ka = k1m*gamma/K5;
Kw = 10^(22.801-4787.3/T -0.010365*T-7.1321*log10(T))/gH/gOH;
K2 = 10^(-107.8871+5151.79/T-0.03252849*T + 38.92561*log10(T) -563713.9/T^2)*gHCO3/gH/gCO3;
KH = 10^(108.3865 - 6919.53*T^(-1) + 0.01985076*T - 40.45154 *log10( T) + 669365*T^( -2));