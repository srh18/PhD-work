classdef shape_solveh
    
    %Class to help solve Thin film down a cylinder of similar size to the
    %wavelength
    
    properties
        n {mustBeInteger,mustBePositive} = 256  %number of grid points INT
        
        % Fluid properties Parameters
        
        Bo {mustBeNonnegative} = 1 % Bond number
        Re {mustBeNonnegative} = 1 % Reynolds number
        ep {mustBeNonnegative} = 1e-1 % Ratio between the fluid thickness and the radius of the cylinder
        R {mustBeNonnegative} = 1 % Radius of cylinder
        
        q = 1/3 % flow rate
        T {mustBeNonnegative} = 100 %time
        delt  {mustBePositive} = 0.05 % Time step
        mass0 = 1 %initial mass
        
        %Wall shape parameters
        
        del {mustBeNonnegative} = 1 % Disturbance amplitude
        L {mustBeNonnegative} = 2*pi % wavelength of disturbance
        l {mustBeNonnegative,mustBeLessThanOrEqual(l,1)} = 0.5 %width of the step
        del2 {mustBeNonnegative} = 0.1 % steepness of the step
        dir {mustBeMember(dir,[-1,1])} = 1 %direction of slope
        ps {mustBeMember(ps,[0,1])} = 0; %whether to use Psuedo spectral differentiation rather than finite differences
        % Asthetic features
        
        periods {mustBePositive,mustBeInteger} = 1 % number of periods to plot INT - ignored for now
        
        %Derived from parameters
        
        z % domain
        t % time domain
        h % mass conserving fluid thickness
        
        eta % wall shape
        h0 %steady state
        
        
        
        % Info about code
        
        integration_time
        init_condition 
        
        %Toggle
        equation = 0 %which equation to use 0 for Benny like, 1 for surface tension dominating 10 for not in conservative form.
        paramtoggle  {mustBeMember(paramtoggle,[1,2,3,4,5,6,7])} = 1 % 1 for Bo, 2 for R, 3 for L, 4 for Re, 5 for ep, 6 for del, 7 for l
        flow {mustBeMember(flow,[0,1])} = 0 % 0 if keeping volume constant - 1 if flow rate is constant
        fdo {mustBeMember(fdo,[2,4])}= 2 % Order of the finite difference scheme used (may add later)
        wall_shape {mustBeMember(wall_shape,[0,1,2,3,4])} = 0 % 0 for cosine, 1 step, 2 sawtooth, 3 for flat, 4 user input
        boundary {mustBeMember(boundary,[0,1])} = 0 % 0 for periodic boundary conditions, 1 for not
        force_mass {mustBeMember(force_mass,[0,1])}= 0 %force mass conservation on integration
        usejac = 0 %whether to use a jacobian
        suppression = 1e-12 %remove noise from fft
        pad = 4
    end
    
    properties(Hidden)
        %additional outputs from fsolve
        Fval
        eflag
        out
        jac
        
        c %peak speed
        ac %wall moving at peak speed
        
        pks
        loc
        W
        P
        goodP
        

        
        zpos
        speed
        
        %string of parameters
        params = {'Bo','R','L','Re','\epsilon','\delta'}
        % Properties of the fluid
        
        sol %ode output
        reltol = 1e-5 %relative tolerance for ode15s
        abstol = 1e-8 %absolute tolerance for ode15s
        hmax % maximum fluid thickness
        hmin % minimum fluid thickness
        diff % range of fluid thickness
        hmaxloc % location of the maximum
        hminloc % location of the minimum
        maxmindiff % distance between the maximum and the minimum
        minmaxdiff % distance between the minimum and the maximum
        notol = 0 %whether to use default tolerances
        npks %number of peaks
        
        peakloc {mustBeMember(peakloc,[0,1,2,3])} = 0
        peakpos
        peakspeed
        
        class
        notes
        %derivatives of eta to speed things up
        etaz
        etazz
        etazzz
        etazzzz
        
        %save time with the Jacobian if these are already calculated
        hz
        hzz
        hzzz
        hzzzz
        
        integration {mustBeMember(integration,[0,1,2])} = 0 % 0 is crank Nicholson 1 implicit 2 explicit
        
        
    end
    properties(Dependent = true,Hidden)
        S %position from origin
        nz %normalised z
         % effective thickness
        mass %integral of h
        h2norm % integral of h^2
        amass % integral of a
        amassgrowth
        a2norm %integral of a^2
        Q % flow rate
        Qint % integral of flow rate
        hdiff %differenence from the steady state
        at %time derivative of amass
        a
        
    end
    
    
    methods
        function obj = shape_solveh(n,L,Bo,Re, ep, del, periods)
            
            if nargin>0
                obj.n = n ;
                obj.L = L;
                
                obj.Bo = Bo;
                obj.Re = Re;
                obj.ep = ep;
                obj.del = del;
                obj.periods = periods;
                
            end
            
            obj.delt = 1/obj.n;
            obj = obj.reset;
            
            
        end
        %dependent function
        
        function value = get.a(obj)
            if obj.equation == 0
                value = (sqrt(2*obj.h*obj.R/obj.ep +  (obj.eta+obj.R/obj.ep).^2)-(obj.R/obj.ep+ obj.eta));
            else
                value = obj.h;
            end
        end
        function value = get.S(obj)
            value = obj.a+ obj.eta;
        end
        function value = get.nz(obj)
            value = obj.z/obj.L;
        end
        function value = get.at(obj)
            value = diff(obj.amass)./diff(value.t);
        end
        function value = get.mass(obj)
            if obj.equation == 10
                value = sum(obj.h+ obj.ep/obj.R*(obj.h.*obj.eta+obj.h.^2/2),2)/obj.n;
            else
                value = sum(obj.h,2)/obj.n;
            end
        end
        function value = get.Q(obj)
            [hz,~,hzzz,~] = obj.getdiv(obj.h);
            if obj.equation == 0
                value = obj.h.^3/3+ obj.ep*(2/15*obj.Re*obj.h.^6.*hz + obj.h.^3/3/obj.Bo.*((hz+obj.etaz)/obj.R^2+ hzzz+ obj.etazzz)-2/3/obj.R*obj.h.^3.*obj.eta - 1/6/obj.R*obj.h.^4);
            elseif obj.equation == 1
                value  = obj.h.^3/3 +  obj.ep*(obj.h.^3/3/obj.Bo.*((hz+obj.etaz)/obj.R^2+ hzzz+ obj.etazzz));
            end
            
        end
        function value = get.Qint(obj)
            value = sum(obj.Q,2)/obj.n;
        end
        function value = get.amass(obj)
            value = sum(obj.a,2)/obj.n;
        end
        function value = get.amassgrowth(obj)
            
            value = sum(-obj.a.^3.*obj.etaz*obj.ep/obj.R/3,2)/obj.n;
            
        end
        function value = get.h2norm(obj)
            value = sum(obj.h.^2,2)/obj.n;
        end
        
        function value = get.hdiff(obj)
            value = obj.h - obj.h0;
        end
        %Functions which also change other parameters when changed
        function obj = set.L(obj,value)
            obj.L = value;
            obj = obj.reset;
        end
        
        function obj = set.equation(obj,value)
            obj.equation = value;
            if value == 2
                obj.wall_shape = 3;
            end
        end
        
   
        function obj = set.n(obj,value)
            obj.n = value;
            obj = obj.reset;
            
        end
        
        function obj = set.del(obj,value)
            obj.del = value;
            obj = obj.get_wall;
        end
        
        function obj = set.del2(obj,value)
            obj.del2 = value;
            obj = obj.get_wall;
        end
        
        function obj = set.eta(obj,value)
            obj.eta = value;
            [obj.etaz,obj.etazz,obj.etazzz,obj.etazzzz] = obj.getdiv(value);
            if length(obj.z)~=length(obj.eta)
                error('Vectors for the wall shape and the domain must be the same length.')
            end
            
        end
        
        function obj = set.l(obj,value)
            if value >1 || value<0
                warning('l must be between 0 and 1 and so has not been changed')
            else
                obj.l = value;
                obj = obj.reset;
            end
        end
        function obj = set.periods(obj,value)
            obj.periods = value;
            obj = obj.reset;
        end
        
        function obj = set.wall_shape(obj,value)
            obj.wall_shape = value;
            if obj.wall_shape <4
                obj = obj.get_wall;
            end
        end
        
        %         function obj = set.T(obj,value)
        %             obj.T = value;
        %             obj.t = 0:obj.delt:obj.T;
        %         end
        %         function obj = set.delt(obj,value)
        %             obj.delt = value;
        %             obj.t = 0:obj.delt:obj.T;
        %         end
        
        function obj = reset(obj)
            %realign z, wall if n or L are changed
            obj.z = 0: obj.L/obj.n: obj.periods*obj.L- obj.L/obj.n;
            obj = obj.get_wall;
            obj.q = 1/3 - 1/12*obj.ep/obj.R;
            
            
        end
        
        function obj = makeascale(obj)
            obj.mass0 = 1+obj.ep/(obj.R*2);
        end
        
        function obj = setp(obj,value)
            if obj.paramtoggle == 1
                obj.Bo = value;
            elseif obj.paramtoggle == 2
                obj.R = value;
            elseif obj.paramtoggle == 3
                obj.L = value;
            elseif obj.paramtoggle == 4
                obj.Re = value;
            elseif obj.paramtoggle == 5
                obj.ep = value;
            elseif obj.paramtoggle == 6
                obj.del = value;
            elseif obj.paramtoggle ==7
                obj.l = value;
            end
        end
        
        function value = getp(obj)
            if obj.paramtoggle == 1
                value = obj.Bo ;
            elseif obj.paramtoggle == 2
                value = obj.R ;
            elseif obj.paramtoggle == 3
                value = obj.L ;
            elseif obj.paramtoggle == 4
                value = obj.Re ;
            elseif obj.paramtoggle == 5
                value = obj.ep;
            elseif obj.paramtoggle == 6
                value = obj.del;
            elseif obj.paramtoggle ==7
                value = obj.l;
            end
        end
        
        
        
        function obj = get_wall(obj)
            %Function to create the wall depending on the wall shape
            if obj.wall_shape == 1
                obj.eta = obj.del/2*(tanh((obj.z/obj.L-(1-obj.l)/2)/obj.del2) - tanh((obj.z/obj.L - (1+obj.l)/2)/obj.del2));
                
            elseif  obj.wall_shape ==2
                obj.eta = (obj.z/obj.L-(1-obj.l)/2)/(obj.l).*(obj.del/2*(tanh((obj.z/obj.L-(1-obj.l)/2)/obj.del2) - tanh((obj.z-obj.L*(1+obj.l)/2)/obj.del2)));
            elseif obj.wall_shape == 0
                
                obj.eta = obj.del*cos(2*pi/obj.L*obj.z);
            elseif obj.wall_shape == 3
                obj.eta = 0*obj.z;
            else
                warning('Wall Shape may no longer match coordinates')
            end
            
        end
        
        function obj = get_ac(obj,c)
            if nargin == 2
                obj.c = c;
            end
            
            phi = mod((obj.z-obj.c*obj.t),obj.L);
            [~,I] = sort(phi,2);
            I = I + obj.n*(0:length(I)-1)';
            ac = obj.a';
            obj.ac = ac(I);
        end
        
        
        %MAIN FUNCTIONS
        
        function hinit = linear_shape(obj)
            % warning uses old equation may be incorrect.
            %analytic solution to first order of a disturbance
            %del*cos(2*pi/L) (where del<<1) used as a initial guess for
            %fsolve
            a0 = 1;
            k = 2*pi/obj.L;
            A = -k  +obj.ep*(2*k./(3.*obj.R));
            B = obj.ep*(-2/15*obj.Re*k.^2+1./(3.*obj.Bo)*(-k.^2/obj.R.^2+k.^4));
            C = obj.ep*(2.*k./(3*obj.R));
            D = (obj.ep./(3.*obj.Bo).*(-k.^2./obj.R.^2+k.^4));
            a = -(A.*C+B.*D)./(A.^2+B.^2);
            b = (B.*C-D.*A)/(A.^2+B.^2);
            
            hinit = 1+ a.*obj.del.*cos(2*pi./obj.L.*obj.z)-b.*obj.del.*sin(2*pi./obj.L.*obj.z);
            
        end
        function [A,lambda] = shift(obj,values, param,linear)
            if nargin>2
            obj.paramtoggle = param;
            end
            if nargin<4
                linear = 0;
            end
            if linear ==1 
                obj.del = 1;
            end
           
            for i = 1:length(values)
                obj = obj.setp(values(i));
                if linear ==1 
                    obj.h = obj.linear_shape;
                else
                obj = obj.get_h;
                obj.h = obj.h0;
                end
                [hmax,hi] = max(obj.h);
                A(i) = (hmax - min(obj.h))/2/obj.del;
                lambda(i) = obj.z(hi);
            end
                
           
        end
                function [A,lambda] = ashift(obj,values, param,linear)
            if nargin>2
            obj.paramtoggle = param;
            end
            if nargin<4
                linear = 0;
            end
            if linear ==1 
                obj.del = 1;
            end
           
            for i = 1:length(values)
                obj = obj.setp(values(i));
                if linear ==1 
                    obj.h = obj.linear_shape;
                else
                obj = obj.get_h;
                obj.h = obj.h0;
                end
                [hmax,hi] = max(obj.a);
                A(i) = (hmax - min(obj.a))/2/obj.del;
                lambda(i) = obj.z(hi)/obj.L;
            end
                
           
        end
        function [A,lambda] = linear_shift(obj) 
             a0 = 1;
            k = 2*pi./obj.L;
            A = -(k  -obj.ep*(4.*k./(6.*obj.R)));
            B = obj.ep*(-2/15.*obj.Re.*k.^2+1./(3.*obj.Bo)*(-k.^2./obj.R.^2+k.^4));
            C = obj.ep*(2./(3.*obj.R));
            D = -(obj.ep./(3.*obj.Bo).*(-k.^2./obj.R.^2+k.^4));
            a = -(A.*C+B.*D)./(A.^2+B.^2);
            b = -(B.*C-D.*A)./(A.^2+B.^2);
            
            A = sqrt(a.^2 +b.^2);
            lambda = 1./k*atan(-b./a);
        end
        function ainit = linear_shapea(obj)
            % warning uses old equation may be incorrect.
            %analytic solution to first order of a disturbance
            %del*cos(2*pi/L) (where del<<1) used as a initial guess for
            %fsolve
            
            
            Rep = obj.R/obj.ep;
            k = 2*pi/obj.L;
            
            a0 = -Rep +sqrt(Rep^2+2*Rep*obj.mass0);
            
            
            A = -(a0^2*k  +obj.ep*(a0^3*k/(3*obj.R)));
            B = obj.ep*(-2/15*obj.Re*a0^6*k^2+a0^3/(3*obj.Bo)*(-k^2/obj.R^2+k^4));
            C = -obj.ep*(a0^3/(3*obj.R)*k);
            D = (obj.ep*a0^3/(3*obj.Bo)*(-k^2/obj.R^2+k^4));
            a = -(A*C+B*D)/(A^2+B^2);
            b = (D*A-B*C)/(A^2+B^2);
            
            ainit = a0+ b*obj.del*sin(2*pi/obj.L*obj.z)+a*obj.del*cos(2*pi/obj.L*obj.z);
        end
        
        function a = small_ep_shape(obj)
            %warning uses old equation may be incorrect
            a = 1 - obj.ep/3*obj.del*(cos(2*pi/obj.L*obj.z)+1/obj.Bo*(-2*pi/obj.L+8*pi^3/obj.L^3)*sin(2*pi/obj.L*obj.z));
        end
        
        function F = hfun(obj,h)
            
            if obj.flow == 1
                [hz,~,hzzz,~] = obj.getdiv(h);
                F = h.^3/3 + obj.ep*(2/15*obj.Re*h.^6.*hz+h.^3/(3*obj.Bo).*((hz+obj.etaz)/obj.R^2+ hzzz+ obj.etazzz)-h.^4/(6*obj.R)-2/(3*obj.R)*h.^3.*obj.eta) - obj.q;
            else
                l = length(h);
                [hz,~,hzzz,~] = obj.getdiv(h(1:l-1));
                H = h(1:l-1);
                h(l);
                if obj.equation == 0
                    
                    F = H.^3/3 + obj.ep*(2/15*obj.Re*H.^6.*hz+H.^3/(3*obj.Bo).*((hz+obj.etaz)/obj.R^2+ hzzz+ obj.etazzz)-H.^4/(6*obj.R)-2/(3*obj.R)*H.^3.*obj.eta) - h(l);
                elseif obj.equation ==1
                    F = H.^3/3+H.^3/(3*obj.Bo).*((hz+obj.etaz)/obj.R^2+ hzzz+ obj.etazzz) - h(l);
                end
                F(l) = sum(H)/obj.n - obj.mass0;
                %F = hz.*h.^2 + obj.ep.*(-4*hz.*h.^3/(12.*obj.R)+2/15.*obj.Re.*(h.^6.*hzz+6.*h.^5.*hz.^2)+1/(3.*obj.Bo).*(h.^3.*((hzz+obj.etazz)/obj.R^2+hzzzz+obj.etazzzz)+3*hz.*h.^2.*((hz+obj.etaz)/obj.R^2+hzzz+obj.etazzz))-2*obj.etaz.*h.^3/(3*obj.R)-2*hz.*h.^2.*obj.eta/obj.R);
            end
        end
        %         function F = hfunflux(obj,h,q)
        %             [hz,hzz,hzzz,hzzzz] = obj.getdiv(h);
        %
        %             F = h.^3/3 + obj.ep*(2/15*obj.Re*h.^6.*hz+h.^3/(3*obj.Bo).*((hz+obj.etaz)/obj.R^2+ hzzz+ obj.etazzz)-h.^4/(12*obj.R)-2/(3*obj.R)*h.^3.*obj.eta) - q;
        %         end
        %
        %
        
            function [az,azz,azzz,azzzz] = getdiv(obj,a)
                dim = size(a);
                if obj.ps ==1
                    %DIFF_PSEUDO_SPECTRAL Uses the pseudo-spectral method to differentiate
                    %   Detailed explanation goes here
                    
                    %suppression = 1e-12;
                    
                    
                    % Transform into fourier space
                    yF = fft(a);
                    
                    N = obj.n/2;
                    
                    % Determine k in matlab form
                    %k = fftshift([0,-N+1:N-1])';
                    
                    
                    % Prior suppression
                    yF(abs(yF)<obj.suppression) = 0;
                    %yF = [yF,zeros(1,3*obj.n)]*4;
                    % Apply pseudo-spectral differentiation
                    padzeros = zeros(1,(obj.pad-1)*obj.n);
                    
                    yF = [yF(1:N),padzeros, yF(N+1:obj.n)]*obj.pad;
                    k = [0:N-1,padzeros, 0, 1-N:-1] * 2*pi/obj.L;
                    % Posterior suppression
                    % dyF(abs(dyF) < suppression*N*2) = 0 ;
                    % dyF(abs(dyF) < suppression*max(abs(dyF))) = 0 ;
                    
                    % Transform back into real space
                    dyF = (1i*k).^1.*yF;
                    az = real(ifft(dyF));
                    %az = real(ifft(dyF,obj.n*obj.pad));
                    az = az(1:obj.pad:end);
                    dyF = (1i*k).^2.*yF;
                    azz = real(ifft(dyF));
                    %azz = real(ifft(dyF,obj.n*obj.pad));
                    azz = azz(1:obj.pad:end);
                    dyF = (1i*k).^3.*yF;
                    
                    %azzz = real(ifft(dyF,obj.n*obj.pad));
                    azzz = real(ifft(dyF));
                    azzz = azzz(1:obj.pad:end);
                    dyF = (1i*k).^4.*yF;
                    azzzz = real(ifft(dyF));
                    %azzzz = real(ifft(dyF,obj.n*obj.pad));
                    azzzz = azzzz(1:obj.pad:end);
                else
                    if dim(1) == 1
                        
                        an2 = a([end-1 end 1:end-2]);
                        an1 = a([end 1:end-1]);
                        a1 = a([2:end 1]);
                        a2 = a([3:end 1 2]);
                    else
                        
                        an2 = a([end-1 end 1:end-2],:);
                        an1 = a([end 1:end-1],:);
                        a1 = a([2:end 1],:);
                        a2 = a([3:end 1 2],:);
                    end
                    
                    % if obj.fdo ==2
                    az = obj.n/obj.L*(-1/2*an1+1/2*a1);
                    azz =(obj.n/obj.L)^2*(an1-2*a+ a1);
                    azzz= (obj.n/obj.L)^3*(-1/2*an2+an1-a1+1/2*a2);
                    azzzz =(obj.n/obj.L)^4 *(an2-4*an1+6*a-4*a1+a2);
                    %             elseif obj.fdo == 4
                    %                 an3 = a([end-2 end-1 end 1:end-3]);
                    %                 a3 = a([4:end 1 2 3]);
                    %                 az = (obj.n/obj.L)*(1/12*an2-2/3*an1+2/3*a1-1/12*a2);
                    %                 azz = (obj.n/obj.L)^2*(-1/12*an2+4/3*an1-5/2*a+4/3*a1-1/12*a2);
                    %                 azzz = (obj.n/obj.L)^3*(1/8*an3-an2+13/8*an1-13/8*a1+a2-1/8*a3);
                    %                 azzzz = (obj.n/obj.L)^4*(-1/6*an3+2*an2-13/2*an1+28/3*a-13/2*a1+2*a2-1/6*a3);
                    %             end
                end
            end
            
            
        
        
        function obj = get_h(obj,optionon,hinit)
            %performs fsolve
            options = optimoptions('fsolve','Display', 'none');
            if nargin <=2
                %hinit = obj.linear_shape;
                hinit = 1+0*obj.z;
                
            end
            if nargin >1
                if optionon ==1
                    options = optimoptions('fsolve','Display', 'iter','PlotFcn',@optimplotfirstorderopt);
                end
            end
            if obj.flow == 0
                hinit = [hinit obj.q];
            end
            [obj.h0,obj.Fval,obj.eflag,obj.out,obj.jac] = fsolve(@obj.hfun,hinit,options);
            if obj.eflag <= 0
                warning('no steady state solution found')
            end
            if obj.flow == 0
                obj.q = obj.h0(end);
                obj.h0 = obj.h0(1:end-1);
            end
            
        end
        function plot_param(obj, values, param,etatoggle)
            %Vary the amplitude of the disturbance
            obj.paramtoggle = param;
            hold on
            if nargin <4
                etatoggle = 0;
                if nargin ==1
                    
                    values = obj.getp;
                end
            end
            for val = values
                obj = obj.setp(val);
                
                obj = obj.get_h;
                if obj.eflag < 1
                    warning('No solution found for %s = %g\n',obj.params{obj.paramtoggle},val)
                else
                    obj.h = obj.h0;
                    plot(obj.a+etatoggle*obj.eta,obj.nz,'DisplayName',sprintf('$%s = %g$',obj.params{obj.paramtoggle},val))
                end
            end
            if etatoggle ==1
                plot(obj.eta,obj.z/obj.L,'DisplayName','Wall')
            end
            legend
            flip_y
            xlabel('$a$')
            ylabel('$\frac{z}{L}$')
        end
        
        
        
        
        
        
        
        
        function obj = get_fluid_features(obj)
            %             if obj.eflag<1
            %                 obj.amax = 0;
            %                 obj.amin = 0;
            %                 obj.diff = 0;
            %                 obj.amaxloc = 0;
            %                 obj.aminloc = 0;
            %                 obj.minmaxdiff = 0;
            %                 obj.maxmindiff = 0;
            %             else
            [obj.hmax,i] = max(obj.a,[],2);
            [obj.hmin,j] = min(obj.a,[],2);
            obj.diff = obj.hmax-obj.hmin;
            obj.hmaxloc = obj.z(i);
            obj.hminloc = obj.z(j);
            if obj.hmaxloc<obj.hminloc
                obj.maxmindiff = obj.hminloc-obj.hmaxloc;
                obj.minmaxdiff = obj.L-obj.hminloc+obj.hmaxloc;
            else
                obj.maxmindiff = obj.L-obj.hmaxloc +obj.hminloc;
                obj.minmaxdiff = obj.hmaxloc - obj.hminloc;
                %                 end
            end
            obj.hmaxloc = obj.unperiod(obj.hmaxloc);
            obj.hminloc = obj.unperiod(obj.hminloc);
        end
        function val = unperiod(obj,val)
            [~,loc] = findpeaks(-val);
            for i =loc
                val(i:end) = val(i:end) +obj.L;
            end
        end
        
        function obj = peak_speed(obj)
            obj.peakpos = zeros(1,length(obj.h));
            k = 1;
            m = 0;
            for i = 1:length(obj.h)
                [~,j] = findpeaks(obj.h(i,:));
                if length(j) == 1
                    if j< k
                        m = m+1;
                    end
                    obj.peakpos(i) = obj.z(j)+m*obj.L;
                    k = j;
                else
                    l = 1;
                    while l<=length(j)
                        if j(l)>= k
                            obj.peakpos(i) = obj.z(j(l))+m*obj.L;
                            k = j(l);
                            break
                        end
                        l = l + 1 ;
                        if l == length(j) + 1
                            m = m+1;
                            obj.peakpos(i) = obj.z(j(1))+m*obj.L;
                            k = j(1);
                        end
                    end
                end
            end
            obj.peakspeed = diff(obj.peakpos)/obj.delt;
        end
        
        function plot_mass(obj)
            plot(obj.t,obj.mass)
            xlabel('time')
            ylabel('mass')
        end
        
        
        
        
        function plot_minmax(obj,x,param)
            obj.paramtoggle = param;
            max = [];
            min = [];
            diff = [];
            maxloc = [];
            minloc = [];
            maxmindiff = [];
            minmaxdiff = [];
            for X = x
                obj = obj.setp(X);
                
                obj = obj.reset;
                obj = obj.get_h;
                obj = obj.get_fluid_features;
                max = [max obj.hmax];
                min = [min obj.hmin];
                diff = [diff obj.diff];
                maxloc = [maxloc obj.hmaxloc/obj.L];
                minloc = [minloc obj.hminloc/obj.L];
                maxmindiff = [maxmindiff obj.maxmindiff/obj.L];
                minmaxdiff = [minmaxdiff obj.minmaxdiff/obj.L];
            end
            
            figure,clf,hold on
            plot(x,max,'DisplayName', 'Maximum thickness')
            plot(x,min,'DisplayName', 'Minimum thickness')
            plot(x,diff,'DisplayName', 'Range of thickness')
            xlabel(obj.params{obj.paramtoggle})
            ylabel('fluid thickness')
            title('How the fluid thickness is affected')
            legend
            figure,clf,hold on
            plot(x,maxloc,'DisplayName','Location of maximum')
            plot(x,minloc,'DisplayName','Location of minimum')
            xlabel(obj.params{obj.paramtoggle})
            ylabel('Turning point locations')
            title('How the fluids phase is shifted')
            flip_y
            legend
            figure,clf,hold on
            plot(x,maxmindiff,'DisplayName','distance between maximum and minimum')
            plot(x,minmaxdiff,'DisplayName','distance between minimum and maximum')
            legend
            xlabel(obj.params{obj.paramtoggle})
            ylabel('Distance between turning points')
            title('How the wavelength sees nonlinearity')
            
            
            
        end
        
        function evolve_time(obj,T, k )
            % In this strange growth regime I consider the fluid to just
            % instataneously solidify and become the new wall. I'm not
            % convinced this is that far fetched.
            if nargin <=2
                k = 0.1;
            end
            clf, hold on
            flip_y
            obj = obj.get_a;
            
            plot(obj.eta,obj.z,'DisplayName', 'Initial Wall')
            plot(obj.eta+k*obj.a,obj.z,'DisplayName','$t = 1$')
            obj.wall_shape = 3;
            for t = 2:T
                obj.eta = obj.eta+ k*obj.a;
                obj = obj.get_a;
                if obj.eflag < 1
                    warning('No solution found at time step %i. Stopping.',t)
                    break
                end
                plot(obj.eta+k*obj.a,obj.z,'DisplayName',sprintf('$t = %i$',t))
            end
            legend
        end
        function [F, J] = tfun(obj,h,hOld)
            m = obj.n*obj.periods;
            F = zeros(1,m);
            
            [~,~,~,hzzzz] = obj.getdiv(h);
            [hz,hzz,hzzz,~] = obj.getdiv(hOld);
            %[etaz,etazz,etazzz,etazzzz] = obj.getdiv(obj.eta);
            if obj.equation == 2
                F = (h-hOld)./obj.delt+3*hOld.^2.*hz+(hOld.^3.*(hzz)*obj.R)+3*hOld.^2.*hz.*((hz)*obj.R+hzzz)+h.^3.*hzzzz;
            else
                F = (h-hOld)./obj.delt+hOld.^2.*hz+obj.ep.*(hOld.^3./(3*obj.Bo).*((obj.etazz+hzz)./obj.R^2+obj.etazzzz)+hOld.^2.*hz./obj.Bo.*((obj.etaz+hz)/obj.R^2+obj.etazzz+hzzz)-2*hOld.^3.*hz/(3*obj.R)+2/15*obj.Re*(hOld.^6.*hzz+hOld.^5.*hz.^2)+h.^3.*hzzzz/(3*obj.Bo)-2/3*obj.R*(hOld.^3.*obj.etaz+3*hOld.^2.*hz.*obj.eta));
            end%maybe turn a to aold
            %F = (h-hOld)./obj.delt+hOld.^2.*hz+obj.ep.*(hOld.^3./(3*obj.Bo).*((obj.etazz+hzz)./obj.R^2+obj.etazzzz)+hOld.^2.*hz./obj.Bo.*((obj.etaz+hz)/obj.R^2+obj.etazzz+hzzz)-2*hOld.^3.*hz/(3*obj.R)+2/15*obj.Re*(hOld.^6.*hzz+hOld.^5.*hz.^2)+h.^3.*hzzzz/(3*obj.Bo)-2/3*obj.R*(hOld.^3.*obj.etaz+3*hOld.^2.*hz.*obj.eta));%maybe turn a to aold
            
            %F = (a-aold)./obj.delt+aold.^2.*az+obj.ep.*(aold.^3./(3*obj.Bo).*((obj.etazz+azz)./obj.R^2+obj.etazzzz)+aold.^2.*az./obj.Bo.*((obj.etaz+az)/obj.R^2+obj.etazzz+azzz)+2*aold.^3.*az/(3*obj.R)+2/15*obj.Re*(aold.^6.*azz+aold.^5.*az.^2)+aold.^3.*azzzz/(3*obj.Bo))+obj.etaz.*obj.ep.*aold.^3/obj.R;
            if nargout>1
                an2 = h([end-1 end 1:end-2]);
                an1 = h([end 1:end-1]);
                a1 = h([2:end 1]);
                a2 = h([3:end 1 2]);
                
                A = diag(1/obj.delt + obj.ep/obj.Bo*((obj.n)/obj.L)^4.*((an2-4*an1+8*h-4*a1+a2).*h.^2));
                B = obj.ep/(3*obj.Bo)*((obj.n)/obj.L)^4*((diag(h(1:end-2).^3,-2)+diag(h(end-1:end).^3,m-2))-4*(diag(h(1:end-1).^3,-1)+diag(h(end).^3,m-1))-4*(diag(h(2:end).^3,1)+diag(h(1).^3,1-m))+(diag(h(3:end).^3,2)+diag(h(1:2).^3,2-m)));
                J = A + B;
            end
        end
        
        function obj = dynamics(obj,init,T)
            if nargin ==3
                obj.T = T;
            end
            t = 0: obj.delt:obj.T;
            
            options = optimoptions('fsolve','Display', 'none','SpecifyObjectiveGradient',false);
            obj.h = zeros(length(t),obj.n*obj.periods);
            
            [obj.h(1,:),obj.Fval,obj.eflag,obj.out,obj.jac] = fsolve(@(a)obj.tfun(a,init),init,options);
            if obj.eflag < 1
                error('No solution found for t = %g\n',0)
            end
            
            for i = 2:length(t)
                
                [obj.h(i,:),obj.Fval,obj.eflag,obj.out,obj.jac] = fsolve(@(a)obj.tfun(a,obj.h(i-1,:)),obj.h(i-1,:),options);
                if obj.eflag < 1
                    warning('No solution found for t = %g\n',t(i))
                    break
                    
                end
                
            end
        end
        
        function obj = odedyn(obj,init)
            tic
            if nargin == 1
                obj = obj.get_h;
                init = obj.h0 + 0.1*sin(2*pi*obj.periods*obj.nz);
                Fh = fft(init);
                Fh(value.n/2+1) = 0;
                Fh(abs(Fh)<obj.suppression) = 0;
                init = ifft(Fh);
                
            end
            if obj.usejac == 1
                %opts = odeset('Stats','on','Jacobian',@obj.get_Jac);
                sparse = diag(ones(obj.n,1))+ diag(ones(obj.n -1,1),1)+ diag(ones(obj.n -2,1),2)+diag(ones(obj.n-1,1),-1)+diag(ones(obj.n-2,1),-2) + diag(1,1-obj.n) + diag([1,1],2-obj.n) + diag(1,obj.n-1) + diag([1,1],obj.n-2);
                opts = odeset('RelTol',obj.reltol,'AbsTol',obj.abstol,'Stats','on','Jacobian',@obj.get_Jac,'JPattern',sparse);
            else
                
                opts = odeset('RelTol',obj.reltol,'AbsTol',obj.abstol,'Stats','on');
            end
            %opts = odeset('RelTol',obj.reltol,'AbsTol',obj.abstol,'Vectorized','on');
            %             if obj.flow == 0
            %                 M = diag([ones(obj.n,1) 0]);
            %                 init = [init 0];
            %
            %                 opts = odeset('Mass',M,'RelTol',1e-8,'AbsTol',1e-10);
            %             end
            if obj.notol == 0
                %obj.sol = ode15s(@obj.odefun ,[0 obj.T],init,opts);
                [obj.t,obj.h]= ode15s(@obj.odefun ,0:obj.delt:obj.T,init,opts);
                %[obj.t,obj.h]= ode15s(@obj.odefun ,[0, obj.T],init,opts);
                %[obj.t,obj.h]= ode23s(@obj.odefun ,[0,obj.T],init,opts);
            else
                %obj.sol = ode15s(@obj.odefun ,[0 obj.T],init);
                %[obj.t, obj.h] = ode15s(@obj.odefun ,0:0.05:obj.T,init);
                [obj.t, obj.h] = ode15s(@obj.odefun ,[0,obj.T],init);
            end
            
            %obj.t = obj.sol.x;
            %if obj.flow == 0
            %  obj.sol.y = obj.sol.y(1:end-1,:);
            %end
            %obj.h = obj.sol.y';
            obj.integration_time = toc;
            
        end
        function obj = kstest(obj,init)
            opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
            obj.sol = ode15s(@obj.ksfun ,[0 obj.T],init,opts);
            
            obj.t = obj.sol.x;
            obj.h = obj.sol.y';
        end
        function ht = ksfun(obj,t,h)
            h = h';
            [hz,hzz,hzzz,hzzzz ] = obj.getdiv(h);
            ht = -h.*hz -0.1*hzz- 0.1*hzzzz;
            ht = ht';
            
        end
        function ht = odefun(obj,t,h)
            %             if obj.flow == 0
            %                 l = length(h);
            %                 H = h;
            %                 h = h(1:end-1);
            %             end
            h  =h';
            Fh = fft(h);
            Fh(obj.n/2+1) = 0;
            Fh(abs(Fh)<obj.suppression) = 0;
            h = ifft(Fh);
            if obj.force_mass == 1
                Fh = fft(h);
                Fh(1) = obj.n;
                h = ifft(Fh);
            end
            
            [obj.hz,obj.hzz,obj.hzzz,obj.hzzzz ] = obj.getdiv(h);

            if obj.equation == 1
                %small Bond
                ht = -obj.hz.*h.^2 -(h.^3/(3*obj.Bo).*((obj.hzz+obj.etazz)/obj.R^2+obj.hzzzz+obj.etazzzz)+obj.hz.*h.^2/obj.Bo.*((obj.hz+obj.etaz)/obj.R^2+obj.hzzz+obj.etazzz));
            elseif obj.equation == 2
                %
                ht = -3*obj.hz.*h.^2 -(h.^3.*((obj.hzz)*obj.R+obj.hzzzz)+3*obj.hz.*h.^2.*(obj.hz*obj.R+obj.hzzz));
            elseif obj.equation ==4
                ht = -h.*obj.hz +obj.R*obj.hzz;
            elseif obj.equation == 5
                ht = -obj.hzz-obj.R*obj.hzzzz-h.*obj.hz;
            elseif obj.equation == 6
                ht = -3*h.^2.*obj.hz.^2 - h.^3.*obj.hzz - obj.hz.*obj.hzzz-h.*obj.hzzzz;
            elseif obj.equation == 10
                %a equation
                ht = -obj.hz.*h.^2 -obj.ep*(2/15*obj.Re*(h.^6.*obj.hzz+6*h.^5.*obj.hz.^2)+h.^3/(3*obj.Bo).*((obj.hzz+obj.etazz)/obj.R^2+obj.hzzzz+obj.etazzzz)+obj.hz.*h.^2/obj.Bo.*((obj.hz+obj.etaz)/obj.R^2+obj.hzzz+obj.etazzz)+h.^3.*obj.hz/(3*obj.R) + h.^3.*obj.etaz/(3*obj.R));
            else
                
                %main equation
                ht = -obj.hz.*h.^2 -obj.ep*(2/15*obj.Re*(h.^6.*obj.hzz+6*h.^5.*obj.hz.^2)+h.^3/(3*obj.Bo).*((obj.hzz+obj.etazz)/obj.R^2+obj.hzzzz+obj.etazzzz)+obj.hz.*h.^2/obj.Bo.*((obj.hz+obj.etaz)/obj.R^2+obj.hzzz+obj.etazzz)-2*h.^3.*obj.hz/(3*obj.R)-2*obj.hz.*h.^2.*obj.eta/obj.R - 2*h.^3.*obj.etaz/(3*obj.R));
                
            end
            ht = ht';
            
            
        end
        function Jac = get_Jac(obj,t,h)
            
             h = h';
            
            [hz,hzz,hzzz,hzzzz] = obj.getdiv(h);
            ddh = -2*hz.*h -obj.ep*(2/15*obj.Re*(6*h.^5.*hzz+30*h.^4.*hz.^2) +h.^2/(obj.Bo).*((hzz+obj.etazz)/obj.R^2+hzzzz+obj.etazzzz)+2*hz.*h/obj.Bo.*((hz+obj.etaz)/obj.R^2+hzzz+obj.etazzz)-2*h.^2.*hz/(obj.R)-4*hz.*h.*obj.eta/obj.R - 2*h.^2.*obj.etaz/(obj.R));
            ddhz = -h.^2 -obj.ep*(2/15*obj.Re*(12*h.^5.*hz)+h.^2/obj.Bo.*((hz+obj.etaz)/obj.R^2+hzzz+obj.etazzz)+hz.*h.^2/obj.Bo.*(1/obj.R^2)-2*h.^3/(3*obj.R)-2*h.^2.*obj.eta/obj.R);
            ddhzz =  -obj.ep*(2/15*obj.Re*(h.^6)+h.^3/(3*obj.Bo).*((1)/obj.R^2));
            ddhzzz =  -obj.ep*(hz.*h.^2/obj.Bo);
            ddhzzzz = -obj.ep*(h.^3/(3*obj.Bo));
            d0 = ddh - 2*(obj.n/obj.L)^2*ddhzz+6*(obj.n/obj.L)^4*ddhzzzz;
            dn1 = -1/2*(obj.n/obj.L)*ddhz + (obj.n/obj.L)^2*ddhzz + (obj.n/obj.L)^3*ddhzzz - 4*(obj.n/obj.L)^4*ddhzzzz;
            dn2 = -1/2*(obj.n/obj.L)^3*ddhzzz + (obj.n/obj.L)^4*ddhzzzz;
            d1 = 1/2*(obj.n/obj.L)*ddhz+ (obj.n/obj.L)*ddhzz-(obj.n/obj.L)^3*ddhzzz - 4*(obj.n/obj.L)^4*ddhzzzz;
            d2 = 1/2*(obj.n/obj.L)^3*ddhzzz+(obj.n/obj.L)^4*ddhzzzz;
            Jac = diag(d0)+ diag(d1(1:end-1),1)+ diag(d2(1:end-2),2)+diag(dn1(2:end),-1)+diag(dn2(3:end),-2) + diag(d1(end),1-obj.n) + diag(d2(end-1:end),2-obj.n) + diag(dn1(1),obj.n-1) + diag(dn2(1:2),obj.n-2);
            
        end
        
        function animate(obj,speed,wall,c,clearf)
            if nargin <= 4
                clearf = 0 ;
            end
            if nargin <= 3
                c = 0;
            end
            if nargin <= 2
                wall = 1;
                
            end
            if nargin <=1
                speed = 1;
            end
            data = obj.a;
            l = length(data(:,1));
            %
            eta = repmat(obj.eta,length(obj.t),1);
            if c~= 0
                phi = mod((obj.z-c*obj.t),obj.L);
                [~,I] = sort(phi,2);
                I = I + obj.n*(0:length(I)-1)';
                data = data';
                data = data(I);
                
                eta = eta';
                eta = eta(I);
                
            end
                
            %obj = obj.follow_peak;
                   
                grid
            for i = 1:speed:l
                if clearf == 0
                   clf
                end
                hold on
                
                if (obj.wall_shape ~= 3 && wall ~=0)
                    plot(obj.z,eta(i,:))
                end
                
                    
                plot(obj.z,data(i,:)+wall*eta(i,:))
                %plot(1,mod(obj.zpos(i),obj.L),'>')
                
                
                
                
                
                
                title(sprintf('$t =%g$',obj.t(i)))
                
                %ylim([wall*-obj.del+(1-wall)*0.5,1.5+wall*obj.del])
                if clearf ==0
                pause(0.01)
                end
                pause(0.01)
            end
        end
        
        
        function stabcheck(obj,del)
            figure(10), clf, hold on, title("Amplitude Changes with time")
            figure(11), clf, hold on, title('Position of maxima/minima')
            figure(12), clf, hold on , title('Drift of results')
            for d = del
                obj.del = d;
                obj = obj.get_a;
                obj = obj.dynamics(100,obj.a+ 0.1*sin(obj.z));
                obj = obj.get_fluid_features;
                t = 0:obj.delt:100;
                figure(10)
                plot(t,obj.diff, 'DisplayName',sprintf('$\\delta = %g$',d))
                figure(11)
                plot(t,obj.amaxloc, 'DisplayName',sprintf('maxima $ \\delta = %g$',d))
                fig = gca;
                plot(t,obj.aminloc,'--','Color', fig.Children(1).Color, 'DisplayName',sprintf('minima$ \\delta = %g$',d))
                figure(12)
                plot(t,obj.amax, 'DisplayName',sprintf('$\\delta = %g$',d))
            end
        end
        function plot_features(obj,m)
            if nargin ==1
                m= 0;
            end
            obj = obj.get_fluid_features;
            figure(10+m),  hold on, title("Amplitude Changes with time")
            figure(11+m),  hold on, title('Position of maxima/minima')
            figure(12+m),  hold on , title('Drift of results')
            d = obj.del;
            
            figure(10+m)
            
            plot(obj.t,obj.diff, 'DisplayName',sprintf('$\\delta = %g$',d))
            figure(11+m)
            plot(obj.t,obj.hmaxloc/obj.L, 'DisplayName',sprintf('maxima $ \\delta = %g$',d))
            fig = gca;
            %plot(t,obj.hminloc,'--','Color', fig.Children(1).Color, 'DisplayName',sprintf('minima$ \\delta = %g$',d))
            figure(12+m)
            plot(obj.t,obj.hmax, 'DisplayName',sprintf('$\\delta = %g$',d))
        end
        
        function big_wall_effect(obj,Bo,R)
            if nargin ==1
                Bo = obj.Bo;
                R = obj.R;
            end
            clf, hold on
            plot(obj.eta,obj.z,'DisplayName','Wall')
            for b = Bo
                for r = R
                    for i =1:obj.n
                        [etaz,etazz,~,etazzzz] = obj.getdiv(obj.eta,i);
                        wall(i) = 1/(3*b)*(etazz/r^2 +etazzzz) +etaz/r;
                    end
                    plot(wall,obj.z, 'DisplayName',sprintf('R = %g, Bo = %g',r,b))
                    
                end
            end
            legend
        end
        function [F,J] = small_fluid_disturbance_fun(obj,a,aold)
            
            F = zeros(1,obj.n);
            
            for i = 1:obj.n
                [az,azz,~,azzzz] = obj.getdiv(a,i);
                [aoz,aozz,~,aozzzz] = obj.getdiv(a,i);
                [etaz,~,etazzz,~] = obj.getdiv(obj.eta,i);
                
                F(i) = (a(i)-aold(i))/obj.delt+1/2*(2*az*a(i)+2/(3*obj.R)* az+2/15*azz+1/(3*obj.Bo)*((azz+3*az*etaz)/obj.R^2 + azzzz+3*az*etazzz)+3*a(i)*etaz/obj.R)+1/2*(2*aoz*aold(i)+2/(3*obj.R)* aoz+2/15*aozz+1/(3*obj.Bo)*((aozz+3*aoz*etaz)/obj.R^2 + aozzzz+3*aoz*etazzz)+3*aold(i)*etaz/obj.R);
            end
            if nargout>1
                an2 = a([end-1 end 1:end-2]);
                an1 = a([end 1:end-1]);
                a1 = a([2:end 1]);
                a2 = a([3:end 1 2]);
                etan2 = obj.eta([end-1 end 1:end-2]);
                etan1 = obj.eta([end 1:end-1]);
                eta1 = obj.eta([2:end 1]);
                eta2 = obj.eta([3:end 1 2]);
                etaz = obj.n/obj.L*1/2*(eta1-etan1);
                etazzz = (obj.n/obj.L)^3 * (-1/2*etan2+etan1-eta1+1/2*eta2);
                
                A = diag(1/obj.delt + 1/2*(2*obj.n/obj.L)*(a1-an1) -4/15*obj.Re*(obj.n/obj.L)^2+1/(3*obj.Bo)*(-2*(obj.n/obj.L)^2/obj.R^2+6*(obj.n/obj.L)^4)+3*etaz/obj.R);
                B = 1/2*((1/(3*obj.Bo)*(obj.n/obj.L)^4)*((diag(ones(1,obj.n-2),-2)+diag([1,1],obj.n-2))+((-a-1/(3*obj.R)-1/(2*obj.Bo*obj.R^2)*etaz-1/(2*obj.Bo)*etazzz)*(obj.n/obj.L)+(2*obj.Re/15+1/(3*obj.Bo*obj.R^2))*(obj.n/obj.L)^2)*(diag(ones(1,obj.n-1),-1)+diag(1,obj.n-1))+((a+1/(3*obj.R)+1/(2*obj.Bo*obj.R^2)*etaz+1/(2*obj.Bo)*etazzz)*(obj.n/obj.L)+(2*obj.Re/15+1/(3*obj.Bo*obj.R^2))*(obj.n/obj.L)^2)*(diag(ones(1,obj.n-1),1)+diag(1,1-obj.n))+(1/(3*obj.Bo)*(obj.n/obj.L)^4)*(diag(ones(1,obj.n-2),2)+diag([1,1],2-obj.n))));
                J = A + B;
            end
            
        end
        function obj = small_dynamics(obj,T,init)
            t = 0: obj.delt:T;
            
            options = optimoptions('fsolve','Display', 'none','SpecifyObjectiveGradient',false);
            obj.a = zeros(length(t),obj.n);
            
            [obj.a(1,:),obj.Fval,obj.eflag,obj.out,obj.jac] = fsolve(@(a)obj.small_fluid_disturbance_fun(a,init),init,options);
            if obj.eflag < 1
                error('No solution found for t = %g\n',0)
            end
            
            for i = 2:length(t)
                
                [obj.a(i,:),obj.Fval,obj.eflag,obj.out,obj.jac] = fsolve(@(a)obj.small_fluid_disturbance_fun(a,obj.a(i-1,:)),obj.a(i-1,:),options);
                if obj.eflag < 1
                    warning('No solution found for t = %g\n',t(i))
                    break
                    
                end
                
            end
            
        end
        
        function h2varyamp(obj)
            del = 0.05:0.05:2;
            for i = 1:length(del)
                obj.del = del(i);
                obj = obj.get_h;
                if obj.eflag<= 0
                    break
                end
                obj.h = obj.h0;
                h2(i) = obj.h2norm;
            end
            l = length(h2);
            plot(del(1:l),h2)
        end
        function plot_h2varyamp(obj,vals,param)
            if nargin ==3
                obj.paramtoggle =param;
            end
            hold on
            i = 1;
            for val = vals
                obj = obj.setp(val);
                obj.h2varyamp;
                text{i} = sprintf('$%s = %g$',obj.params{obj.paramtoggle},val);
                i = i+1;
            end
            xlabel('Amplitude')
            ylabel('$||h||_2$')
            title(sprintf('$||h||_2$ as increasing amplitude for different %s',obj.params{obj.paramtoggle}))
            
            legend(text)
        end
        function obj = find_n_peaks(obj)
            for i = 1:length(obj.t)
                pks = findpeaks(obj.h(i,:));
                obj.npks(i) = length(pks);
            end
        end
        
        function obj = new_get_a(obj,optionon,ainit)
            %performs fsolve
            options = optimoptions('fsolve','Display', 'none');
            if nargin <=2
                ainit = obj.linear_shape;
                
            end
            if nargin >1
                if optionon ==1
                    options = optimoptions('fsolve','Display', 'iter','PlotFcn',@optimplotfirstorderopt);
                end
            end
            
            [obj.a,obj.Fval,obj.eflag,obj.out,obj.jac] = fsolve(@obj.new_afun,ainit,options);
            
        end
        
        function mass = controlmass(obj,q)
            
            hinit = obj.linear_shape;
            options = optimoptions('fsolve','Display', 'none');
            h = fsolve(@(h) obj.hfunflux(h,q),hinit,options);
            mass = sum(h,2)/obj.n;
            
        end
        function [q,h2,range] = steadynorms(obj,range,param)
            if nargin<3
                param = obj.paramtoggle;
                if nargin <2
                    range = 0.1:0.1:1.5;
                end
            end
            obj.paramtoggle = param;
            for i = 1:length(range)
                obj = obj.setp(range(i));
                obj = obj.get_q;
                q(i) = obj.q;
                obj = obj.get_h;
                obj.h = obj.h0;
                h2(i) = obj.h2norm;
            end
        end
        
        
        
        function obj = get_q(obj,mass)
            if nargin <2
                mass = 1;
            end
            [obj.q,~,eflag,~] = fzero(@(q) obj.controlmass(q) -mass,1/3);
            if eflag<1
                obj.q = -1;
            end
        end
        
        function [F1,F2] = check_h0(obj)
            [hz,hzz,hzzz,hzzzz] = obj.getdiv(obj.h0);
            F1 = obj.h0.^3/3 + obj.ep*(2/15*obj.Re*obj.h0.^6.*hz+obj.h0.^3/(3*obj.Bo).*((hz+obj.etaz)/obj.R^2+ hzzz+ obj.etazzz)-obj.h0.^4/(12*obj.R)-2/(3*obj.R)*obj.h0.^3.*obj.eta);
            F2 =   hz.*obj.h0.^2 + obj.ep.*(-4*hz.*obj.h0.^3/(12.*obj.R)+2/15.*obj.Re.*(obj.h0.^6.*hzz+6.*obj.h0.^5.*hz.^2)+1/(3.*obj.Bo).*(obj.h0.^3.*((hzz+obj.etazz)/obj.R^2+hzzzz+obj.etazzzz)+3*hz.*obj.h0.^2.*((hz+obj.etaz)/obj.R^2+hzzz+obj.etazzz))-2*obj.etaz.*obj.h0.^3/(3*obj.R)-2*hz.*obj.h0.^2.*obj.eta/obj.R);
        end
        function ploth2(obj,diff)
            if nargin == 1
                diff = 0;
            end
            if diff == 0
                plot(obj.t,obj.h2norm)
            else
                plot(obj.t,sum(obj.hdiff.^2/obj.n,2))
            end
        end
        
        function plotz0(obj,m)
            % m = [0 1 2 3 10] plots the position at the maixima, 0 minima and 0
            if nargin ==1
                m = 10;
            end
            
            if m == 10
                plot(obj.t,mean(obj.a,2))
            else
                plot(obj.t(floor(3/4*end):end), obj.a(floor(3/4*end):end,obj.n/4*m+1))
            end
            xlabel('$t$')
            ylabel('$a$')
            title('Fluid thickness over the maximum')
        end
        
        function phaseplot(obj,m1,m2,npeak)
            if nargin ==3
                hold on
                %plot(obj.h(1:floor(3/4*end),obj.n/4*m1+1),obj.h(1:floor(3/4*end),obj.n/4*m2+1),'--','Color',[0.8 0.8 0.8])
                plot(obj.h(floor(3/4*end):end,obj.n/4*m1+1),obj.h(floor(3/4*end):end,obj.n/4*m2+1),'Color',[0 0.4470 0.7410])
                
                hold off
                xlabel('$z_{max}$')
                ylabel('$z_{min}$')
                title(sprintf('Trajectory over the maximum vs minimum for $L = %g\\pi$, $\\delta = %g$',obj.L/pi,obj.del))
            else
                lasttrajects = obj.loc(obj.goodP);
                plot(obj.h(lasttrajects(npeak):end,obj.n/4*m1+1),obj.h(lasttrajects(npeak):end,obj.n/4*m2+1))
                % plot(obj.h(lasttrajects(npeak):end,obj.n/4*m1+1)./obj.mass(lasttrajects(npeak):end),obj.h(lasttrajects(npeak):end,obj.n/4*m2+1)./obj.mass(lasttrajects(npeak):end))
            end
        end
        function phaseplota(obj,m1,m2,npeak)
            if nargin ==3
                hold on
                plot(obj.a(1:floor(3/4*end),obj.n/4*m1+1),obj.a(1:floor(3/4*end),obj.n/4*m2+1),'--','Color',[0.8 0.8 0.8])
                plot(obj.a(floor(3/4*end):end,obj.n/4*m1+1),obj.a(floor(3/4*end):end,obj.n/4*m2+1),'Color',[0 0.4470 0.7410])
                
                hold off
                xlabel('$z_{max}$')
                ylabel('$z_{min}$')
                title(sprintf('Trajectory over the maximum vs minimum for $L = %g\\pi$, $\\delta = %g$',obj.L/pi,obj.del))
            else
                lasttrajects = obj.loc(obj.goodP);
                plot(obj.h(lasttrajects(npeak):end,obj.n/4*m1+1),obj.a(lasttrajects(npeak):end,obj.n/4*m2+1))
                % plot(obj.h(lasttrajects(npeak):end,obj.n/4*m1+1)./obj.mass(lasttrajects(npeak):end),obj.h(lasttrajects(npeak):end,obj.n/4*m2+1)./obj.mass(lasttrajects(npeak):end))
            end
        end
        function peakpoints(obj,m)
            if nargin == 1
                m = 0;
            end
            [pks,loc,W, P ]  = findpeaks(obj.a(:,m*obj.n/4+1));
            hold on
            %plot(obj.t(loc),pks);
            %plot(obj.t(loc),W)
            %plot(obj.t(loc),P)
            meanP = mean(P);
            goodP = P>meanP;
            plot(obj.t(loc(goodP)),pks(goodP))
            goodloc = loc(goodP);
            
            
            
            %plot(obj.t(goodloc(1:end-1)),diff(obj.t(loc(goodP))))
            periods = diff(obj.t(loc(goodP)));
            %lper = mean(periods(end-10:end));
        end
        function obj = get_peak_data(obj)
            
            
            [obj.pks,obj.loc,obj.W, obj.P ]  = findpeaks(obj.a(:,obj.peakloc*obj.n/4+1));
            meanP = mean(obj.P);
            
            obj.goodP = obj.P>meanP;
        end
        
        function followminmax(obj)
            [~,imax] = max(obj.h,[],2);
            [~,imin] = min(obj.h,[],2);
            hold on
            plot(obj.nz(imax(1:floor(3/4*end))),obj.nz(imin(1:floor(3/4*end))),'--','Color',[0.8 0.8 0.8])
            plot(obj.nz(imax(floor(3/4*end):end)),obj.nz(imin(floor(3/4*end):end)),'Color',[0 0.4470 0.7410])
            xlabel('$h_{max}$')
            ylabel('$h_{min}$')
            title(sprintf('Normalised postion of the maximum vs minimum for $L = %g\\pi$, $\\delta = %g$',obj.L/pi,obj.del))
        end
        function followhminmax(obj)
            hmax = max(obj.h,[],2);
            hmin = min(obj.h,[],2);
            hold on
            plot(hmax(1:floor(3/4*end)),hmin(1:floor(3/4*end)),'--','Color',[0.8 0.8 0.8])
            plot(hmax(floor(3/4*end):end),hmin(floor(3/4*end):end),'Color',[0 0.4470 0.7410])
            xlabel('$h_{max}$')
            ylabel('$h_{min}$')
            title(sprintf('Normalised postion of the maximum vs minimum for $L = %g\\pi$, $\\delta = %g$',obj.L/pi,obj.del))
        end
        function pkl = peaklength(obj)
            mid = (max(obj.h(:,1)) + min(obj.h(:,1)))/2;
            toppeaks = obj.pks > mid;
            
            tpos = obj.t(obj.loc(toppeaks));
            try
                pkl = mean(diff(tpos(floor(end/2):end)));
            catch
                obj = obj.get_peak_data;
                saven500(obj)
                pkl = obj.peaklength;
            end
            
        end
        function plot_trajectories(obj,t0, tint,tend,m)
            if nargin < 3
                tint = 1;
                tend = t0+10;
                m = 0.2;
            end
            tvec = [];
            hold on,
            for t = t0:tint:tend
                [tval,i] = min(abs(obj.t-t));
                
                plot(obj.z,obj.a(i,:)+m*(t)/tint)
                
            end
            ylabel(sprintf('$%g(t-%g+a$)',m/tint,t0))
            xlabel('$z$')
            title('Fluid thickness as time progresses')
            
        end
        
        function obj = follow_peak(obj,t,zreset,region)
            if nargin == 1
                t = 0;
                zreset = 0;
                region = 10;
            end
            [ba,i2] = min(abs(obj.t -t));
            l = length(obj.t(1:end));
            
            m = 0 ;
            [~, k ] = max(obj.a(i2,:));
            
            z0 = obj.z(k);
            obj.zpos(1) = 0;
            for i = 1:l-i2
                if k+region<=obj.n && k-region>0
                    [~, j ] = max(obj.a(i+i2,k-region:k+region));
                    k = k-region + j - 1;
                elseif k+region>obj.n
                    [a, j ] = max(obj.a(i+i2,k-region:obj.n));
                    
                    [b,j2] = max(obj.a(i+i2,1:k+region-obj.n));
                    if a>b
                        k = k-region+j-1;
                    else
                        k = j2;
                        m= m+1;
                    end
                    
                    
                elseif k-20<1
                    [a, j ] = max(obj.h(i+i2,k-region+obj.n:obj.n));
                    [b,j2] = max(obj.h(i+i2,1:k+region));
                    if a > b
                        m= m-1;
                        k = k-region+obj.n+j-1;
                    else
                        k = j2;
                    end
                    
                end
                
                
                obj.zpos(i+1) = obj.z(k) + m*obj.L-zreset*z0;
                %                 if k ==1
                %                     plot(obj.t(i+i2),obj.zpos(i+1),'kx')
                %                 elseif k== 251
                %                      plot(obj.t(i+i2),obj.zpos(i+1),'ko')
                %                 end
                obj.speed(i) = (obj.zpos(i+1) - obj.zpos(i))/(obj.t(i2+i)-obj.t(i2+i-1));
            end
            
            plot(obj.t(i2:end),obj.zpos)
            
        end
        function c = get_c(obj,method)
            if method == 0
                obj = obj.follow_peak;
                c = polyfit(obj.t(floor(0.4*end):end),obj.zpos(floor(0.4*end):end),1);
                c = c(1);
            else
                obj = obj.get_peak_data;
                goodpeaks  = obj.loc(floor(0.4*end):end);
                c = obj.L *(length(goodpeaks)-1)/(obj.t(goodpeaks(end))- obj.t(goodpeaks(1)));
            end
        end
        function surfdata(obj,c)
            if nargin ==1
            obj = obj.follow_peak;
            c = polyfit(obj.t(floor(0.2*end):end),obj.zpos(floor(0.2*end):end),1);
            end
            phi = mod((obj.z-c(1)*obj.t),obj.L);
            [phi,I] = sort(phi,2);
            I = I + obj.n*(0:length(I)-1)';
%             for j = 1:length(obj.t)/100
%                 a(100*j,:) = obj.a(100*j,I(j,:));
%             end
            
            clf
            a = obj.a';
            pcolor(phi,obj.t,a(I))
          
             
          

            shading interp
            
            colorbar
            xlabel('$z-ct$')
            ylabel('$t$')
            del = erase(string(obj.del),'.');
            title(sprintf('Fluid Thickness for wall amplitude = %g, c = %g',obj.del, c(1)))
            name = sprintf('../plots/colour/ColourL%gdel%s',obj.L/pi,del);
%             savefig(name)
%             saveas(gcf,name,'epsc')
%             
        end
        
        
        function obj = eliminate_noise(obj)
            
        Fh = fft(obj.h,[],2);
        Fh(obj.n/2+1) = 0;
        obj.h = real(ifft(Fh,[],2));
        end
    end
end
%% Old Code

%         INCORPORATED INTO get_wall
%                 function obj = make_wall_step(obj,l,del2)
%             obj.eta = obj.del/2*(tanh((obj.z-l/2*obj.L)/del2) - tanh((obj.L*(l/2-1)+obj.z)/del2));
%         end
%         function obj = make_wall_saw(obj,l,del2)
%             obj.eta = (obj.z-l/2*obj.L).*obj.del./2.*(tanh((obj.z-l/2*obj.L)/del2) - tanh((obj.L*(l/2-1)+obj.z)/del2));
%         end

%         function F = afun(obj,a)
%
%             %The system of equations for fsolve to solve
%             %Not really set up for flow = 0
%
%             F = 0*(1:obj.n);
%
%             if obj.boundary == 1
%                 F(1) = a(1)-1;
%                 F(obj.n) = a(obj.n)-1;
%                 F(2) = 2*a(2)-3/2*a(1)-1/2*a(3);
%                 F(obj.n - 1) = 3/2*a(obj.n)-2* a(obj.n-1)+1/2*a(obj.n-2);
%             end
%             for m = 2*obj.boundary+1:obj.n-obj.boundary*2
%                 [az,azz,azzz,azzzz] = obj.getdiv(a,m);
%                 if obj.wall_shape >0
%
%                     [etaz,~,etazzz,~] = obj.getdiv(obj.eta,m);
%                 else
%                     etaz = -obj.del*2*pi/obj.L*sin(2*pi/obj.L*obj.z(m));
%                     etazzz = obj.del*8*pi^3/obj.L^3*sin(2*pi/obj.L*obj.z(m));
%                 end
%
%                 if obj.flow ==1
%                     F(m) = a(m)^3/3+obj.ep*(a(m)^4/3+a(m)^3*obj.eta(m)+a(m)^3/(3*obj.Bo)*(az+azzz+etaz+ etazzz)-3/40*obj.Re*a(m)^6*az)-obj.q;
%                 else
%
%                     F(m) = a(m)^2*az + obj.ep*(4/3*a(m)^3*az+3*a(m)^2*az*obj.del*cos(2*pi/obj.L*obj.z(m))-2*pi/obj.L*obj.del*sin(2*pi/obj.L*obj.z(m))*a(m)^3+a(m)^2*az/obj.Bo*(az+azzz+(8*pi^3/obj.L^3-2*pi/obj.L)*obj.del*sin(2*pi/obj.L*obj.z(m)))+a(m)^3/(3*obj.Bo)*(azz+azzzz+(16*pi^4/obj.L^4 - 4*pi^2/obj.L^2)*obj.del*cos(2*pi/obj.L*obj.z(m)))-3/40*obj.Re*(a(m)^6*azz+6*a(m)^5*az^2));
%                     F(obj.n+1+obj.flow*(m-1))  = F(obj.n+1+obj.flow*(m-1))  + (1-obj.flow)*(a(m) + obj.ep*(a(m)^2/2 + a(m)*obj.del* cos(2*pi/obj.L*obj.z(m))))*obj.L/obj.n+obj.flow*(a(m)^3 + obj.ep*(a(m)^4/3+a(m)^3*obj.del*cos(2*pi/obj.L*obj.z(m))+a(m)^3/(3*obj.Bo)*(az+azzz-2*pi/obj.L*obj.del*sin(obj.z(m)*2*pi/obj.L)+8*pi^3/obj.L^3*obj.del*sin(2*pi/obj.L*obj.z(m))-3/40*obj.Re*a(m)^6*az))-obj.q);
%
%                 end
%             end
%
%             if obj.flow == 0
%
%                 F(obj.n+1) = F(obj.n+1)  - (1-obj.flow)*(1+obj.ep/2)*obj.L;
%                 %F = F(obj.n+1:2*obj.n);
%             end
%
%         end
%
%         function [az,azz,azzz,azzzz] = getdiv(obj,a,m)
%             %gets the derivatives for afun to use
%             if m ==1
%                 az = obj.n/obj.L*(-1/2*a(obj.n)+1/2*a(m+1));
%                 azz =(obj.n/obj.L)^2*(a(obj.n)-2*a(m)+ a(m+1));
%                 azzz= (obj.n/obj.L)^3*(-1/2*a(obj.n-1)+a(obj.n)-a(m+1)+1/2*a(m+2));
%                 azzzz =(obj.n/obj.L)^4 *(a(obj.n-1)-4*a(obj.n)+6*a(m)-4*a(m+1)+a(m+2));
%             elseif m==2
%                 az = obj.n/obj.L*(-1/2*a(m-1)+1/2*a(m+1));
%                 azz =(obj.n/obj.L)^2*(a(m-1)-2*a(m)+ a(m+1));
%                 azzz= (obj.n/obj.L)^3*(-1/2*a(obj.n)+a(m-1)-a(m+1)+1/2*a(m+2));
%                 azzzz =(obj.n/obj.L)^4 *(a(obj.n)-4*a(m-1)+6*a(m)-4*a(m+1)+a(m+2));
%             elseif m ==obj.n
%                 az = obj.n/obj.L*(-1/2*a(m-1)+1/2*a(1));
%                 azz =(obj.n/obj.L)^2*(a(m-1)-2*a(m)+ a(1));
%                 azzz= (obj.n/obj.L)^3*(-1/2*a(m-2)+a(m-1)-a(1)+1/2*a(2));
%                 azzzz =(obj.n/obj.L)^4 *(a(m-2)-4*a(m-1)+6*a(m)-4*a(1)+a(2));
%             elseif m==obj.n-1
%                 az = obj.n/obj.L*(-1/2*a(m-1)+1/2*a(m+1));
%                 azz =(obj.n/obj.L)^2*(a(m-1)-2*a(m)+ a(m+1));
%                 azzz= (obj.n/obj.L)^3*(-1/2*a(m-2)+a(m-1)-a(m+1)+1/2*a(1));
%                 azzzz =(obj.n/obj.L)^4 *(a(m-2)-4*a(m-1)+6*a(m)-4*a(m+1)+a(1));
%             else
%                 az = obj.n/obj.L*(-1/2*a(m-1)+1/2*a(m+1));
%                 azz =(obj.n/obj.L)^2*(a(m-1)-2*a(m)+ a(m+1));
%                 azzz= (obj.n/obj.L)^3*(-1/2*a(m-2)+a(m-1)-a(m+1)+1/2*a(m+2));
%                 azzzz =(obj.n/obj.L)^4 *(a(m-2)-4*a(m-1)+6*a(m)-4*a(m+1)+a(m+2));
%             end
%         end

%         function compare_constant(obj)
%             %Compare the method of keeping the mass constant, with the
%             %method of keeping the flow rate constant
%             obj = obj.reset;
%             obj.flow = 0;
%             obj = obj.get_a;
%             clf, hold on
%             if obj.eflag <1
%                 warning('No solution found for constant mass')
%             else
%                 plot(obj.a, obj.z)
%             end
%             obj.flow = 1;
%             obj = obj.get_a;
%             if obj.eflag <1
%                 warning('No solution found for constant flow rate')
%             else
%                 plot(obj.a, obj.z)
%             end
%
%         end
%         function [F, J] = tfunbad(obj,a,aold)
%             F = zeros(1,obj.n);
%             for i = 1:obj.n
%                 [~,~,~,azzzz] = obj.getdiv(a,i);
%                 [az,azz,azzz,~] = obj.getdiv(aold,i);
%                 [etaz,etazz,etazzz,etazzzz] = obj.getdiv(obj.eta,i);
%                 F(i) = (a(i)-aold(i))/obj.delt+aold(i)^2*az+obj.ep*(aold(i)^3/(3*obj.Bo)*(etazz+azz+etazzzz)+aold(i)^2*az/obj.Bo*(etaz+az+etazzz+azzz)-aold(i)^4/6-3/40*obj.Re*aold(i)^6*az+a(i)^3*azzzz/(3*obj.Bo));
%             end
%             if nargout>1
%                 an2 = a([end-1 end 1:end-2]);
%                 an1 = a([end 1:end-1]);
%                 a1 = a([2:end 1]);
%                 a2 = a([3:end 1 2]);
%          com
%                 A = diag(1/obj.delt + obj.ep/obj.Bo*((obj.n)/obj.L)^4.*((an2-4*an1+8*a-4*a1+a2).*a.^2));
%                 B = obj.ep/(3*obj.Bo)*((obj.n)/obj.L)^4*((diag(a(1:end-2).^3,-2)+diag(a(end-1:end).^3,obj.n-2))-4*(diag(a(1:end-1).^3,-1)+diag(a(end).^3,obj.n-1))-4*(diag(a(2:end).^3,1)+diag(a(1).^3,1-obj.n))+(diag(a(3:end).^3,2)+diag(a(1:2).^3,2-obj.n)));
%                 J = A + B;
%             end
%         end
%
%
%


