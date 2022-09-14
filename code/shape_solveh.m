classdef shape_solveh
    
    %Class to help solve Thin film down a cylinder of similar size to the
    %wavelength
    
    properties
        n {mustBeInteger,mustBePositive} = 256  %number of grid points INT
        nr = 100
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
        gamma {mustBeNonnegative} = 1e-7 % Growth rate of wall
        del {mustBeNonnegative} = 1 % Disturbance amplitude
        L {mustBeNonnegative} = 2*pi % wavelength of disturbance
        l {mustBeNonnegative,mustBeLessThanOrEqual(l,1)} = 0.5 %width of the step
        del2 {mustBeNonnegative} = 0.1 % steepness of the step
        dir {mustBeMember(dir,[-1,1])} = 1 %direction of slope
        dir2 {mustBeMember(dir2,[-1,1])} = 1 %direction of step

        ps {mustBeMember(ps,[0,1])} = 1; %whether to use Psuedo spectral differentiation rather than finite differences
        % Asthetic features
        
        periods {mustBePositive,mustBeInteger} = 1 % number of periods to plot INT - ignored for now
        
        %Derived from parameters
        r %domain
        
        z % domain
        t % time domain
        h % mass conserving fluid thickness
        
        eta % wall shape
        h0 %steady state
        etat = 0 %if h is actually eta
        H = 1e-5 %mean thickness
        
        f %growth rate of wall
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
        filter_oscillation = 1
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
        peakh
        peakp
        
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
        a0
        w
        u
        U
        pz
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
                function value = get.a0(obj)
                    if isempty(obj.h0)
                        obj = obj.get_h;
                    end
                    
            if obj.equation == 0
                value = (sqrt(2*obj.h0*obj.R/obj.ep +  (obj.eta+obj.R/obj.ep).^2)-(obj.R/obj.ep+ obj.eta));
            else
                value = obj.h0;
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
                value  = obj.h.^3/3 +  (obj.h.^3/3/obj.Bo.*((hz+obj.etaz)/obj.R^2+ hzzz+ obj.etazzz));
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
            value = sqrt(sum(obj.h.^2,2)/obj.n);
        end
        
        function value = get.pz(obj)
            [hz,~,hzzz,~] = obj.getdiv(obj.h);
            
            value = -1/obj.Bo.*((hz+obj.etaz)/obj.R^2+ hzzz+ obj.etazzz);
        end
        function value = get.w(obj)
        if obj.equation == 0
                [az,~,azzz,~] = obj.getdiv(obj.a);
            value = obj.a.^2.*(obj.r - obj.r.^2/2) + obj.ep*(obj.Re*obj.a.^5.*az.*(obj.r.^4/24 -obj.r.^3/6 + obj.r/3) + 1/obj.Bo*obj.a.^2.*((az+obj.etaz)/obj.R^2 + azzz+ obj.etazzz).*(obj.r-obj.r.^2/2) + obj.a.^3/obj.R.*(obj.r.^3/6 - obj.r.^2/2 + obj.r/2));    
            elseif obj.equation == 1
                value = obj.h.^2.*(obj.pz-1).*(obj.r.^2/2-obj.r);
            end
            
                
        end
        
        function value = get.u(obj)
            [hz,~,~,~] = obj.getdiv(obj.h);
            
            if obj.equation == 0
            value = obj.a.^2;
            elseif obj.equation ==1
                [pzz,~,~,~] = obj.getdiv(obj.pz);
                value = obj.etaz.*obj.h.^2.*(obj.pz-1).*(obj.r.^2/2-obj.r) +1/2*hz.*obj.h.^2.*(obj.pz-1).*obj.r.^2-obj.h.^3.*pzz.*(obj.r.^3/6-obj.r.^2/2);
            end
        end
        function value = get.U(obj)
            [hz,~,~,~] = obj.getdiv(obj.a);
            if obj.equation ==0
                value = 3*obj.a.^2.*hz.*(obj.r.^3/6-obj.r.^2/2);
            elseif obj.equation == 1
                [pzz,~,~,~] = obj.getdiv(obj.pz);
                value = -(3*obj.a.^2.*hz.*(obj.pz - 1)+obj.a.^3.*pzz).*(obj.r.^3/6-obj.r.^2/2);
            end
            
        end

        function obj = set.W(obj,r)
        W = (obj.h).^2.*(1-obj.pz).*(r-r.^2/2);
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
        function obj = set.nr(obj,value)
            obj.nr = value;
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
            obj.q = 1/3 + 1/3*obj.ep/obj.R*(1-obj.equation);
            obj.r = (0:1/obj.nr:1)';
            
            
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
                obj.eta = obj.del*(tanh((obj.z/obj.L-(1-obj.l)/2)/obj.del2) - tanh((obj.z/obj.L - (1+obj.l)/2)/obj.del2))-2*obj.del*obj.l;
                
            elseif  obj.wall_shape ==2
                obj.eta = obj.dir*((obj.z/obj.L-(1-obj.dir2*obj.l)/2)/(obj.l).*(obj.del*(tanh((obj.z/obj.L-(1-obj.l)/2)/obj.del2) - tanh((obj.z/obj.L -(1+obj.l)/2)/obj.del2))) -obj.dir2*obj.del*obj.l);
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
            if obj.equation ==0
                a0 = 1;
                k = 2*pi/obj.L;
                R = obj.R;
                Re = obj.Re;
                Bo = obj.Bo;
                eps = obj.ep;
                %             A = -k  +obj.ep*(2*k./(3.*obj.R));
                %             B = obj.ep*(-2/15*obj.Re*k.^2+1./(3.*obj.Bo)*(-k.^2/obj.R.^2+k.^4));
                %             C = obj.ep*(2.*k./(3*obj.R));
                %             D = (obj.ep./(3.*obj.Bo).*(-k.^2./obj.R.^2+k.^4));
                %             a = -(A.*C+B.*D)./(A.^2+B.^2);
                %             b = (B.*C-D.*A)/(A.^2+B.^2);
                %
                %             hinit = 1+ a.*obj.del.*cos(2*pi./obj.L.*obj.z)-b.*obj.del.*sin(2*pi./obj.L.*obj.z);
%                 
%                 a = (5*eps.*(20*eps.*R.^2.*Bo.^2-30*Bo.^2.*R.^3+eps.*k.^2.*(k.^2.*R.^2-1).*((5*k.^2-2*Bo.*Re).*R.^2-5)))./(225*Bo.^2.*R.^4-300*Bo.^2.*eps.*R.^3+eps.^2.*(100*R.^2.*Bo.^2+k.^2.*(5+(2*Bo.*Re-5*k.^2).*R.^2).^2));
%                 
%                 theta = atan(Bo.*k.*R.^2.*(15*k.^2.*R.^2-15-4*eps.*Re.*R.*Bo)./(30*Bo.^2.*R.^3-eps.*(20*Bo.^2.*R.^2+k.^2.*(k.^2.*R.^2-1).*((5*k.^2-2.*Bo.*Re).*R.^2-5))))+(1-sign(a))/2*pi;
%                 A = 5.*eps.*sqrt((4*R.^2.*Bo.^2+(k-k.^3.*R.^2).^2)./(225*Bo.^2.*R.^4-300*eps.*Bo.^2.*R.^3+eps.^2.*(100*R.^2.*Bo.^2+k.^2.*(5 + (2*Bo.*Re-5*k.^2).*R.^2).^2)));
% %                 
                a = (5*eps.*(10*eps.*R.^2.*Bo.^2+30*Bo.^2.*R.^3-eps.*k.^2.*(k.^2.*R.^2-1).*((5*k.^2-2*Bo.*Re).*R.^2-5)))./(225*Bo.^2.*R.^4+150*Bo.^2.*eps.*R.^3+eps.^2.*(25*R.^2.*Bo.^2+k.^2.*(5+(2*Bo.*Re-5*k.^2).*R.^2).^2));
                
                theta = atan(Bo.*k.*R.*(-15*k.^2.*R.^3+15*R+eps.*(15+(-15*k.^2+4*Re.*Bo).*R.^2))./(30*Bo.^2.*R.^3+eps.*(10*Bo.^2.*R.^2-k.^2.*(k.^2.*R.^2-1).*((5*k.^2-2.*Bo.*Re).*R.^2-5))))+(1-sign(a))/2*pi;
                A = 5.*eps.*sqrt((4*R.^2.*Bo.^2+(k-k.^3.*R.^2).^2)./(225*Bo.^2.*R.^4+150*eps.*Bo.^2.*R.^3+eps.^2.*(25*R.^2.*Bo.^2+k.^2.*(5 + (2*Bo.*Re-5*k.^2).*R.^2).^2)));
                
                hinit = obj.mass0 + obj.del*A*cos(k.*obj.z-theta);
            elseif obj.equation ==1
                theta = atan((3*obj.Bo*obj.L^3*obj.R^2)/(2*pi*(4*pi^2-obj.L^2)));
                if theta<0
                    theta = theta + 2*pi;
                end
                
                hinit = 1 + obj.del*-cos(theta)*cos(2*pi/obj.L*obj.z-theta);
            end
            

            
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
        function h = small_ep(obj)
            h = -1/3*(obj.eta/obj.R + 1/obj.Bo*(obj.etaz/obj.R^2+obj.etazzz));
         
        end
        
        function F = hfun(obj,h)
            
            if obj.flow == 1
                [hz,~,hzzz,~] = obj.getdiv(h);
                if obj.equation ==0
                    F = h.^3/3 + obj.ep*(2/15*obj.Re*h.^6.*hz+h.^3/(3*obj.Bo).*((hz+obj.etaz)/obj.R^2+ hzzz+ obj.etazzz)-h.^4/(6*obj.R)-2/(3*obj.R)*h.^3.*obj.eta) - obj.q;
                elseif obj.equation ==1
                    F = h.^3/3.*(1 + 1/obj.Bo*((hz+obj.etaz)/obj.R^2 + hzzz  + obj.etazzz))-obj.q;
                end
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
                if obj.equation == 0
                    m = -(obj.eta+obj.R./obj.ep)+sqrt((obj.eta+obj.R/obj.ep).^2 +2*obj.R*H/obj.ep);
                else
                    m=H;
                end
                F(l) = sum(m)/obj.n - 1;
                
                
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
                    yF = fft(a,[],2);
                    
                    N = obj.n/2;
                    
                    % Determine k in matlab form
                    %k = fftshift([0,-N+1:N-1])';
                    
                    
                    % Prior suppression
                    yF(abs(yF)<obj.suppression) = 0;
                    %yF = [yF,zeros(1,3*obj.n)]*4;
                    % Apply pseudo-spectral differentiation
                    padzeros = zeros(dim(1),(obj.pad-1)*obj.n);
                    
                    yF = [yF(:,1:N),padzeros, yF(:,N+1:obj.n)]*obj.pad;
                    k = [0:N-1,padzeros(1,:), 0, 1-N:-1] * 2*pi/obj.L;
                    % Posterior suppression
                    % dyF(abs(dyF) < suppression*N*2) = 0 ;
                    % dyF(abs(dyF) < suppression*max(abs(dyF))) = 0 ;
                    
                    % Transform back into real space
                    dyF = (1i*k).^1.*yF;
                    az = real(ifft(dyF,[],2));
                    %az = real(ifft(dyF,obj.n*obj.pad));
                    az = az(:,1:obj.pad:end);
                    dyF = (1i*k).^2.*yF;
                    azz = real(ifft(dyF,[],2));
                    %azz = real(ifft(dyF,obj.n*obj.pad));
                    azz = azz(:,1:obj.pad:end);
                    dyF = (1i*k).^3.*yF;
                    
                    %azzz = real(ifft(dyF,obj.n*obj.pad));
                    azzz = real(ifft(dyF,[],2));
                    azzz = azzz(:,1:obj.pad:end);
                    dyF = (1i*k).^4.*yF;
                    azzzz = real(ifft(dyF,[],2));
                    %azzzz = real(ifft(dyF,obj.n*obj.pad));
                    azzzz = azzzz(:,1:obj.pad:end);
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
            
            
        
        
        function obj = get_h(obj,optionon,hinit,q)
            %performs fsolve
            %options = optimoptions('fsolve','Display', 'none','FiniteDifferenceType','central');
            options = optimoptions('fsolve','Display', 'none','FiniteDifferenceType','central','FunctionTolerance',1e-10,'OptimalityTolerance',1e-14,'StepTolerance',1e-10);
            if nargin<4
                q = 1/3*(1+obj.ep*(obj.equation-1));
            if nargin <=2
                hinit = obj.linear_shape;
                %hinit = ones(1,obj.n);
                
            end
            end
            if nargin >1
                if optionon ==1
                    options = optimoptions('fsolve','FiniteDifferenceType','central','FunctionTolerance',1e-10,'OptimalityTolerance',1e-6,'StepTolerance',1e-10,'Display','Iter');
                %'Display', 'iter','PlotFcn',@optimplotfirstorderopt
                end
            end
            if obj.flow == 0
                %always starting with 1/3
                hinit = [hinit q];
            end
            [obj.h0,obj.Fval,obj.eflag,obj.out,obj.jac] = fsolve(@obj.hfun,hinit,options);
            if obj.eflag <= 0
                warning('no steady state solution found')
            end
            if obj.flow == 0
                obj.q = obj.h0(end);
                obj.h0 = obj.h0(1:end-1);
            end
            obj.mass0 = sum(obj.h0)/obj.n;
            
        end
        
        function obj = new_get_h(obj,optionon,hinit)
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
        function err = mass_err(obj)
            err = mean(abs(obj.mass-1));
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
                Fh(obj.n/2+1) = 0;
                Fh(abs(Fh)<obj.suppression) = 0;
                init = ifft(Fh);
                
            end
            if obj.usejac == 1
                %opts = odeset('Stats','on','Jacobian',@obj.get_Jac);
                sparse = diag(ones(obj.n,1))+ diag(ones(obj.n -1,1),1)+ diag(ones(obj.n -2,1),2)+diag(ones(obj.n-1,1),-1)+diag(ones(obj.n-2,1),-2) + diag(1,1-obj.n) + diag([1,1],2-obj.n) + diag(1,obj.n-1) + diag([1,1],obj.n-2);
                opts = odeset('RelTol',obj.reltol,'AbsTol',obj.abstol,'Stats','on','Jacobian',@obj.get_Jac,'JPattern',sparse);
            else
                
                opts = odeset('RelTol',obj.reltol,'AbsTol',obj.abstol,'Stats','on','BDF','on');
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
                [obj.t,obj.h]= ode15s(@(t,h) obj.odefun(t,h) ,0:obj.delt:obj.T,init,opts);
                %[obj.t,obj.h]= ode15s(@obj.odefun ,[0, obj.T],init,opts);
                %[obj.t,obj.h]= ode23s(@obj.odefun ,[0,obj.T],init,opts);
            else
                %obj.sol = ode15s(@obj.odefun ,[0 obj.T],init);
                %[obj.t, obj.h] = ode15s(@obj.odefun ,0:0.05:obj.T,init);
                [obj.t, obj.h] = ode15s(@(t,h) obj.odefun(t,h) ,[0,obj.T],init);
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
%             if obj.force_mass == 1
%                 mass = h(end);
%                 h = h(1:end-1);
%             end
            if obj.filter_oscillation == 1
                Fh = fft(h);
                Fh(obj.n/2+1) = 0;
                Fh(abs(Fh)<obj.suppression) = 0;
                h = ifft(Fh);
            end
            
            
                

            
            [hz,hzz,hzzz,hzzzz ] = obj.getdiv(h);

            if obj.equation == 1
                %small Bond
                ht = -hz.*h.^2 -(h.^3/(3*obj.Bo).*((hzz+obj.etazz)/obj.R^2+hzzzz+obj.etazzzz)+hz.*h.^2/obj.Bo.*((hz+obj.etaz)/obj.R^2+hzzz+obj.etazzz));
            elseif obj.equation == 2
                %
                ht = -3*hz.*h.^2 -(h.^3.*((hzz)*obj.R+hzzzz)+3*hz.*h.^2.*(hz*obj.R+hzzz));
            elseif obj.equation ==4
                ht = -h.*hz +obj.R*hzz;
            elseif obj.equation == 5
                ht = -hzz-obj.R*hzzzz-h.*hz;
            elseif obj.equation == 6
                ht = -3*h.^2.*hz.^2 - h.^3.*hzz - hz.*hzzz-h.*hzzzz;
            elseif obj.equation == 10
                %a equation
                ht = -hz.*h.^2 -obj.ep*(2/15*obj.Re*(h.^6.*hzz+6*h.^5.*hz.^2)+h.^3/(3*obj.Bo).*((hzz+obj.etazz)/obj.R^2+hzzzz+obj.etazzzz)+hz.*h.^2/obj.Bo.*((hz+obj.etaz)/obj.R^2+hzzz+obj.etazzz)+h.^3.*hz/(3*obj.R) + h.^3.*obj.etaz/(3*obj.R));
            elseif obj.equation == 20
                
                ht = -0.1*hz.*h.^2;
            else
                
                %main equation
                ht = -hz.*h.^2 -obj.ep*(2/15*obj.Re*(h.^6.*hzz+6*h.^5.*hz.^2)+h.^3/(3*obj.Bo).*((hzz+obj.etazz)/obj.R^2+hzzzz+obj.etazzzz)+hz.*h.^2/obj.Bo.*((hz+obj.etaz)/obj.R^2+hzzz+obj.etazzz)-2*h.^3.*hz/(3*obj.R)-2*hz.*h.^2.*obj.eta/obj.R - 2*h.^3.*obj.etaz/(3*obj.R));
                
            end
            if obj.force_mass ==1
                Fht = fft(ht);
                Fht(1) = 0;
                ht = real(ifft(Fht));
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
        function M = animate_new(obj,t0,tint,tend,filename,wall,periods,c)
            if nargin<8
                c = 0;
                if nargin<7
                    periods = 1;
                    if nargin<6
                        wall = 0;
                        if nargin<5
                            filename = 0;
                        if nargin<4
                            tend = obj.t(end);
                            if nargin<3
                                tint = obj.delt;
                                if nargin<2
                                    t0 = obj.t(1);
                                end
                            end
                        end
                        end
                    end
                end
            end
            record = ischar(filename);
            if record
                v = VideoWriter(strcat('../video/',filename),'MPEG-4');
                open(v);
            end
            i0 = floor((t0-obj.t(1))/obj.delt)+1;
            int = floor(tint/obj.delt);
            if int<1
                int =1;
            end
            iend = floor((tend-t0)/obj.delt)+i0;
                
          	ivec = i0:int:iend;
            
            clf
            
            ax = gca;
            ax.NextPlot = 'replaceChildren';
            obj.z = 0:obj.L/obj.n:periods*obj.L-obj.L/obj.n;
            obj = obj.get_wall;
            obj.h = repmat(obj.h,1,periods);
            
            if wall == 1
                ymin = min(obj.eta);
            else
                ymin = min(min(obj.a(ivec,:),[],2));
            end
            ylim([ymin-0.05,max(max(obj.a(ivec,:)+wall*obj.eta,[],2))]+0.05)
            xlim([0,periods*obj.L])
            loops = length(ivec);
            M(loops) = struct('cdata',[],'colormap',[]);
            
            for i = 1:loops
                if wall == 1
                    plot(obj.z,obj.eta)
                    hold on;
                end
                plot(obj.z,obj.a(ivec(i),:)+wall*obj.eta)
                            ylim([ymin-0.05,max(max(obj.a(ivec,:)+wall*obj.eta,[],2))]+0.05)
            xlim([0,periods*obj.L])
                title(sprintf('$t = %g$',obj.t(ivec(i))))
                drawnow
                
   
                M(i) = getframe(gcf);
                if record
                writeVideo(v,M(i));
                end
                hold off
            end
            if record
            close(v);
            end
            
        end
        
        function animate(obj,tint,wall,c,clearf)
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
                tint = 1;
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
            for i = 1:tint:l
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
            figure(10+m),clf,  hold on, title("Amplitude Changes with time")
            figure(11+m),clf,  hold on, title('Position of maxima/minima')
            figure(12+m),clf,  hold on , title('Drift of results')
            d = obj.del;
            
            figure(10+m)
            
            plot(obj.t,obj.diff, 'DisplayName',sprintf('$\\delta = %g$',d))
            figure(11+m)
            plot(obj.t,obj.hmaxloc/obj.L, 'DisplayName',sprintf('maxima $ \\delta = %g$',d))
            fig = gca;
            hold on
            plot(obj.t,obj.hminloc/obj.L,'--','Color', fig.Children(1).Color, 'DisplayName',sprintf('minima$ \\delta = %g$',d))
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
            del = 0.05:0.05:4;
            for i = 1:length(del)
                obj.del = del(i);
                obj = obj.get_h;
                if obj.eflag<= 0 | obj.q<0
                    break
                end
                obj.h = obj.h0;
                h2(i) = obj.h2norm;
            end
            l = length(h2);
            plot(del(1:l),h2)
        end
              function qvaryamp(obj)
            del = 0.05:0.05:4;
            for i = 1:length(del)
                obj.del = del(i);
                obj = obj.get_h;
                if obj.eflag<= 0 | obj.q<0
                    break
                end
                obj.h = obj.h0;
                q(i) = obj.q;
            end
            l = length(q);
            plot(del(1:l),q)
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
            xlabel('$t$')
            %ylabel('$||h||_2$')
            ylabel('$||H||_2$')
            title(sprintf('$||h||_2$ for $\\delta =%g$, $L= %g\\pi$',obj.del,obj.L/pi))
        end

        function plotz0(obj,m)
            % m = [0 1 2 3 10] plots the position at the maixima, 0 minima and 0
            if nargin ==1
                m = 0;
            end
            
            if m == 10
                plot(obj.t,mean(obj.a,2))
            else
                plot(obj.t, obj.a(:,obj.n/4*m+1))
            end
            xlabel('$t$')
            ylabel('$h$')
            maxstring = {'maximum', '0 after max','minimum','0 after min'};
            title(sprintf('Fluid thickness over the wall %s for $\\delta =%g$, $L= %g\\pi$',maxstring{m+1},obj.del,obj.L/pi))
        end
        
        function plotQ(obj)
            plot(obj.t,obj.Qint)
            xlabel('$t$')
            ylabel('$Q$')
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
        function phaseplotE(obj,t,earlytraj)
            if nargin <3
                earlytraj = 0;
            end
            if nargin == 1
                t = obj.T/2;
            end
            nt = floor(t/obj.delt);
                if earlytraj == 1
                    hold on 
                    plot(obj.h2norm(1:nt),obj.h2norm(2:nt+1),'--','Color',[0.8 0.8 0.8])
                end
                plot(obj.h2norm(nt:end-1),obj.h2norm(nt+1:end),'Color',[0 0.4470 0.7410])
                
                hold off
                xlabel('$E(t)$')
                ylabel('$E(t+\delta t$')
                title(sprintf('Energy compared between consecutive time steps for $L = %g\\pi$, $\\delta = %g$',obj.L/pi,obj.del))
            
        end
        function phaseplotEt(obj,t,earlytraj)
            if nargin <3
                earlytraj = 0;
            end
            if nargin == 1
                t = obj.t(1)+obj.T/2;
            end
            nt = floor((t-obj.t(1))/obj.delt +1);
                if earlytraj == 1
                    hold on 
                    plot(obj.h2norm(1:nt),(obj.h2norm(2:nt+1)-obj.h2norm(1:nt))/obj.delt,'--','Color',[0.8 0.8 0.8])
                end
                plot(obj.h2norm(nt:end-1),(obj.h2norm(nt+1:end)-obj.h2norm(nt:end-1))/obj.delt,'Color',[0 0.4470 0.7410])
                
                hold off
                xlabel('$||H(t)||_2$')
                ylabel('$\frac{d||H||_2}{dt}$')
                title(sprintf('Phase plot for $L = %g\\pi$, $\\delta = %g$',obj.L/pi,obj.del))
            
        end
            
            
        function phaseplotEQ(obj,t,earlytraj)
            if nargin <3
                earlytraj = 0;
            end
            if nargin == 1
                t = obj.t(1)+obj.T/2;
            end
            nt = floor((t-obj.t(1))/obj.delt +1);
                if earlytraj == 1
                    hold on 
                    plot(obj.h2norm(1:nt),obj.Qint(1:nt),'--','Color',[0.8 0.8 0.8])
                end
                plot(obj.h2norm(nt:end),obj.Qint(nt:end),'Color',[0 0.4470 0.7410])
                

                xlabel('$||h(t)||_2$')
                ylabel('$Q$')
                title(sprintf('Phase plot for $L = %g\\pi$, $\\delta = %g$',obj.L/pi,obj.del))
            
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
            
            for i = 1:length(obj.h)
            [pks,loc]  = findpeaks(obj.a(i,:));
            obj.pks{i} = pks;
            obj.loc{i} = loc;
            end
            
        end
        function plot_peaks(obj,t0)
            obj = obj.get_peak_data;
            if nargin ==1
                t0i = 1;
            else
            t0i = floor(t0/obj.delt);
            end
            hold on 
            for i = t0i:length(obj.h)
                plot(obj.z(obj.loc{i}),obj.pks{i},'.','color','#0072BD')
            end
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
        
        function obj = follow_peak(obj,n,t,zreset,region)
            if nargin<3
                t = obj.t(1);
                zreset = 0;
                region = 10;
            end
            if nargin == 1
                
                n = 16;
                
            end
            obj2 = obj;
            obj2.n = obj.n*n;
            obj2.ps = 0;
            [ba,i2] = min(abs(obj.t -t));
            l = length(obj.t(1:end));
            
            m = 0 ;
            a = interp1([obj.z,obj.L],[obj.a(i2,:),obj.a(i2,1)],0:obj.L/obj.n/n:obj.L - obj.L/obj.n/n,'spline','extrap');
            [h, k ] = max(a);
            [~,hzz,~,~ ] = obj2.getdiv(a);
            z0 = (k-1)/n/obj.n*obj.L;
            obj.zpos(1) = z0;
            obj.peakh(1) = h;
            
            obj.peakp(1)  = -1/obj2.Bo*((obj2.eta(k)+h)/obj2.R^2+obj2.etazz(k)+ hzz(k));
            for i = 1:l-i2
                a = interp1([obj.z,obj.L],[obj.a(i+i2,:),obj.a(i+i2,1)],0:obj.L/obj.n/n:obj.L - obj.L/obj.n/n,'spline','extrap');
                [~,hzz,~,~ ] = obj2.getdiv(a);
                if k+region<=n*obj.n && k-region>0
                    [h, j ] = max(a(k-region:k+region));
                    k = k-region + j - 1;
                    
                elseif k+region>obj.n
                    [ma, j ] = max(a(k-region:obj.n*n));
                    
                    [mb,j2] = max(a(1:k+region-n*obj.n));
                    if ma>mb
                        k = k-region+j-1;
                        h = ma;
                    else
                        k = j2;
                        m= m+1;
                        h = mb;
                    end
                    
                    
                elseif k-region<1
                    [ma, j ] = max(a(k-region+n*obj.n:n*obj.n));
                    [mb,j2] = max(a(1:k+region));
                    if ma > mb
                        m= m-1;
                        k = k-region+obj.n+j-1;
                        h = ma;
                    else
                        k = j2;
                        h= mb;
                    end
                    
                end
                
                obj.peakh(i+1) = h;
                obj.zpos(i+1) = (k-1)/obj.n/n*obj.L+m*obj.L;
                obj.peakp(i+1)  = -1/obj2.Bo*((obj2.eta(k)+h)/obj2.R^2+obj2.etazz(k)+ hzz(k));
                %                 if k ==1
                %                     plot(obj.t(i+i2),obj.zpos(i+1),'kx')
                %                 elseif k== 251
                %                      plot(obj.t(i+i2),obj.zpos(i+1),'ko')
                %                 end
                
                
            end
            %obj.zpos = smooth(smooth(obj.zpos,101),101);
            obj.speed = [(obj.zpos(2)-obj.zpos(1))/obj.delt,(obj.zpos(3:end)-obj.zpos(1:end-2))/(2*obj.delt),(obj.zpos(end)-obj.zpos(end-1))/obj.delt];
            
            %plot(obj.t(i2:end),obj.zpos)
            
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
        function surfdata(obj,c,shift)
            if nargin <=2
                shift = 0;
            end
            if nargin ==1
            obj = obj.follow_peak;
            c = polyfit(obj.t(floor(0.5*end):end),obj.zpos(floor(0.5*end):end),1);
            end
            phi = mod((obj.z-c(1)*obj.t+shift),obj.L);
            [phi,I] = sort(phi,2);
            I = I + obj.n*(0:length(I)-1)';
%             for j = 1:length(obj.t)/100
%                 a(100*j,:) = obj.a(100*j,I(j,:));
%             end
            
            
            a = obj.a';
            
            pcolor(phi,obj.t,a(I))
            %pcolor(obj.t,phi',a(I)')
          
             
          

            shading interp
            
            colorbar
            xlabel('$z-ct$')
            ylabel('$t$')
            xlim([0,obj.L])
            del = erase(string(obj.del),'.');
            if obj.del == 1
                title(sprintf('Fluid Thickness for wall $\\eta =\\cos %.2g z$, $c = %g$',2*pi/obj.L ,c(1)))
            else
                
            title(sprintf('Fluid Thickness for wall $\\eta =%g\\cos %.2g z$, $c = %g$',obj.del,2*pi/obj.L ,c(1)))
            end
            name = sprintf('../plots/colour/ColourL%gdel%s',obj.L/pi,del);
%             savefig(name)
%             saveas(gcf,name,'epsc')
%             
        end
        
        
        function obj = eliminate_noise(obj)
            
        Fh = fft(obj.h,[],2);
        Fh(:,obj.n/2+1) = 0;
        
        Fh(abs(Fh)<obj.suppression) = 0;
        %yF = [yF,zeros(1,3*obj.n)]*4;
        % Apply pseudo-spectral differentiation
        
        
        
        obj.h = real(ifft(Fh,[],2));
        
        end
        function plot_h2peaks(obj)
            [~,hloc] = findpeaks(obj.h2norm);
            hloc = hloc(hloc>2000);
            hold on 
            for i = hloc
                plot(obj.z,obj.h(i,:))
            end
        end
        function [V,d] = Floquet(obj,h)
            if nargin ==1
                obj = obj.get_h;
            else
                obj.h0 = h;
            end
            [h0z,h0zz,h0zzz,h0zzzz] = obj.getdiv(obj.h0);
            h0 = obj.h0;
            eta = obj.eta;
            etaz = obj.etaz;
            etazz = obj.etazz;
            etazzz = obj.etazzz;
            etazzzz = obj.etazzzz;
            if obj.equation == 1
            H = 2*h0.*h0z.*(1 + 1/obj.Bo*((h0z+etaz)/obj.R^2 + h0zzz+etazzz))+h0.^2/obj.Bo.*((h0zz+etazz)/obj.R^2 +h0zzzz + etazzzz);
            Hz = h0.^2.*(1 + 1/obj.Bo*((2*h0z+etaz)/obj.R^2 + h0zzz+etazzz));
            Hzz = h0.^3/(3*obj.Bo*obj.R^2);
            Hzzz = (h0.^2.*h0z)/(3*obj.Bo);
            Hzzzz = h0.^3/3/obj.Bo;
            elseif obj.equation ==0
            H = 2*h0.*h0z.*(1 + obj.ep/obj.Bo*((h0z+etaz)/obj.R^2 + h0zzz+etazzz))+h0.^2*obj.ep/obj.Bo.*((h0zz+etazz)/obj.R^2 +h0zzzz + etazzzz)+obj.ep*(2/15*obj.Re*30*h0.^4.*h0z.^2 -2*h0.^2.*h0z/obj.R-4*h0.*h0z.*eta/obj.R-2*h0.^2.*etaz/obj.R);
            Hz = h0.^2.*(1 + obj.ep/obj.Bo*((2*h0z+etaz)/obj.R^2 + h0zzz+etazzz))+obj.ep*2/15*obj.Re*12*h0.^5.*h0z-obj.ep*2/3/obj.R*h0.^3 -2*obj.ep/obj.R*h0.^2.*eta;
            Hzz = h0.^3*obj.ep/(3*obj.Bo*obj.R^2)+obj.ep*2/15*obj.Re*h0.^6;
            Hzzz = (h0.^2.*h0z)*obj.ep/(3*obj.Bo);
            Hzzzz = h0.^3/3*obj.ep/obj.Bo;
            end
            
            d0 = H - 2*(obj.n/obj.L)^2*Hzz+6*(obj.n/obj.L)^4*Hzzzz;
            dn1 = -1/2*(obj.n/obj.L)*Hz + (obj.n/obj.L)^2*Hzz + (obj.n/obj.L)^3*Hzzz- 4*(obj.n/obj.L)^4*Hzzzz;
            dn2 = -1/2*(obj.n/obj.L)^3*Hzzz + (obj.n/obj.L)^4*Hzzzz;
            d1 = 1/2*(obj.n/obj.L)*Hz+ (obj.n/obj.L)^2*Hzz-(obj.n/obj.L)^3*Hzzz - 4*(obj.n/obj.L)^4*Hzzzz;
            d2 = 1/2*(obj.n/obj.L)^3*Hzzz+(obj.n/obj.L)^4*Hzzzz;
            Jac = diag(d0)+ diag(d1(1:end-1),1)+ diag(d2(1:end-2),2)+diag(dn1(2:end),-1)+diag(dn2(3:end),-2) + diag(d1(end),1-obj.n) + diag(d2(end-1:end),2-obj.n) + diag(dn1(1),obj.n-1) + diag(dn2(1:2),obj.n-2);
            [V,D] = eig(Jac);
            d = diag(D);
            index = d<0;
            d = d(index);
            V = real(V(:,index));
            
            
        end
            
        function n = most_unstable(obj)
            [V,d] = obj.Floquet();
            [val,loc] = max(-real(d));
            [pval,ploc] = findpeaks(V(:,loc));
            n = length(pval);
        end
        
        function n = num_unstable(obj,h)
            if nargin ==2
            [V,d] = obj.Floquet(h);
            else
                [V,d] = obj.Floquet();
            end
            
           %
           n = floor(length(d)/2);
        end
        function c = init_speed(obj)
            [~,d] = obj.Floquet();
            c = -1i*d*obj.L/2/pi;
            c = abs(real(c));
        end
        function [k,kflat] = init_growth(obj)
           [~,d] = obj.Floquet();
            
            k = real(d);
            kflat = (imag(d).^4-imag(d).^2)/3/obj.Bo;
        end
        function c = get_c_spectrally(obj,delay)
            if nargin ==1
                delay = 200;
            end
           
           c = [loc(loc2)-0.5,loc(loc2), loc(loc2)+0.5]/(obj.T-delay)*obj.L;
        end
        
        function plot_period(obj,t,delt,n)
            if nargin<4
                n = 1; 
            end
            if nargin<3
                delt = 0.1;
            end
            if nargin <2
                t = 300;
            end
            
            [hmax,hloc0] = findpeaks(obj.h2norm);
            [hmin,hmloc] = findpeaks(-obj.h2norm);
            %[hmax2,hloc2] = findpeaks(obj.Qint);
            %[hmin2,hmloc2] = findpeaks(-obj.Qint);
            hloc = hloc0(hloc0>(t-obj.t(1))/obj.delt);
            hmax1 = hmax(hloc0>(t-obj.t(1))/obj.delt);
            hm1 = hloc(1);
            hm2 = hloc(1+n);
            hl = hmloc(hmloc<hm2&hmloc>hm1);
            %h2m = hloc2(hloc2<hm2&hloc2>hm1);
            %h2l = hmloc2(hmloc2<hm2&hmloc2>hm1);
            obj = obj.get_h;
            hold on,
            title(sprintf('$h$ from between $t =%g$ and $t=%g$ at $t = %g$ intervals for wall $\\eta = %g\\cos %.2gz$',obj.t(hm1),obj.t(hm2),delt,obj.del,2*pi/obj.L))
            
            for i = hm1:floor(delt/obj.delt):hm2-1
                plot(obj.nz,obj.a(i,:),'color',[0.8,0.8,0.8],'linewidth',0.5,'HandleVisibility','off')
            end
            set(gca,'ColorOrderIndex',1)
            plot(obj.nz,obj.a(hm1,:),'linewidth',2,'DisplayName',sprintf('maximum $||h||_2 = %.4g$',hmax1(1)))
            plot(obj.nz,obj.a(hl,:),'linewidth',2,'DisplayName',sprintf('minimum $||h||_2 = %.4g$',-hmin(hmloc<hm2&hmloc>hm1)))
            %plot(obj.nz,obj.h(h2m,:),'linewidth',2,'DisplayName',sprintf('maximum $Q$ = %.4g',hmax2(hloc2<hm2&hloc2>hm1)))
            %plot(obj.nz,obj.h(h2l,:),'linewidth',2,'DisplayName',sprintf('minimum $Q$ = %.4g',-hmin2(hmloc2<hm2&hmloc2>hm1)))
            %plot(value.nz,value.h(floor((hm1+hl)/2),:),'linewidth',2,'HandleVisibility','off')
            %plot(value.nz,value.h(floor((hm2+hl)/2),:),'linewidth',2,'HandleVisibility','off')
%             plot(obj.nz,sum(obj.h(hloc(1):hloc(end)-1,:))/(hloc(end)-1-hloc(1)),'linewidth',2,'DisplayName','time averaged thickness')
%             plot(obj.nz,obj.h0,'linewidth',2,'DisplayName','steady state thickness')
            legend()
            
            grid
            xticks([0 0.25 0.5 0.75 1])
            xlabel('$\frac{z}{L}$')
            ylabel('$h$')
        end
         function plot_all(obj)
            
            for i = hm1:hm2-1
                plot(obj.nz,obj.h(i,:),'color',[0.8,0.8,0.8],'linewidth',0.5,'HandleVisibility','off')
            end
            set(gca,'ColorOrderIndex',1)
            plot(obj.nz,obj.h(hm1,:),'linewidth',2,'DisplayName',sprintf('maximum $||h||_2 = %.4g$',hmax1(1)))
            plot(obj.nz,obj.h(hl,:),'linewidth',2,'DisplayName',sprintf('minimum $||h||_2 = %.4g$',-hmin(hmloc<hm2&hmloc>hm1)))
            plot(obj.nz,obj.h(h2m,:),'linewidth',2,'DisplayName',sprintf('maximum $Q$ = %.4g',hmax2(hloc2<hm2&hloc2>hm1)))
            plot(obj.nz,obj.h(h2l,:),'linewidth',2,'DisplayName',sprintf('minimum $Q$ = %.4g',-hmin2(hmloc2<hm2&hmloc2>hm1)))
            %plot(value.nz,value.h(floor((hm1+hl)/2),:),'linewidth',2,'HandleVisibility','off')
            %plot(value.nz,value.h(floor((hm2+hl)/2),:),'linewidth',2,'HandleVisibility','off')
%             plot(obj.nz,sum(obj.h(hloc(1):hloc(end)-1,:))/(hloc(end)-1-hloc(1)),'linewidth',2,'DisplayName','time averaged thickness')
%             plot(obj.nz,obj.h0,'linewidth',2,'DisplayName','steady state thickness')
            legend('NumColumns',2)
            
            grid
            xticks([0 0.25 0.5 0.75 1])
            xlabel('$\frac{z}{L}$')
            ylabel('$h$')
         end
         function [npks,time_periodic,c,T] = get_peak_info(obj,t,display)
             if nargin<3
                 display = 0;
                 if nargin<2
                 t = obj.t(floor(end/2));
                 end
             end
             n0 = floor((t-obj.t(1))/obj.delt)+1;
             hp = obj.h(n0:end,1);
             lim =  (min(max( obj.h(floor((t-obj.t(1))/obj.delt)+1:end,:)))+1)/2;
             ht = obj.h(floor((t-obj.t(1))/obj.delt)+1:end,1+obj.n/2);
             [pkp,locp] = findpeaks(hp);
             locp = locp(pkp>lim)+n0-1;
             pkp = pkp(pkp>lim);
             [pkt,loct] = findpeaks(ht);
             %loct = loct(pkt>lim)+n0-1;
             %pkt = pkt(pkt>lim);
             

             if length(locp)<4
                 c =0;
                 T = 0;
                 npks = 0;
                 time_periodic = -1;
                 return
             end
             if display ==1
                 figure
                 disp(pkp)
                plot(obj.h(floor((locp(end-1)+locp(end))/2),:))
             end
             [npk,~] = findpeaks(obj.h(floor((locp(end-1)+locp(end))/2),:));
             npks = sum(npk>lim);
             
             cp = [];
             ct = [];
             time_periodic = [];
             T = [];
             for i=1:npks
                 t_pks = diff(locp(i:npks:end));
                 if display == 1
                 figure
                 plot(obj.z,obj.h(locp(i:npks:end),:))
                 
                 figure
                 plot(t_pks)
                 end
                 tp_test = sum(abs(t_pks - mean(t_pks))>1);
                 if tp_test == 0
                     time_periodic = [time_periodic ,1];
                     cp = [cp, obj.L/(mean(t_pks)*obj.delt)];
                     T = [T, mean(t_pks)*obj.delt];
                 elseif tp_test<length(t_pks)/2
                     tp_test2 = sum(abs(t_pks(floor(end/2):end) - mean(t_pks(floor(end/2):end)))>1);
                     if tp_test2 == 0 
                         time_periodic = [time_periodic ,0.9];
                         cp = [cp, obj.L/(mean(t_pks(floor(end/2):end))*obj.delt)];
                         T = [T, mean(t_pks)*obj.delt];
                     else
                         time_periodic = [time_periodic ,0.5];
                         cp = [cp,0];
                         T = [T, mean(t_pks)*obj.delt];
                     end
                 else
                     time_periodic = [time_periodic ,0];
                     cp = [cp,0];
                     T = [T, mean(t_pks)*obj.delt];
                 end
               
             end
             c = mean(cp);
             time_periodic = mean(time_periodic);
             T = mean(T);
         end
         function hm = mean_h(obj,fit_periods)
             if nargin<2
                 fit_periods = 0;
             end
             l = length(obj.t);
             if fit_periods == 1
             [pks,loc] = findpeaks(obj.h2norm);
             if length(loc)>1
             loc = loc(abs(pks - pks(end))<1e-4);
             if length(loc)>1
                 i1 = loc(1);
                 iend = loc(end)-1;
             else
                 i1 = floor(l/2);
                 iend = l;
             end
             else
                 i1 = floor(l/2);
                 iend = l;
             end
             else
                 i1 = floor(l/2);
                 iend = l;
             end
             hm = sum(obj.h(i1:iend,:))/(obj.t(iend)-obj.t(i1))*obj.delt;
             
         end
         function [dif,Afac,pshift,tshift,pdiff,tdiff] = steady_mean_diff(obj)
             znew = 0:obj.L/obj.n/4:obj.L - obj.L/obj.n/4;
             mh = obj.mean_h(1);
             if length(obj.h0)~=obj.n
                obj = obj.get_h;
             end
             h0 = obj.h0;
             dif = norm(mh-h0,1)/obj.n;
             mh = interp1(0:obj.L/obj.n:obj.L,[mh mh(1)],0:obj.L/obj.n/4:obj.L - obj.L/obj.n/4 );
             h0 = interp1(0:obj.L/obj.n:obj.L,[h0 h0(1)],0:obj.L/obj.n/4:obj.L - obj.L/obj.n/4 );
             [mpk,mpl] = max(mh);
             [hpk,hpl] = max(h0);
             [mtg,mtl] = max(-mh);
             [htg,htl] = max(-h0);
             mtg = -mtg;
             htg = -htg;
             pshift = (mpl-hpl)/obj.n/4;
             tshift = (mtl-htl)/obj.n/4;
             Afac = (mpk - mtg)/(hpk - htg);
             pdiff = (mpk - hpk);
             tdiff = (mtg - htg);
         end
         function obj = CO2_growth(obj,R0,del,Q,L,gamma)
             if nargin<4
             gamma = 1e11;
             Q = 1e-9;
             L = 5e-3;
             if nargin<3
                 del = 0.1;
                 if nargin<2
                     R0 = 0.01;
                 end
             end
             end
             
             tic
             obj.del = del;
              f0 = ones(1,obj.n)*1e-7;
              c0 = 1+1/2*f0.*(obj.r.^2)*1e-3 +obj.ep*sin(obj.nz*2*pi);
             [c0,f0] = obj.get_conc2(c0,f0,R0,Q,L);
             obj = obj.getndparams(R0,Q,L);
             eta0 = obj.eta*obj.H;
             y0 = [eta0 R0];
             
             %opts = odeset('RelTol',obj.reltol,'AbsTol',obj.abstol,'Stats','on','BDF','on','Events',@obj.EventsFcn);
             opts = odeset('RelTol',obj.reltol,'AbsTol',obj.abstol,'Stats','on','BDF','on');
           [t,y] = ode15s(@(t,y) obj.growth_CO2_fun(t,y,Q,L,gamma,c0,f0),0:obj.delt:obj.T,y0 ,opts);
           obj.etat = 1;
           obj.t = t;
           obj.h = y;
           obj.integration_time = toc;
         end
         function dy = growth_CO2_fun(obj,t,y,Q,L,gamma,c0,f0)
            eta = y(1:obj.n)';
            R = y(obj.n+1);
            
             obj = obj.getndparams2(R,Q,L);
             obj.eta = eta/obj.H;
             [~,f] = obj.get_conc2(c0,f0,R,Q,L);
             etat = gamma*obj.H*(f-mean(f));
             Rt = obj.H*mean(f)*gamma;
             dy = [etat Rt]';
             
             
             
         end
         function h = simple_growth_function(obj,t,eta)
             obj.eta = eta';
             obj = obj.get_h;
             if obj.equation == 0
                 obj.h = obj.h0;
                 h0 = obj.a;
             else
                 h0 = obj.h0;
             end

             h = obj.gamma*(h0'-1);
%              if obj.eflag <= 0
%               stop = 0;
%              else
%                  stop = 1;
%              end
%              save('breakfunction','stop')
%              
             end
         
         function obj = simple_growth(obj)
             tic
             %opts = odeset('RelTol',obj.reltol,'AbsTol',obj.abstol,'Stats','on','BDF','on','Events',@obj.EventsFcn);
             opts = odeset('RelTol',obj.reltol,'AbsTol',obj.abstol,'Stats','on','BDF','on');
           [t,y] = ode15s(@(t,y) obj.simple_growth_function(t,y),0:obj.delt:obj.T,obj.eta ,opts);
           obj.etat = 1;
           obj.t = t;
           obj.h = y;
           obj.integration_time = toc;
         end
                  function obj = simple_growth45(obj)
             tic
             opts = odeset('RelTol',obj.reltol,'AbsTol',obj.abstol,'Stats','on','BDF','on');
             
           [t,y] = ode45(@(t,y) obj.simple_growth_function(t,y),[0:obj.delt:obj.T],obj.eta );
           obj.etat = 1;
           obj.t = t;
           obj.h = y;
           obj.integration_time = toc;
         end
         
         
         function [value,isterminal,direction] = EventsFcn(obj,t,y)
             load('breakfunction','stop')
         value = stop;

         % The value that we want to be zero
         isterminal = 1;  % Halt integration 
         direction = 0;   % The zero can be approached from either direction
         end
         
         
         function obj = fast_growth(obj)
             tic
             obj = obj.get_h;
              opts = odeset('RelTol',obj.reltol,'AbsTol',obj.abstol,'Stats','on','BDF','on','Events',@obj.EventsFcn);
           [t,y] = ode15s(@(t,y) obj.fast_growth_function(t,y),0:obj.delt:obj.T,[obj.h0,obj.eta] ,opts);
           obj.etat = 1;
           obj.t = t;
           obj.h = y;
           obj.integration_time = toc;
           
         end
         
         function dt = fast_growth_function(obj,t,y)
             y = y';
             l = length(y);
             h = y(1:end/2);
             eta = y(l/2+1:end);
             obj.eta = eta;
             [hz,hzz,hzzz,hzzzz ] = obj.getdiv(h);

            if obj.equation == 1
                %small Bond
                ht = -hz.*h.^2 -(h.^3/(3*obj.Bo).*((hzz+obj.etazz)/obj.R^2+hzzzz+obj.etazzzz)+hz.*h.^2/obj.Bo.*((hz+obj.etaz)/obj.R^2+hzzz+obj.etazzz));
            elseif obj.equation == 2
                %
                ht = -3*hz.*h.^2 -(h.^3.*((hzz)*obj.R+hzzzz)+3*hz.*h.^2.*(hz*obj.R+hzzz));
            elseif obj.equation ==4
                ht = -h.*hz +obj.R*hzz;
            elseif obj.equation == 5
                ht = -hzz-obj.R*hzzzz-h.*hz;
            elseif obj.equation == 6
                ht = -3*h.^2.*hz.^2 - h.^3.*hzz - hz.*hzzz-h.*hzzzz;
            elseif obj.equation == 10
                %a equation
                ht = -hz.*h.^2 -obj.ep*(2/15*obj.Re*(h.^6.*hzz+6*h.^5.*hz.^2)+h.^3/(3*obj.Bo).*((hzz+obj.etazz)/obj.R^2+hzzzz+obj.etazzzz)+hz.*h.^2/obj.Bo.*((hz+obj.etaz)/obj.R^2+hzzz+obj.etazzz)+h.^3.*hz/(3*obj.R) + h.^3.*obj.etaz/(3*obj.R));
            elseif obj.equation == 20
                
                ht = -0.1*hz.*h.^2;
            else
                
                %main equation
                ht = -hz.*h.^2 -obj.ep*(2/15*obj.Re*(h.^6.*hzz+6*h.^5.*hz.^2)+h.^3/(3*obj.Bo).*((hzz+obj.etazz)/obj.R^2+hzzzz+obj.etazzzz)+hz.*h.^2/obj.Bo.*((hz+obj.etaz)/obj.R^2+hzzz+obj.etazzz)-2*h.^3.*hz/(3*obj.R)-2*hz.*h.^2.*obj.eta/obj.R - 2*h.^3.*obj.etaz/(3*obj.R));
                
            end
            dt = [1/obj.gamma*ht  h]';
              
         end
         function [hout,etaout,tout] = Euler_growth(obj,res)
             if nargin <2
                 res = 1;
             end
             eta0 = obj.eta;
             obj = obj.get_h;
             h0 = obj.h0;
             etaold = eta0;
             h = h0;
             q = obj.q;
             tout = 0:obj.delt:obj.T;
             t = 0:obj.delt/res:obj.T;
             etaout = [etaold];
             hout = [h];
             
             for i = 2:length(t)
                 etan = etaold +obj.delt/res*(h-1);
                 obj.eta = etan;
                 obj = obj.get_h(0,h,q);
                 h = obj.h0;
                 q = obj.q;
                 if mod(i-1,res) == 0
                 hout = [hout;h];
                 etaout = [etaout;etan];
                 etaold = etan;
                 end
                 
                 
                 
             end
             
         end
         
         function [hout,etaout,tout] = RK_growth(obj,res)
             tic
             if nargin <2
                 res = 1;
             end
             eta0 = obj.eta;
             obj = obj.get_h;
             h0 = obj.h0;
             etao = eta0;
             h = h0;
             q = obj.q;
             tout = 0:obj.delt:obj.T;
             t = 0:obj.delt/res:obj.T;
             etaout = [etao];
             hout = [h];
             
             for i = 2:length(t)
                 k1 = h-1;
                 obj.eta = etao+obj.delt/res*k1/2;
                 obj = obj.get_h(0,h,q);
                 k2 = obj.h0-1;
                 obj.eta = etao+obj.delt/res*k2/2;
                 obj = obj.get_h(0,h,q);
                 k3 = obj.h0-1;
                 obj.eta = etao+obj.delt/res*k3;
                 obj = obj.get_h(0,h,q);
                 k4 = obj.h0-1;
                 etan = etao +obj.delt/res/6*(k1+2*k2+2*k3+k4);
                 obj.eta = etan;
                 obj = obj.get_h(0,h,q);
                 h = obj.h0;
                 q = obj.q;
                 if mod(i-1,res) == 0
                 hout = [hout;h];
                 etaout = [etaout;etan];
                 etao = etan;
                 end
                 
                 
                 
             end
             toc
             
         end
         function plot_fourier_modes(obj)
             y = fft(obj.h(floor(end/2):end,:),[],2);
             y = abs(y(:,2:11));
             %clf
             hold on 
             errorbar(1:10,sum(y)/(obj.T-obj.t(floor(end/2)))*obj.delt,mean(y) - min(y),max(y)-mean(y),'_','MarkerSize',50)
             xlim([0.5 10.5])
%              plot(max(y),'_')
%              plot(min(y),'_')
title('Mean value of Fourier modes, with range')
xlabel('$k$')
ylabel('$|H(k)|$')
         end
         function [c1,c2,f] = get_conc(obj,c10,c20)
             obj = obj.get_h;
             obj.h = obj.h0;
             options = optimoptions('fsolve','Display','iter');
             cinit = [c10;c20];
             [c] = fsolve(@obj.cfun,cinit,options);
             c1 = c(1:obj.nr+1,:);
             c2 = c(obj.nr+2:end-1,:);
             f = -0.9e-9*5*obj.nr*(-3/2*c1(1,:) +2*c1(2,:)-1/2*c1(3,:))./(obj.H*obj.a);
         end
         function F = cfun(obj,c)
             h0 = 1e-4;
             D1 = 0.9e-9;
             D2 = 1.3e-9;
             g = 9.81;
             nu = 1e-6;
             Pe2 = h0^3*g/(nu*D2);
             Pe1 = 1.3/0.9*Pe2;
             C1 = 5;
             C2 = 2.7e-2;
             km = 2e-4;
             kp = 1e-2;
             rhoc  = 3.69e-5;
             w = obj.w(2:end-1,:);
             U = obj.U(2:end-1,:);
             c1 = c(1:obj.nr+1,:);
             c2 = c(obj.nr+2:end,:);
             %f = c(end,:);
             c11 = c1(3:end,:);
             c10 = c1(2:end-1,:);
             c1n1 = c1(1:end-2,:);
             c1r = obj.nr*(c11 -c1n1)/2 ;
             c1rr = obj.nr^2*(c1n1+c11-2*c10);
             [c1z,~,~,~] =obj.getdiv(c10); 
             
             c21 = c2(3:end,:);
             c20 = c2(2:end-1,:);
             c2n1 = c2(1:end-2,:);
             c2r = obj.nr*(c21 -c2n1)/2 ;
             c2rr = obj.nr^2*(c2n1+c21-2*c20);
             [c2z,~,~,~] =obj.getdiv(c20); 

             F1 = c1rr +obj.ep/obj.R*c1r.*obj.a -obj.ep*Pe1*(U.*c1r.*obj.a+w.*c1z.*obj.a.^2);
             F2 = c2rr +obj.ep/obj.R*c1r.*obj.a+ kp*h0^2/D2*(2*C1*km/(kp*C2)*c10 -c20).*obj.a.^2 - obj.ep*Pe2*(U.*c2r.*obj.a + w.*c2z.*obj.a.^2);
             %bc1 = D1*C1*obj.nr*(-3/2*c1(1,:) +2*c1(2,:)-1/2*c1(3,:)) + h0*obj.a.*f;
             bc1 = D1*C1*obj.nr*(-3/2*c1(1,:) +2*c1(2,:)-1/2*c1(3,:)) - D2*C2*obj.nr*(3/2*c2(end,:) -2*c2(end-1,:)+1/2*c2(end-2,:)).*(1 + obj.ep/obj.R.*obj.a)+ mean(D2*C2*obj.nr*(3/2*c2(end,:) -2*c2(end-1,:)+1/2*c2(end-2,:)).*(1 + obj.ep/obj.R.*obj.a));
             bc2 = obj.nr*(-3/2*c2(1,:) +2*c2(2,:)-1/2*c2(3,:));
             bc3 = obj.nr*(3/2*c1(end,:) -2*c1(end-1,:)+1/2*c1(end-2,:));
             %bc4 = D2*C2*obj.nr*(3/2*c2(end,:) -2*c2(end-1,:)+1/2*c2(end-2,:)) + h0*obj.a.*(1-obj.ep/obj.R).*f;
             bc5 = c2(end,:) -1;
             %F = [F1;F2;bc1;bc2;bc3;bc4;bc5];
             F = [F1;F2;bc1;bc2;bc3;bc5];
             
             
         end
         function F = cfunfd(obj,c)
             h0 = 9.9106e-06;
             D1 = 0.9e-9;
             D2 = 1.3e-9;
             g = 9.81;
             nu = 1e-6;
             Pe2 = h0^3*g/(nu*D2);
             Pe1 = 1.3/0.9*Pe2;
             C1 = 5;
             C2 = 2.7e-2;
             km = 2e-4;
             kp = 1e-2;
             w = obj.w(2:end-1,:);
             U = obj.U(2:end-1,:);
             c1 = c(1:obj.nr+1,:);
             c2 = c(obj.nr+2:end,:);
             %f = c(end,:);
             c11 = c1(3:end,:);
             c10 = c1(2:end-1,:);
             c1n1 = c1(1:end-2,:);
             c1r = obj.nr*(c11 -c1n1)/2 ;
             c1rr = obj.nr^2*(c1n1+c11-2*c10);
             %[c1z,~,~,~] =obj.getdiv(c10); 
             cz1 = [c1(:,2:end) c1(:,1)];
             czn1 = [c1(:,end) c1(:,1:end-1)];
             c1z = 1/2*obj.n/obj.L*(cz1-czn1);
             
             c21 = c2(3:end,:);
             c20 = c2(2:end-1,:);
             c2n1 = c2(1:end-2,:);
             c2r = obj.nr*(c21 -c2n1)/2 ;
             c2rr = obj.nr^2*(c2n1+c21-2*c20);
             %[c2z,~,~,~] =obj.getdiv(c20); 
                          cz2 = [c2(:,2:end) c2(:,1)];
             czn2 = [c2(:,end) c2(:,1:end-1)];
             c2z = 1/2*obj.n/obj.L*(cz2-czn2);
             
             c1z = c1z(2:end-1,:);
             c2z = c2z(2:end-1,:);
             F1 = c1rr +obj.ep/obj.R*c1r.*obj.a -obj.ep*Pe1*(U.*c1r.*obj.a+w.*c1z.*obj.a.^2);
             F2 = c2rr +obj.ep/obj.R*c2r.*obj.a+ kp*h0^2/D2*(2*C1*km/(kp*C2)*c10 -c20).*obj.a.^2 - obj.ep*Pe2*(U.*c2r.*obj.a + w.*c2z.*obj.a.^2);
             %bc1 = D1*C1*obj.nr*(-3/2*c1(1,:) +2*c1(2,:)-1/2*c1(3,:)) + h0*obj.a.*f;
             bc1 = D1*C1*obj.nr*(-3/2*c1(1,:) +2*c1(2,:)-1/2*c1(3,:)) - D2*C2*obj.nr*(3/2*c2(end,:) -2*c2(end-1,:)+1/2*c2(end-2,:))*(1 + obj.ep/obj.R);
             bc2 = obj.nr*(-3/2*c2(1,:) +2*c2(2,:)-1/2*c2(3,:));
             bc3 = obj.nr*(3/2*c1(end,:) -2*c1(end-1,:)+1/2*c1(end-2,:));
             %bc4 = D2*C2*obj.nr*(3/2*c2(end,:) -2*c2(end-1,:)+1/2*c2(end-2,:)) + h0*obj.a.*(1-obj.ep/obj.R).*f;
             bc5 = c2(end,:) -1;
             %F = [F1;F2;bc1;bc2;bc3;bc4;bc5];
             F = [F1;F2;bc1;bc2;bc3;bc5];
             
             
         end
         function obj = getndparams(obj,R1,Q,L)
L = L/(2*pi);
obj.H = (3*Q*1e-6./(2*pi*9.81*R1)).^(1/3);
obj.R = R1./L;
obj.ep = obj.H./L;
obj.Bo = 997*L.^2*9.81/0.072./(obj.equation*(obj.ep -1)+1);
obj.Re = 3*Q./(2*pi*R1*1e-6);
obj.L = 2*pi;
         end
function obj = getndparams2(obj,R1,Q,L)
%L = L/(2*pi);
obj.H = (3*Q*1e-6./(2*pi*9.81*R1)).^(1/3);
obj.R = 1;
obj.ep = obj.H./R1;
obj.Bo = 997*R1.^2*9.81/0.072./(obj.equation*(obj.ep -1)+1);
obj.Re = 3*Q./(2*pi*R1*1e-6);
obj.L = L/R1;
end
                  function [c2,f] = get_conc2(obj,c20,f0,R,Q,L)
                      if nargin<=3
                          R = 0.1;
                          Q = 1e-9;
                          L = 1e-2;
                      end
                      eta = obj.eta;
                      obj = obj.getndparams2(R,Q,L);
                      obj.eta = eta;
                                   obj = obj.get_h;
             obj.h = obj.h0;
             options = optimoptions('fsolve','Display','none');
             cinit = [c20;f0];
             [c,~,eflag] = fsolve(@obj.c2fun,cinit,options);
             if eflag<0
                 c= 0;
             end
             c2 = c(1:obj.nr+1,:);
             
             f = c(end,:);
                  end
         function [c1,c2,f] = get_conc3(obj,c10,c20,f0,R,Q,L)
                      if nargin<=4
                          R = 0.1;
                          Q = 1e-9;
                          L = 1e-2;
                      end
                      eta = obj.eta;
                      obj = obj.getndparams2(R,Q,L);
                      obj.eta = eta;
                                   obj = obj.get_h;
             obj.h = obj.h0;
             options = optimoptions('fsolve','Display','none','FiniteDifferenceType','central','FunctionTolerance',1e-8,'StepTolerance',1e-8);
             cinit = [c10;c20;f0];
             [c] = fsolve(@obj.c3fun,cinit,options);
             c1 = c(1:obj.nr+1,:);
             c2 = c(obj.nr+2:end-1,:)
             
             f = c(end,:);
         end
          function F = c2fun(obj,c)
             h0 = obj.H;
             D1 = 0.9e-9;
             D2 = 1.3e-9;
             g = 9.81;
             nu = 1e-6;
             Pe2 = h0^3*g/(nu*D2);
             Pe1 = 1.3/0.9*Pe2;
             C1 = 5;
             C2 = 2.7e-2;
             km = 2e-4;
             kp = 1e-2;
             rhoc  = 3.69e-5;
            
             
             w = obj.w(2:end-1,:);
             U = obj.U(2:end-1,:);
             c1 = ones(1,obj.n);
             c2 = c(1:obj.nr+1,:);
             %c2 = c(obj.nr+2:end-1,:);
             f = c(end,:);
             
             
             c21 = c2(3:end,:);
             c20 = c2(2:end-1,:);
             c2n1 = c2(1:end-2,:);
             c2r = obj.nr*(c21 -c2n1)/2 ;
             c2rr = obj.nr^2*(c2n1+c21-2*c20);
             [c2z,~,~,~] =obj.getdiv(c20); 
             
             F2 = c2rr +obj.ep/obj.R*c2r.*obj.a+ kp*h0^2/D2*(2*C1*km/(kp*C2)*c1 -c20).*obj.a.^2 - obj.ep*Pe2*(U.*c2r.*obj.a + w.*c2z.*obj.a.^2);
             %bc1 = D1*C1*obj.nr*(-3/2*c1(1,:) +2*c1(2,:)-1/2*c1(3,:)) + h0*obj.a.*f;
             %bc1 = D1*C1*obj.nr*(-3/2*c1(1,:) +2*c1(2,:)-1/2*c1(3,:)) - D2*C2*obj.nr*(3/2*c2(end,:) -2*c2(end-1,:)+1/2*c2(end-2,:))*(1 + obj.ep/obj.R);
             bc2 = obj.nr*(-3/2*c2(1,:) +2*c2(2,:)-1/2*c2(3,:));
             %bc3 = obj.nr*(3/2*c1(end,:) -2*c1(end-1,:)+1/2*c1(end-2,:));
             bc4 = D2*C2*obj.nr*(3/2*c2(end,:) -2*c2(end-1,:)+1/2*c2(end-2,:)) + h0*obj.a.*(1-obj.ep.*obj.a/obj.R)./rhoc.*f;
             bc5 = c2(end,:) -1;
%             
%              F2 = c2rr + k0*(k1*(c1+k2) -c20).*obj.h.^2 - Pe2*(U.*c2r.*obj.h + w.*c2z.*obj.h.^2);
%             
%              bc2 = obj.nr*(-3/2*c2(1,:) +2*c2(2,:)-1/2*c2(3,:));
%     
%              bc4 = obj.nr*(3/2*c2(end,:) -2*c2(end-1,:)+1/2*c2(end-2,:)) + v2*obj.h.*f;
%              bc5 = c2(end,:) -CO2;
             F = [F2;bc2;bc4;bc5];
             
             
         end
         function [c2,f] = get_co2(obj,c20)
             options = optimoptions('fsolve','Display','iter');
             
             [c] = fsolve(@obj.co2fun,c20,options);
             c2 = c(1:obj.nr+1,:);
             
             f = -1e-9/1e-5*2.7e-2./obj.h.*obj.nr.*(3/2*c(end,:) -2*c(end-1,:)+1/2*c(end-2,:));
         end
         function F = c3fun(obj,c)
             h0 = obj.H;
             D1 = 0.9e-9;
             D2 = 1.3e-9;
             g = 9.81;
             nu = 1e-6;
             Pe2 = h0^3*g/(nu*D2);
             Pe1 = 1.3/0.9*Pe2;
             C1 = 5;
             C2 = 2.7e-2;
             km = 2e-4;
             kp = 1e-2;
             rhoc  = 3.69e-5;
            
             
             w = obj.w(2:end-1,:);
             U = obj.U(2:end-1,:);
%                           Fc = fft(c');
%              Fc(33,:) = 0;
%              Fc(log10(abs(Fc))<-4) = 0;
%              c = real(ifft(Fc))';
             c1 = c(1:obj.nr+1,:);

             c2 = c(obj.nr+2:end-1,:);
             %c2 = c(obj.nr+2:end-1,:);
             f = c(end,:);
             c11 = c1(3:end,:);
             c10 = c1(2:end-1,:);
             
             c1n1 = c1(1:end-2,:);
             c1r = obj.nr*(c11 -c1n1)/2 ;
             c1rr = obj.nr^2*(c1n1+c11-2*c10);
             [c1z,~,~,~] =obj.getdiv(c10); 
             
             c21 = c2(3:end,:);
             c20 = c2(2:end-1,:);
             c2n1 = c2(1:end-2,:);
             c2r = obj.nr*(c21 -c2n1)/2 ;
             c2rr = obj.nr^2*(c2n1+c21-2*c20);
             [c2z,~,~,~] =obj.getdiv(c20); 
             
             %F2 = c2rr +obj.ep/obj.R*c2r.*obj.a+ kp*h0^2/D2*(2*C1*km/(kp*C2)*c1 -c20).*obj.a.^2 - obj.ep*Pe2*(U.*c2r.*obj.a + w.*c2z.*obj.a.^2);
             %bc1 = D1*C1*obj.nr*(-3/2*c1(1,:) +2*c1(2,:)-1/2*c1(3,:)) + h0*obj.a.*f;
             %bc1 = D1*C1*obj.nr*(-3/2*c1(1,:) +2*c1(2,:)-1/2*c1(3,:)) - D2*C2*obj.nr*(3/2*c2(end,:) -2*c2(end-1,:)+1/2*c2(end-2,:))*(1 + obj.ep/obj.R);
             %bc2 = obj.nr*(-3/2*c2(1,:) +2*c2(2,:)-1/2*c2(3,:));
             %bc3 = obj.nr*(3/2*c1(end,:) -2*c1(end-1,:)+1/2*c1(end-2,:));
             bc4 = D2*C2*obj.nr*(3/2*c2(end,:) -2*c2(end-1,:)+1/2*c2(end-2,:)) + h0*obj.a.*(1-obj.ep.*obj.a/obj.R)./rhoc.*f;
             %bc5 = c2(end,:) -1;
%             

F1 = c1rr +obj.ep/obj.R*c1r.*obj.a -obj.ep*Pe1*(U.*c1r.*obj.a+w.*c1z.*obj.a.^2);
             F2 = c2rr +obj.ep/obj.R*c1r.*obj.a+ kp*h0^2/D2*(2*C1*km/(kp*C2)*c10 -c20).*obj.a.^2 - obj.ep*Pe2*(U.*c2r.*obj.a + w.*c2z.*obj.a.^2);
             %bc1 = D1*C1*obj.nr*(-3/2*c1(1,:) +2*c1(2,:)-1/2*c1(3,:)) + h0*obj.a.*f;
             bc1 = D1*C1*obj.nr*(-3/2*c1(1,:) +2*c1(2,:)-1/2*c1(3,:)) +  h0*obj.a./rhoc.*(f-mean(f)); %+ h0*(obj.a - 1)/rhoc*mean(f);
             bc2 = obj.nr*(-3/2*c2(1,:) +2*c2(2,:)-1/2*c2(3,:));
             bc3 = obj.nr*(3/2*c1(end,:) -2*c1(end-1,:)+1/2*c1(end-2,:));
             %bc4 = D2*C2*obj.nr*(3/2*c2(end,:) -2*c2(end-1,:)+1/2*c2(end-2,:)) + h0*obj.a.*(1-obj.ep/obj.R).*f;
             bc5 = c2(end,:) -1;
             %F = [F1;F2;bc1;bc2;bc3;bc4;bc5];
             F = [F1;F2;bc1;bc2;bc3;bc4;bc5];
%              F2 = c2rr + k0*(k1*(c1+k2) -c20).*obj.h.^2 - Pe2*(U.*c2r.*obj.h + w.*c2z.*obj.h.^2);
%             
%              bc2 = obj.nr*(-3/2*c2(1,:) +2*c2(2,:)-1/2*c2(3,:));
%     
%              bc4 = obj.nr*(3/2*c2(end,:) -2*c2(end-1,:)+1/2*c2(end-2,:)) + v2*obj.h.*f;
%              bc5 = c2(end,:) -CO2;
             
             
             
         end

         
         function F = co2fun(obj,c)
             D = 1.3e-9;
             h0 = 9.9106e-06;
             C1 = 2*5;
             C2 = 2.7e-2;
             k1 = 2e-4;
             k2 = 1e-2;
             c1 = ones(1,obj.n);
             c2 = c(1:obj.nr+1,:);
             %c2 = c(obj.nr+2:end-1,:);
             %f = c(end,:);
             
             
             c21 = c2(3:end,:);
             c20 = c2(2:end-1,:);
             c2n1 = c2(1:end-2,:);
             c2r = obj.nr*(c21 -c2n1)/2 ;
             c2rr = obj.nr^2*(c2n1+c21-2*c20);
             [c2z,~,~,~] =obj.getdiv(c20); 
             
             F2 = C2*D/h0^2.*c2rr +(C1*k1*c1 -C2*k2*c20).*obj.h.^2 ;
            
             bc2 = obj.nr*(-3/2*c2(1,:) +2*c2(2,:)-1/2*c2(3,:));
    
             %bc4 = D*C2/h0*obj.nr*(3/2*c2(end,:) -2*c2(end-1,:)+1/2*c2(end-2,:)) + v2*obj.h.*f;
             bc5 = c2(end,:) -1;
             F = [F2;bc2;bc5];
         end
         
         function [c] = get_calc(obj,c20)
             options = optimoptions('fsolve','Display','iter');
             
             [c] = fsolve(@obj.calcfun,c20,options);
         end
         
         function F = calcfun(obj,c)
             
             h0 = obj.H;
             D1 = 0.9e-9;
             D2 = 1.3e-9;
             g = 9.81;
             nu = 1e-6;
             Pe2 = h0^3*g/(nu*D2);
             Pe1 = 1.3/0.9*Pe2;

             C1 = 5;
             C2 = 2.7e-2;
             k1 = 2e-4;
             k2 = 1e-2;
             w = obj.w(2:end-1,:);
             U = obj.U(2:end-1,:);
             c1 = c;
             %c2 = c(obj.nr+2:end-1,:);
             %f = c(end,:);
             
             
             c11 = c1(3:end,:);
             c10 = c1(2:end-1,:);
             c1n1 = c1(1:end-2,:);
             c1r = obj.nr*(c11 -c1n1)/2 ;
             c1rr = obj.nr^2*(c1n1+c11-2*c10);
             [c1z,~,~,~] =obj.getdiv(c10); 
          
             F2 = c1rr +obj.ep*(c1r/obj.R.*obj.a-Pe1*(U.*c1r.*obj.a+w.*c1z.*obj.a.^2));
            
             bc2 = D1*C1/h0*obj.nr*(-3/2*c1(1,:) +2*c1(2,:)-1/2*c1(3,:)) + obj.a.*obj.f;
    
             bc4 = obj.nr*(3/2*c1(end,:) -2*c1(end-1,:)+1/2*c1(end-2,:)) ;
             
             F = [F2;bc2;bc4];
         end
         
         function [eta,R] = stalactite_grow(obj,del0,R0,Q,L)
             gamma = 1e12;
             f0 = ones(1,obj.n)*1e-7;
             
             c0 = 1+1/2*f0.*(1-obj.r.^2)*1e-3 +obj.ep*sin(obj.nz*2*pi);
         delt = 0.01;
         T = 1.2;
         obj.del = del0;
         obj = obj.getndparams2(R0,Q,L);
         etao = obj.eta*obj.H;
         
         R = [R0];
         t = 0:delt:T;
         
         
         eta = [etao];
         obj.eta = obj.H*obj.eta;
         for i =1:length(t)
             obj = obj.getndparams2(R0,Q,L);
             obj.eta = etao/obj.H;
             [c,f] = obj.get_conc2(c0,f0,R0,Q,L);
             
             etan = obj.eta+delt*(f-mean(f))*gamma;
             etan = obj.H*etan;
             eta = [eta;etan];
             obj.eta = etan;
             etao = etan;
             Rn = R0+delt*gamma*mean(f)*obj.H;
             
             R = [R;Rn];
             R0 = Rn;
         
             c0 = c;
             f0 = f;
         end
         end
         function [c0,c1,Rt,etat] = stalactite_stab(obj,c0,c1,Rt,etat,R,Q,L)
             
              if nargin<=5
                R = 0.01;
                          Q = 1e-9;
                          L = 1e-2;
              end
                      if nargin ==1
                      c0 = 1+1/2*(obj.r.^2)*1e-3;
                      c1 = 0*c0;
                      Rt = 1;
                      etat = 0;
                      end
                      
                      
                      obj = obj.getndparams2(R,Q,L);
                      
                                   obj = obj.get_h;
             obj.h = obj.h0;
              
             options = optimoptions('fsolve','Display','none');
             cinit = [c0;Rt;c1;etat];
             [c] = fsolve(@obj.stalactite_stab_fun,cinit,options);
             c0 = c(1:end/2-1);
             c1 = c(end/2+1:end-1);
             Rt = c(end/2);
             etat = c(end);
             
             
             
         end
         function [c0,c1,ca,Rt,etat] = stalactite_stab_ca(obj,c0,c1,ca,Rt,etat,R,Q,L)
             
              if nargin<=5
                          R = 0.005;
                          Q = 1e-9;
                          L = 3e-3;
              end
                      if nargin ==1
                      c0 = 1+1/2*(obj.r.^2)*1e-3;
                      c1 = 0*c0;
                      ca = c1;
                      Rt = 1;
                      etat = 0;
                      
                      end
                      
                      
                      obj = obj.getndparams2(R,Q,L);
                      
                                   
              
             options = optimoptions('fsolve','Display','none');
             cinit = [c0;Rt;c1;etat;ca];
             [c] = fsolve(@obj.stalactite_stab_fun_ca,cinit,options);
             
             nr = length(obj.r); 
             c0 =c(1:nr);
             Rt = c(nr+1);
             c1 = c(nr+2:2*nr +1);
             etat = c(2*nr + 2);
             ca = c(2*nr+3:end);
             
             
             
         end
         function F = stalactite_stab_fun(obj,c)
              h0 = obj.H;
             D1 = 0.9e-9;
             D2 = 1.3e-9;
             g = 9.81;
             nu = 1e-6;
             Pe2 = h0^3*g/(nu*D2);
             Pe1 = 1.3/0.9*Pe2;
             C1 = 5;
             C2 = 2.7e-2;
             km = 2e-4;
             kp = 1e-2;
             rhoc  = 3.69e-5;
             
              k = 2*pi/obj.L;
                R = obj.R;
                Re = obj.Re;
                Bo = obj.Bo;
                eps = obj.ep;
                a = (5*eps.*(10*eps.*R.^2.*Bo.^2+30*Bo.^2.*R.^3-eps.*k.^2.*(k.^2.*R.^2-1).*((5*k.^2-2*Bo.*Re).*R.^2-5)))./(225*Bo.^2.*R.^4+150*Bo.^2.*eps.*R.^3+eps.^2.*(25*R.^2.*Bo.^2+k.^2.*(5+(2*Bo.*Re-5*k.^2).*R.^2).^2));
                
                theta = atan(Bo.*k.*R.*(-15*k.^2.*R.^3+15*R+eps.*(15+(-15*k.^2+4*Re.*Bo).*R.^2))./(30*Bo.^2.*R.^3+eps.*(10*Bo.^2.*R.^2-k.^2.*(k.^2.*R.^2-1).*((5*k.^2-2.*Bo.*Re).*R.^2-5))))+(1-sign(a))/2*pi;
                A = 5.*eps.*sqrt((4*R.^2.*Bo.^2+(k-k.^3.*R.^2).^2)./(225*Bo.^2.*R.^4+150*eps.*Bo.^2.*R.^3+eps.^2.*(25*R.^2.*Bo.^2+k.^2.*(5 + (2*Bo.*Re-5*k.^2).*R.^2).^2)));
                 phi = atan(R*A*tan(theta)/(R*A - eps*sec(theta)));
                B = A*R*sin(theta)/((R+eps)*sin(phi));
                h = B*exp(-1i*phi);
                
                

                
             c0 =c(1:end/2-1);
             Rt = c(end/2);
             c1 = c(1+end/2:end-1);
             etat = c(end);
             c0n1 = c0(1:end-2);
             c00  = c0(2:end-1);
             c01 = c0(3:end);
             c0r = obj.nr*(c01 -c0n1)/2 ;
             c0rr = obj.nr^2*(c0n1+c01-2*c00);
             c1n1 = c1(1:end-2);
             c10  = c1(2:end-1);
             c11 = c1(3:end);
             c1r = obj.nr*(c11 -c1n1)/2 ;
             c1rr = obj.nr^2*(c1n1+c11-2*c10);
             r = linspace(0,1,length(c0))';
             w0 = r - r.^2/2 + eps/R*(r/2 - r.^2/2+r.^3/6);
             w0i = r.^2/2- r.^3/6 + eps/R*(r.^2/4 - r.^3/6 + r.^4/24);
             w1i = h*(r.^2-r.^3/3)+eps/24.*(12*h.*(3*r.^2/2-r.^3+3*r.^4/4)/R + 12i*(1+h).*k.*(r.^3/3-r.^2).*(R^2*k^2-1)/(Bo*R^2)+1i*h*k*Re.*(4*r.^2-r.^4+r.^5/5));
             w0 = w0(2:end-1);
             w1i = w1i(2:end-1);
             w0i = w0i(2:end-1);
             F1 = c0rr+eps/R*c0r+kp*h0^2/D2*(2*C1*km/(kp*C2) -c00);
             F2 = c0(end)-1;
             F3 = obj.nr*(-3/2*c0(1,:) +2*c0(2,:)-1/2*c0(3,:));
             F4 = c1rr+eps/R*(c1r+h*c0r)+kp*h0^2/D2*(2*h*(2*C1*km/(kp*C2) -c00)-c10)-1i*k*eps*Pe2*(w0.*c10-(h*w0i+w1i).*c0r);
             F5 = c1(end);
             F6 = obj.nr*(-3/2*c1(1,:) +2*c1(2,:)-1/2*c1(3,:));
             F7 =   D2*C2*obj.nr*(3/2*c0(end,:) -2*c0(end-1,:)+1/2*c0(end-2,:)) + h0*(1-obj.ep/obj.R)./rhoc.*Rt;

             F8 = D2*C2*obj.nr*(3/2*c1(end,:) -2*c1(end-1,:)+1/2*c1(end-2,:)) + h0./rhoc.*((1-obj.ep/obj.R).*(Rt*h+etat)-obj.ep/obj.R*Rt*h);
             
             F = [F1;F2;F3;F4;F5;F6;F7;F8];
         end
         function F = stalactite_stab_fun_ca(obj,c)
              h0 = obj.H;
             D1 = 0.9e-9;
             D2 = 1.3e-9;
             g = 9.81;
             nu = 1e-6;
             Pe2 = h0^3*g/(nu*D2);
             Pe1 = 1.3/0.9*Pe2;
             C1 = 5;
             C2 = 2.7e-2;
             km = 2e-4;
             kp = 1e-2;
             rhoc  = 3.69e-5;
             
              k = 2*pi/obj.L;
                R = obj.R;
                Re = obj.Re;
                Bo = obj.Bo;
                eps = obj.ep;
                a = (5*eps.*(10*eps.*R.^2.*Bo.^2+30*Bo.^2.*R.^3-eps.*k.^2.*(k.^2.*R.^2-1).*((5*k.^2-2*Bo.*Re).*R.^2-5)))./(225*Bo.^2.*R.^4+150*Bo.^2.*eps.*R.^3+eps.^2.*(25*R.^2.*Bo.^2+k.^2.*(5+(2*Bo.*Re-5*k.^2).*R.^2).^2));
                
                theta = atan(Bo.*k.*R.*(-15*k.^2.*R.^3+15*R+eps.*(15+(-15*k.^2+4*Re.*Bo).*R.^2))./(30*Bo.^2.*R.^3+eps.*(10*Bo.^2.*R.^2-k.^2.*(k.^2.*R.^2-1).*((5*k.^2-2.*Bo.*Re).*R.^2-5))))+(1-sign(a))/2*pi;
                A = 5.*eps.*sqrt((4*R.^2.*Bo.^2+(k-k.^3.*R.^2).^2)./(225*Bo.^2.*R.^4+150*eps.*Bo.^2.*R.^3+eps.^2.*(25*R.^2.*Bo.^2+k.^2.*(5 + (2*Bo.*Re-5*k.^2).*R.^2).^2)));
                 phi = atan(R*A*tan(theta)/(R*A - eps*sec(theta)));
                B = A*R*sin(theta)/((R+eps)*sin(phi));
                h = B*exp(-1i*phi);

             nr = length(obj.r); 
             c0 =c(1:nr);
             Rt = c(nr+1);
             c1 = c(nr+2:2*nr +1);
             etat = c(2*nr + 2);
             ca = c(2*nr+3:end);
             
             c0n1 = c0(1:end-2);
             c00  = c0(2:end-1);
             c01 = c0(3:end);
             c0r = obj.nr*(c01 -c0n1)/2 ;
             c0rr = obj.nr^2*(c0n1+c01-2*c00);
             c1n1 = c1(1:end-2);
             c10  = c1(2:end-1);
             c11 = c1(3:end);
             c1r = obj.nr*(c11 -c1n1)/2 ;
             c1rr = obj.nr^2*(c1n1+c11-2*c10);
                          can1 = ca(1:end-2);
             ca0  = ca(2:end-1);
             ca1 = ca(3:end);
             car = obj.nr*(ca1 -can1)/2 ;
             carr = obj.nr^2*(can1+ca1-2*ca0);
             r = linspace(0,1,length(c0))';
             w0 = r - r.^2/2 + eps/R*(r/2 - r.^2/2+r.^3/6);
             w0i = r.^2/2- r.^3/6 + eps/R*(r.^2/4 - r.^3/6 + r.^4/24);
             w1i = h*(r.^2-r.^3/3)+eps/24.*(12*h.*(3*r.^2/2-r.^3+3*r.^4/4)/R + 12i*(1+h).*k.*(r.^3/3-r.^2).*(R^2*k^2-1)/(Bo*R^2)+1i*h*k*Re.*(4*r.^2-r.^4+r.^5/5));
             w0 = w0(2:end-1);
             w1i = w1i(2:end-1);
             w0i = w0i(2:end-1);
             F1 = c0rr+eps/R*c0r+kp*h0^2/D2*(2*C1*km/(kp*C2) -c00);
             F2 = c0(end)-1;
             F3 = obj.nr*(-3/2*c0(1,:) +2*c0(2,:)-1/2*c0(3,:));
             F4 = c1rr+eps/R*(c1r+h*c0r)+kp*h0^2/D2*(2*h*(2*C1*km/(kp*C2) -c00)+2*C1*km/(kp*C2)*ca0-c10)-1i*k*eps*Pe2*(w0.*c10-(h*w0i+w1i).*c0r);
             F5 = c1(end);
             F6 = obj.nr*(-3/2*c1(1,:) +2*c1(2,:)-1/2*c1(3,:));
             F7 =   D2*C2*obj.nr*(3/2*c0(end,:) -2*c0(end-1,:)+1/2*c0(end-2,:)) + h0*(1-obj.ep/obj.R)./rhoc.*Rt;

             F8 = D2*C2*obj.nr*(3/2*c1(end,:) -2*c1(end-1,:)+1/2*c1(end-2,:)) + h0/rhoc.*((1-obj.ep/obj.R).*(Rt*h+etat)-obj.ep/obj.R*Rt*h);
             F9 = carr+eps/R*(car)-1i*k*eps*Pe1*(w0.*ca0);
             F10 =  obj.nr*(3/2*ca(end,:) -2*ca(end-1,:)+1/2*ca(end-2,:));
             F11 = D1*C1*obj.nr*(-3/2*ca(1,:) +2*ca(2,:)-1/2*ca(3,:)) + h0/rhoc*(Rt*h+etat);
             F = [F1;F2;F3;F4;F5;F6;F7;F8;F9;F10;F11];
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


