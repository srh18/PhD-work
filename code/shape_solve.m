classdef shape_solve
    
    %Class to help solve Thin film down a cylinder of similar size to the
    %wavelength
    
    properties
        n {mustBeInteger,mustBePositive}= 400  %number of grid points INT
        
        % Fluid properties Parameters
        
        
        Bo {mustBeNonnegative} = 1 % Bond number
        Re {mustBeNonnegative}= 1 % Reynolds number
        ep {mustBeNonnegative}= 1e-3 % Ratio between the fluid thickness and the radius of the cylinder
        R {mustBeNonnegative} = 1 % Radius of cylinder
        
        q {mustBeNonnegative}= 1/3 % flow rate
        
        %Wall shape parameters
        
        del {mustBeNonnegative}= 1 % Disturbance amplitude
        L {mustBeNonnegative}= 2*pi % wavelength of disturbance
        l {mustBeNonnegative,mustBeLessThanOrEqual(l,1)} = 0.5 %width of the step
        del2 {mustBeNonnegative} = 0.1 % steepness of the step
        
        
        % Asthetic features
        
        periods {mustBePositive,mustBeInteger} = 1 % number of periods to plot INT - ignored for now
        
        %Derived from parameters
        
        z % domain
        a % fluid thickness
        eta % wall shape
        
        %Toggles
        paramtoggle  {mustBeMember(paramtoggle,[1,2,3,4,5,6,7])} = 1 % 1 for Bo, 2 for R, 3 for L, 4 for Re, 5 for ep, 6 for del, 7 for l
        flow {mustBeMember(flow,[0,1])} = 1 % 0 if keeping volume constant - 1 if flow rate is constant
        fdo {mustBeMember(fdo,[2,4])}= 2 % Order of the finite difference scheme used (may add later)
        wall_shape {mustBeMember(wall_shape,[0,1,2,3,4])} = 0 % 0 for cosine, 1 step, 2 sawtooth, 3 for flat, 4 user input
        boundary {mustBeMember(boundary,[0,1])} = 0 % 0 for periodic boundary conditions, 1 for not
        
        %dynamic parameters
        delt {mustBePositive} = 0.001 % Time step
        
    end
    
    properties(Hidden)
        %additional outputs from fsolve
        Fval
        eflag
        out
        jac
        %string of parameters
        params = {'Bo','R','L','Re','\\epsilon','\\delta'}
        % Properties of the fluid
        
        amax % maximum fluid thickness
        amin % minimum fluid thickness
        diff % range of fluid thickness
        amaxloc % location of the maximum
        aminloc % location of the minimum
        maxmindiff % distance between the maximum and the minimum
        minmaxdiff % distance between the minimum and the maximum
        
        %derivatives of eta to speed things up
        etaz
        etazz
        etazzz
        etazzzz
        integration {mustBeMember(integration,[0,1,2])} = 0 % 0 is crank Nicholson 1 implicit 2 explicit
    end
    properties(Dependent = true,Hidden)
        S %position from origin
        nz %normalised z
        h % effective thickness
        mass %integral of h
        h2norm % integral of h^2
        amass % integral of a
        a2norm %integral of a^2
        Q % flow rate
        Qint % integral of flow rate
    end
    
    
    methods
        function obj = shape_solve(n,L,Bo,Re, ep, del, periods)
            
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
        function value = get.S(obj)
            value = obj.a+ obj.eta;
        end
        function value = get.nz(obj)
            value = obj.z/obj.L;
        end
        function value = get.h(obj)
            value = obj.a + obj.ep/obj.R*(obj.a.*obj.eta+ obj.a.^2/2);
        end
        function value = get.mass(obj)
            value = sum(obj.h,2)/obj.n*obj.L;
        end
        function value = get.h2norm(obj)
            value = sum(obj.h.^2,2)/obj.n*obj.L;
        end
        function value = get.a2norm(obj)
            value = sum(obj.a.^2,2)/obj.n*obj.L;
        end
        function value = get.amass(obj)
            value = sum(obj.a,2)/obj.n*obj.L;
        end
        %dont work
%         function value = get.Q(obj) 
%             hz= 0*obj.h;
%             hzzz = 0 *obj.h;
%             for i = 1:length(obj.h)
%                 
%                 [hz(i,:),~,hzzz(i,:),~] = obj.getdiv(obj.h(i,:));
%             end
%             [etaz,~,etazzz,~] = obj.getdiv(obj.eta);
%             
%             value = obj.h.^3/3+obj.ep*(2/15*obj.Re*obj.h.^6.*hz+obj.h.^3/(3*obj.Bo)*((hz+etaz)/obj.R^2+hzzz+etazzz)-s/3*obj.h.^3.*obj.eta-1/6*obj.h^4);
%         end
%         function value = get.Qint(obj)
%             value = trapz(obj.Q,2);
%         end
            
            
        
        
        %Functions which also change other parameters when changed
        function obj = set.L(obj,value)
            obj.L = value;
            obj = obj.reset;
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
        
        function obj = set.wall_shape(obj,value)
            obj.wall_shape = value;
            if obj.wall_shape <4
                obj = obj.get_wall;
            end
        end
        
        function obj = reset(obj)
            %realign z, wall if n or L are changed
            obj.z = 0: obj.L/obj.n: obj.periods*obj.L- obj.L/obj.n;
            obj = obj.get_wall;
            
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
                obj.eta = (obj.z/obj.L-(1-obj.l)/2)/(obj.l*obj.L).*(obj.del/2*(tanh((obj.z/obj.L-(1-obj.l)/2)/obj.del2) - tanh((obj.z-obj.L*(1+obj.l)/2)/obj.del2)));
            elseif obj.wall_shape == 0
                
                obj.eta = obj.del*cos(2*pi/obj.L*obj.z);
            elseif obj.wall_shape == 3
                obj.eta = 0*obj.z;
            else
                warning('Wall Shape may no longer match coordinates')
            end
            
        end
        
        
        %MAIN FUNCTIONS
        
        function ainit = linear_shape(obj)
            % warning uses old equation may be incorrect.
            %analytic solution to first order of a disturbance
            %del*cos(2*pi/L) (where del<<1) used as a initial guess for
            %fsolve
            A = -(2*pi/obj.L)^3+(1-9*obj.Bo*obj.Re/40)*2*pi/obj.L;
            B = 4*obj.Bo + 3*obj.Bo/obj.ep;
            g1 = - 3*obj.Bo;
            g2 = 2*pi/obj.L*(1 - (2*pi/obj.L)^2) ;
            
            a = 1/(A^2+B^2)*(A*g1+B*g2);
            b = 1/(A^2+B^2)*(B*g1 - A* g2);
            
            
            ainit = 1+ a*obj.del*sin(2*pi/obj.L*obj.z)+b*obj.del*cos(2*pi/obj.L*obj.z);
            
        end
        
        function a = small_ep_shape(obj)
            %warning uses old equation may be incorrect
            a = 1 - obj.ep/3*obj.del*(cos(2*pi/obj.L*obj.z)+1/obj.Bo*(-2*pi/obj.L+8*pi^3/obj.L^3)*sin(2*pi/obj.L*obj.z));
        end
        
                function F = afun(obj,a)
            [az,azz,azzz,azzzz] = obj.getdiv(a);
            [etaz,etazz,etazzz,etazzzz] = obj.getdiv(obj.eta);
            F = az.*a.^2 + obj.ep.*(4*az.*a.^3/(6.*obj.R)+2/15.*obj.Re.*(a.^6.*azz+6.*a.^5.*az.^2)+1/(3.*obj.Bo).*(a.^3.*((azz+etazz)/obj.R^2+azzzz+etazzzz)+3*az.*a.^2.*((az+etaz)/obj.R^2+azzz+etazzz))+etaz.*a.^3/(3*obj.R));
            
        end
            function [az,azz,azzz,azzzz] = getdiv(obj,a)
                
             an2 = a([end-1 end 1:end-2]);
                an1 = a([end 1:end-1]);
                a1 = a([2:end 1]);
                a2 = a([3:end 1 2]);
                az = obj.n/obj.L*(-1/2*an1+1/2*a1);
                azz =(obj.n/obj.L)^2*(an1-2*a+ a1);
                azzz= (obj.n/obj.L)^3*(-1/2*an2+an1-a1+1/2*a2);
                azzzz =(obj.n/obj.L)^4 *(an2-4*an1+6*a-4*a1+a2);
            end
        

        
        function obj = get_a(obj,optionon,ainit)
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
            
            [obj.a,obj.Fval,obj.eflag,obj.out,obj.jac] = fsolve(@obj.afun,ainit,options);
            
        end
        function plot_param(obj, values, param,etatoggle)
            %Vary the amplitude of the disturbance
            obj.paramtoggle = param;
            clf, hold on
            if nargin <4
                etatoggle = 0;
                if nargin ==1
                    
                    values = obj.getp;
                end
            end
            for val = values
                obj = obj.setp(val);
                
                obj = obj.get_a;
                if obj.eflag < 1
                    warning('No solution found for %s = %g\n',obj.params{obj.paramtoggle},val)
                else
                    plot(obj.a+etatoggle*obj.eta,obj.nz,'DisplayName',sprintf('$%s = %g$',obj.params{obj.paramtoggle},val))
                end
            end
            if etatoggle ==1
                plot(obj.eta,obj.z/obj.L,'DisplayName','Wall')
            end
            legend
            flip_y
        end
        

        
        function plot_del(obj,dels,etatoggle)
            %Vary the amplitude of the disturbance
            clf, hold on
            if nargin <3
                etatoggle = 0;
                if nargin ==1
                    dels = obj.del;
                end
            end
            for del = dels
                obj.del = del;
                
                obj = obj.get_a;
                if obj.eflag < 1
                    warning('No solution found for del = %g\n',del)
                else
                    plot(obj.a+etatoggle*obj.eta,obj.z,'DisplayName',sprintf('$\\delta = %g$',del))
                end
            end
            if etatoggle ==1
                plot(obj.eta,obj.z/obj.L,'DisplayName','Wall')
            end
            legend
            flip_y
        end
        
        function plot_L(obj,Ls,etatoggle)
            %Vary the wavelength of the disturbance
            clf, hold on
            if nargin <3
                etatoggle = 0;
                if nargin ==1
                    Ls = obj.L;
                end
            end
            
            
            for L = Ls
                obj.L = L;
                
                obj = obj.get_a;
                if obj.eflag < 1
                    warning('No solution found for L = %g\n',L)
                else
                    plot(obj.a+etatoggle*obj.eta,obj.z/obj.L,'DisplayName',sprintf('$L = %g$',L))
                    
                end
            end
            if etatoggle ==1
                plot(obj.eta,obj.z/obj.L,'DisplayName','Wall')
            end
            legend
            flip_y
        end
        
        function plot_Re(obj,Res,etatoggle)
            %vary the Reynold's number
            clf, hold on
            if nargin <3
                etatoggle = 0;
                if nargin ==1
                    Res = obj.Re;
                end
            end
            for Re = Res
                obj.Re = Re;
                obj = obj.get_a;
                if obj.eflag < 1
                    warning('No solution found for Re = %g\n',Re)
                else
                    plot(obj.a+etatoggle*obj.eta,obj.z/obj.L,'DisplayName',sprintf('$Re = %g$',Re))
                end
            end
            if etatoggle ==1
                plot(obj.eta,obj.z/obj.L,'DisplayName','Wall')
            end
            legend
            flip_y
        end
        
        function plot_Bo(obj,Bos,etatoggle)
            %Vary the Bond number
            clf, hold on
            if nargin <3
                etatoggle = 0;
                if nargin ==1
                    Bos = obj.Bo;
                end
            end
            for Bo = Bos
                obj.Bo = Bo;
                obj = obj.get_a;
                if obj.eflag < 1
                    warning('No solution found for Bo = %g\n',Bo)
                else
                    plot(obj.a+etatoggle*obj.eta,obj.z/obj.L,'DisplayName',sprintf('$Bo = %g$',Bo))
                end
            end
            if etatoggle ==1
                plot(obj.eta,obj.z/obj.L,'DisplayName','Wall')
            end
            legend
            flip_y
        end
        
        function plot_ep(obj,eps,etatoggle)
            % vary the ratio beween the fluid thickness and the
            % radius
            clf, hold on
            if nargin <3
                etatoggle = 0;
                if nargin ==1
                    eps = obj.ep;
                end
            end
            
            for ep = eps
                obj.ep = ep;
                obj = obj.get_a;
                if obj.eflag < 1
                    warning('No solution found for ep = %g\n',ep)
                else
                    plot(obj.a+etatoggle*obj.eta,obj.z/obj.L,'DisplayName',sprintf('$\\epsilon = %g$',ep))
                end
            end
            if etatoggle ==1
                plot(obj.eta,obj.z/obj.L,'DisplayName','Wall')
            end
            legend
            flip_y
        end
        
        function plot_q(obj,qs,etatoggle)
            % vary the ratio beween the fluid thickness and the
            % radius
            hold on
            if nargin <3
                etatoggle = 0;
                if nargin ==1
                    qs = obj.q;
                end
            end
            for q = qs
                obj.q = q;
                
                obj = obj.get_a;
                if obj.eflag < 1
                    warning('No solution found for q = %g\n',q)
                else
                    plot(obj.a+etatoggle*obj.eta,obj.z/obj.L,'DisplayName',sprintf('$q= %g$',q))
                end
            end
            if etatoggle ==1
                plot(obj.eta,obj.z/obj.L,'DisplayName','Wall')
            end
            legend
            flip_y
        end
        
        function plot_l(obj,ls)
            clf, hold on
            fig = gca;
            if nargin ==1
                ls = obj.l;
            end
            for l = ls
                obj.l = l;
                obj = obj.get_a;
                if obj.eflag < 1
                    warning('No solution found for l = %g\n',l)
                else
                    plot(obj.a+obj.eta,obj.z/obj.L,'DisplayName',sprintf('$l= %g$',l))
                    colour = fig.Children(1).Color;
                    plot(obj.eta,obj.z/obj.L,'--','Color', colour,'DisplayName',sprintf('Wall with $l= %g$',l))
                    
                end
            end
            legend
            flip_y
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
            [obj.amax,i] = max(obj.a,[],2);
            [obj.amin,j] = min(obj.a,[],2);
            obj.diff = obj.amax-obj.amin;
            obj.amaxloc = obj.z(i);
            obj.aminloc = obj.z(j);
            if obj.amaxloc<obj.aminloc
                obj.maxmindiff = obj.aminloc-obj.amaxloc;
                obj.minmaxdiff = obj.L-obj.aminloc+obj.amaxloc;
            else
                obj.maxmindiff = obj.L-obj.amaxloc +obj.aminloc;
                obj.minmaxdiff = obj.amaxloc - obj.aminloc;
                %                 end
            end
        end
        
        function plot_minmax(obj,x,param)
            obj.paramtoggle = param
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
                obj = obj.get_a;
                obj = obj.get_fluid_features;
                max = [max obj.amax];
                min = [min obj.amin];
                diff = [diff obj.diff];
                maxloc = [maxloc obj.amaxloc/obj.L];
                minloc = [minloc obj.aminloc/obj.L];
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
        function [F, J] = tfun(obj,a,aold)
            F = zeros(1,obj.n);
            
                [~,~,~,azzzz] = obj.getdiv(a);
                [az,azz,azzz,~] = obj.getdiv(aold);
                %[etaz,etazz,etazzz,etazzzz] = obj.getdiv(obj.eta);
                F = (a-aold)./obj.delt+aold.^2.*az+obj.ep.*(aold.^3./(3*obj.Bo).*((obj.etazz+azz)./obj.R^2+obj.etazzzz)+aold.^2.*az./obj.Bo.*((obj.etaz+az)/obj.R^2+obj.etazzz+azzz)+2*aold.^3.*az/(3*obj.R)+2/15*obj.Re*(aold.^6.*azz+aold.^5.*az.^2)+a.^3.*azzzz/(3*obj.Bo))+obj.etaz.*obj.ep.*aold.^3/obj.R;%maybe turn a to aold
                %F = (a-aold)./obj.delt+aold.^2.*az+obj.ep.*(aold.^3./(3*obj.Bo).*((obj.etazz+azz)./obj.R^2+obj.etazzzz)+aold.^2.*az./obj.Bo.*((obj.etaz+az)/obj.R^2+obj.etazzz+azzz)+2*aold.^3.*az/(3*obj.R)+2/15*obj.Re*(aold.^6.*azz+aold.^5.*az.^2)+aold.^3.*azzzz/(3*obj.Bo))+obj.etaz.*obj.ep.*aold.^3/obj.R;
            if nargout>1
                an2 = a([end-1 end 1:end-2]);
                an1 = a([end 1:end-1]);
                a1 = a([2:end 1]);
                a2 = a([3:end 1 2]);
                
                A = diag(1/obj.delt + obj.ep/obj.Bo*((obj.n)/obj.L)^4.*((an2-4*an1+8*a-4*a1+a2).*a.^2));
                B = obj.ep/(3*obj.Bo)*((obj.n)/obj.L)^4*((diag(a(1:end-2).^3,-2)+diag(a(end-1:end).^3,obj.n-2))-4*(diag(a(1:end-1).^3,-1)+diag(a(end).^3,obj.n-1))-4*(diag(a(2:end).^3,1)+diag(a(1).^3,1-obj.n))+(diag(a(3:end).^3,2)+diag(a(1:2).^3,2-obj.n)));
                J = A + B;
            end
        end

        function obj = dynamics(obj,T,init)
            t = 0: obj.delt:T;
            
            options = optimoptions('fsolve','Display', 'none','SpecifyObjectiveGradient',true);
            obj.a = zeros(length(t),obj.n);
            
            [obj.a(1,:),obj.Fval,obj.eflag,obj.out,obj.jac] = fsolve(@(a)obj.tfun(a,init),init,options);
            if obj.eflag < 1
                error('No solution found for t = %g\n',0)
            end
            
            for i = 2:length(t)
                
                [obj.a(i,:),obj.Fval,obj.eflag,obj.out,obj.jac] = fsolve(@(a)obj.tfun(a,obj.a(i-1,:)),obj.a(i-1,:),options);
                if obj.eflag < 1
                    warning('No solution found for t = %g\n',t(i))
                    break
                    
                end
                
            end
        end
        function animate(obj,speed,wall)

            if nargin <= 2
                wall = 1;
            end
            data = obj.a;
            if iscell(data)
                l = length(data{1});
                m = length(data);
            else
                l = length(data);
                m= 0;
            end
            %
            % xmax = max(max(data));
            % xmin = min(min(data));
            for i = 1:speed:l
                clf, hold on
                flip_y
                if (obj.wall_shape ~= 3 && wall ~=0)
                    plot(obj.eta,obj.z)
                end
                if m==0
                    
                    plot(data(i,:)+wall*obj.eta,obj.z)
                else
                    for j = 1:m
                        plot(data{j}(i,:),obj.z)
                    end
                end
                plot_n_periods
                plot_minx,plot_maxx
                %xlim([xmin, xmax])
                
                title(sprintf('$t =%g$',obj.delt*i))
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
        function plot_features(obj,T)
            obj = obj.get_fluid_features;
            figure(10),  hold on, title("Amplitude Changes with time")
            figure(11),  hold on, title('Position of maxima/minima')
            figure(12),  hold on , title('Drift of results')
            d = obj.del;
            t = 0:obj.delt:T;
            figure(10)
            
                plot(t,obj.diff, 'DisplayName',sprintf('$\\delta = %g$',d))
                figure(11)
                plot(t,obj.amaxloc, 'DisplayName',sprintf('maxima $ \\delta = %g$',d))
                fig = gca;
                plot(t,obj.aminloc,'--','Color', fig.Children(1).Color, 'DisplayName',sprintf('minima$ \\delta = %g$',d))
                figure(12)
                plot(t,obj.amax, 'DisplayName',sprintf('$\\delta = %g$',d))
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