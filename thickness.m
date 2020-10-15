classdef thickness
    properties
        L = 2*pi %wavelength
        A = 1 %amplitude
        B = 1 %Bond number like non dimensional value
        n = 100 %number of points to plot
        z %domain
        eta %shape of the wall
        etaz %first derivative
        etazzz %third derivative
        a2 %fluid thickness at 2nd order (where it first varies)
        delta = 0.1 % small parameter
        a0 = 1 %first thickness
        a % version of a2 where we want the wavelength to be normalised
        k = 1 % wavenumber
        nz % normalised z
    end
    methods
        function obj = thickness(L,A,B,n)
            
            if nargin>0
                obj.L = L;
                obj.A = A;
                obj.B = B;
                obj.n = n;
            end
            obj = obj.get_z;
            obj = obj.get_thickness;
            obj.nz = -1/2:1/obj.n:1/2;
            
            
        end
        function obj = get_z(obj)
            obj.z =-obj.L/2:obj.L/obj.n:obj.L/2;
        end
        function obj = get_thickness(obj)
            obj.eta = obj.A.*sin(2*pi./obj.L.*obj.z);
            obj.etaz = 2*pi./obj.L.*obj.A.*cos(2*pi./obj.L.*obj.z);
            obj.etazzz = -(2*pi./obj.L).^3.*obj.A.*cos(2*pi./obj.L.*obj.z);
            obj.a2 = obj.a0/3*(2*obj.etaz.^2- obj.eta - obj.B.*obj.etazzz);
            
            obj.a = obj.A/3*obj.a0*(2.*obj.k.^2.*obj.A.*cos(obj.nz*2*pi).^2-sin(2*pi*obj.nz) +obj.B.*obj.k.^3.*cos(2*pi*obj.nz));
        end
        function plot_B(obj,B)
            clf, hold on
            
            for b = B
                obj.B = b;
                obj = obj.get_thickness;
                
                plot(obj.a2,obj.z,'DisplayName',sprintf('$B = %g$',b))
            end
            plot(obj.eta,obj.z,'--','DisplayName', 'Wall Shape')
            flip_y
            legend
            xlabel('$a_2$')
            ylabel('$z$')
            
            
        end
        function plot_L(obj,L)
            
            clf, hold on
            
            for l = L
                obj.L = l;
                obj = obj.get_z;
                obj = obj.get_thickness;
                
                plot(obj.a2,obj.z./obj.L,'DisplayName',sprintf('$\\lambda = %g$',l))
            end
            plot(obj.eta,obj.z/obj.L,'--','DisplayName', 'Wall Shape')
            flip_y
            legend
            xlabel('$a_2$')
            ylabel('$z$')
            
            
        end
        function plot_A(obj,A)
            clf, hold on
            
            for a = A
                obj.A = a;
                obj = obj.get_thickness;
                
                plot(obj.a2,obj.z,'DisplayName',sprintf('$A = %g$',a))
            end
            plot(obj.eta,obj.z,'--','DisplayName', 'Wall Shape')
            flip_y
            xlabel('$a_2$')
            ylabel('$z$')
            legend
        end
        
        function surf_B(obj)
            [obj.B,obj.z] = meshgrid(0.1:0.01:5,obj.z);
            obj = obj.get_thickness;
            clf
            surf(obj.B,obj.z,obj.a2);
            xlabel('Bond number');
            ylabel('$z$');
            zlabel('Thickness');
            shading interp;
            
        end
        function surf_A(obj)
            [obj.A,obj.z] = meshgrid(0.5:0.01:2,obj.z);
            obj = obj.get_thickness;
            clf
            surf(obj.A,obj.z,obj.a2);
            xlabel('Amplitude');
            ylabel('$z$');
            zlabel('Thickness');
            shading interp;
        end
        function surf_L(obj)
            [k,z] = meshgrid(0.1:0.01:3,-0.5:0.01:0.5);
            a = obj.A/3*obj.a0*(2.*k.^2.*obj.A.*cos(z*2*pi).^2-sin(2*pi*z) +obj.B.*k.^3.*cos(2*pi*z));
            surf(k,z,a)
            xlabel('$k$');
            ylabel('$\frac{kz}{2\pi}$');
            zlabel('Thickness');
            shading interp;
        end
        function plot_end_root(obj,B)
            clf, hold on 
            for b=B
                obj.B = b;
                obj.k = 2*obj.A/b;
                obj = obj.get_thickness;
                plot(obj.a,obj.z,'DisplayName',sprintf('$B = %g$',b))
            end
            plot(obj.eta,obj.z,'--','DisplayName', 'Wall Shape')
            flip_y
            legend
            xlabel('$a_2$')
            ylabel('$z$')
            
        end
        function plot_end_rootk(obj,k)
            clf, hold on 
            for K=k
                obj.k = K;
                obj.A = 1/2*obj.B*K;
                obj = obj.get_thickness;
                plot(obj.a,obj.z,'DisplayName',sprintf('$k = %g$',K))
            end
            plot(obj.eta,obj.z,'--','DisplayName', 'Wall Shape')
            flip_y
            legend
            xlabel('$a_2$')
            ylabel('$z$')
            
        end
        function plot_end_rootA(obj,A)
            clf, hold on 
            for a=A
                obj.A = a;
                obj.B = 2*a/obj.k;
                obj = obj.get_thickness;
                plot(obj.a,obj.z,'DisplayName',sprintf('$A = %g$',a))
            end
            plot(obj.eta,obj.z,'--','DisplayName', 'Wall Shape')
            flip_y
            legend
            xlabel('$a_2$')
            ylabel('$z$')
            
        end
        function turning_root(obj,k)
            clf, hold on 
            for K= k
                obj.k = K;
                obj.A = 1/(4*K^2);
                obj.B = 1/(2*K^3);
                obj = obj.get_thickness;
                plot(obj.a,obj.nz,'DisplayName',sprintf('$k = %g$',K))
            end
            plot(sin(2*pi*obj.nz),obj.nz,'--','DisplayName', 'Wall Shape')
            flip_y
            legend
            xlabel('$a_2$')
            ylabel('$z$')
        end
        
        function early_root(obj,k)
            clf, hold on
            for K=k
                obj.k = K;
                obj.B = (obj.A*2*sqrt(2)*K^2 -1)/K^3;
                if obj.B<0
                    continue
                end
                obj = obj.get_thickness;
                plot(obj.a,obj.nz,'DisplayName',sprintf('$k = %g$, $B =%.3g$',K,obj.B))
            end
            plot(sin(2*pi*obj.nz),obj.nz,'--','DisplayName', 'Wall Shape')
            flip_y
            legend
            xlabel('$a_2$')
            ylabel('$z$')
        end
                
            
    end
end