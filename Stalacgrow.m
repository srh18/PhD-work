classdef Stalacgrow
    properties
        z = 0:pi/100:2*pi
        etai
        detai
        ddetai
        dddetai 
        B = 1
        pco2 = 1
        D = 1
        D1 = 1
        Pe = 1
        Pe1 = 1
        k1 = 1
        k2 = 1
        g0
        c0
        c2
        g2
        cinit = 1
        eta0
        t = 1
        eta2
        vm = 1
        a2
        
        
        
        
    end
    
    methods
        function obj = Stalacgrow
        obj.etai = sin(obj.z);
        obj.detai = cos(obj.z);
        obj.ddetai = - sin(obj.z);
        obj.dddetai = - cos(obj.z);
        obj.eta0 = obj.etai;
        obj = obj.get_thickness;
        obj = obj.get_chem;
        
        end
        function obj = get_thickness(obj)
        obj.a2 = 1/3*(2*obj.detai.^2 - obj.eta0 - obj.B*obj.dddetai);
        end
        function obj = get_chem(obj)
            obj.g0 = obj.pco2;
            obj.c0 = obj.cinit - 3*obj.D1*obj.g0*obj.k1*obj.z/(obj.Pe*obj.D);
        end
        function obj = get_wall(obj)
            obj = obj.get_thickness;
            obj.eta0 = obj.D1*obj.pco2*obj.k1*obj.vm*obj.t+obj.eta0;
            obj.eta2 = -obj.D1*(obj.k1*obj.pco2*(2*obj.a2-1/3*obj.k1-obj.detai.^2)-2*obj.k2*obj.c0)*obj.t- obj.D1*obj.pco2*obj.k1*obj.vm*obj.t^2;
        end
        function plot_wall(obj)
            for t = 0:0.1:10
                obj.t = t;
                obj = obj.get_wall;
                clf, hold on 
                plot(obj.eta0,obj.z)
                
                plot(obj.eta0 +0.01*obj.eta2,obj.z)
                ax = gca;
                ax.YDir ='reverse';
                pause(0.1)

            end
        end
    end
        
            
        
        
        
end
