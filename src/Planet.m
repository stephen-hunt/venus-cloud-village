classdef Planet
    %PLANET Definition of a planet and its parameters
    
    properties
        R;     %m
        g_o;   %m/s^2
        H;     %m
        rho_o; %kg/m^3
        c;
    end
    
    methods
        function obj = Planet(R, g_o, H, rho_o, c)
            %PLANET Construct an instance of this class
            obj.R     = R;
            obj.g_o   = g_o;
            obj.H     = H;
            obj.rho_o = rho_o;
            obj.c     = c;
        end
        
        function rho = rho(self, h)
            %RHO Density at a given altitude
            rho = self.rho_o * exp(-h/self.H);
        end
        
        function r = r(self, h)
            %r Distance from center of planet
            r = self.R + h;
        end
    end
end

