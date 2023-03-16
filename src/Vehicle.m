classdef Vehicle
    %VEHICLE Definition of a entry vehicle and its parameters
    
    properties
        m;                   %kg
        Cd;
        D_Vehicle;           %m
        S_Vehicle;           %m^2
        BC;                  %kg/m^2
        epsilon;
        Cx;
        D_flotation;         %m
        S_flotation;         %m^2
        Cd_flotation;
        D_flotation_initial; %m
    end
    
    methods
        function obj = Vehicle(m, D_Vehicle, C_D, epsilon, Cx, D_flotation, D_flotation_initial, Cd_flotation)
            %VEHICLE Construct an instance of this class
            obj.m = m;
            obj.Cd = C_D;
            obj.D_Vehicle = D_Vehicle;
            obj.S_Vehicle = pi * (D_Vehicle/2)^2;
            obj.BC = m/(C_D*obj.S_Vehicle);
            obj.epsilon = epsilon;
            obj.Cx = Cx;
            obj.D_flotation = D_flotation;
            obj.S_flotation = pi * (D_flotation/2)^2;
            obj.Cd_flotation = Cd_flotation;
            obj.D_flotation_initial = D_flotation_initial;
        end
    end
end

