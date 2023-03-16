classdef EDLSim
    %EDLSIM Simulation of an EDL Sequence
    
    properties
        % Constants
        G_E = 9.8 % m/s^2
        KtoC = -273.15;
        TIME_STEP = .5; %s
        
        % Sim variables
        planet;
        vehicle;
        h_e;
        LtD;
        FPA;
        trim;
        V_e;
        termDescentAlt;
        floatDeployAlt;
        landAlt;
    end
    
    methods
        function obj = EDLSim(h_e, V_e, LtD, FPA, trim, planet, vehicle, termDescentAlt, floatDeployAlt, landAlt)
            %EDLSim Construct an instance of this class
            obj.planet         = planet;
            obj.vehicle        = vehicle;
            obj.h_e            = h_e;
            obj.LtD            = LtD;
            obj.FPA            = deg2rad(FPA);
            obj.trim           = deg2rad(trim);
            obj.V_e            = V_e;
            obj.termDescentAlt = termDescentAlt;
            obj.landAlt        = landAlt;
            obj.floatDeployAlt = floatDeployAlt;
        end
        
        function plotVAvsH(self)
            %plotVAvsH Plot the velocity and acceleratoin as a function of
            %altitude
            
            figure(1);
            
            % Calculate Velocity and Acceleration during hypersonic
            % entry
            [hv, V] = calculateVelocity(self);
            [ha, a] = calculateAcceleration(self);
            
            % Reverse order of results
            V = fliplr(V);
            a = fliplr(a);
            hv = fliplr(hv);
            ha = fliplr(ha);
            
            % Calculate Velocity and Acceleration during terminal
            % descent
            [ht, Vt, at] = calculateTerminalDescent(self);
            
            % Get acceleration in g's
            at = (at - self.planet.g_o) ./ self.G_E;
            
            % Append both hypersonic and terminal descent
            hv = [hv, ht.'];
            ha = [ha, ht.'];
            V = [V, Vt.'];
            a = [a, at.'];
            
            % Print final velocity
            fprintf("Final Velocity: %d\n", V(end));
            fprintf("Maximum Acceleration: %d\n", min(a));
            
            % Plot the results
            title('Atmospheric Entry and Descent');
            
            yyaxis left;
            plot(hv./1000, V);
            xlabel('Altitude (km)');
            ylabel('Velocity (m/s)');
            axis([self.landAlt/1000 self.h_e/1000 0 8000])
            
            yyaxis right;
            plot(ha./1000, a);
            ylabel('Acceleration (g''s)');
            axis([self.landAlt/1000 self.h_e/1000 -3 .5])
            set(gca, 'XDir','reverse');
        end
        
        function plotDRDvsH(self)
            %plotDRDvsH plot the down range distance vs altitude
            
            figure(2);
            
            % Calculate the DRD
            [h, DRD] = calculateDRD(self);
            
            % Plot the results
            ln = plot(h, DRD./1000);
            
            % Set plot values
            ln.Color = [0.85 0.33 0.10];
            title('Down Range Distance vs Altitude');
            xlabel('Altitude (km)');
            ylabel('Down Range Distance (km)');
            ax = gca;
            ax.YRuler.Exponent = 0;
            ax.YColor = [0.85 0.33 0.10];
            set(gca, 'XDir','reverse');
        end
        
        function plotqvsH(self)
            %plotqvsH plot heat flux vs altitude
            
            figure(3);
            
            % Calculate the heat flux
            [h, q] = calculateq(self);
            
            % Plot the results
            ln = plot(h, q);
            
            % Set plot values
            ln.Color = [0.85 0.33 0.10];
            title('Heat Flux vs Altitude');
            xlabel('Altitude (km)');
            ylabel('Convective Heat Flux (W/cm^2)');
            ax = gca;
            ax.YColor = [0.85 0.33 0.10];
            set(gca, 'XDir','reverse');
        end
        
        function printPeakWallTemp(self)
            %printPeakWallTemp print the peak wall temperature of the
            %vehicle
            
            % Don't display this figure
            f = figure('visible', 'off');
            
            % Calculate the max heat flux
            [h, q] = calculateq(self);
            q_max = max(q);
            
            % Calculate and print the resulting wall temperature
            Tw = 1000 * (q_max / (5.67 * self.vehicle.epsilon))^.25;
            fprintf("q_max = %d\n", q_max);
            fprintf("Tw = %d degrees C\n", Tw + self.KtoC);
            
            close(f);
        end
        
        function [h, V] = calculateVelocity(self)
            %calculateVelocity get velocity values vs altitude
            fp = fplot(@(x) self.v(x), [self.termDescentAlt self.h_e]);
            V = fp.YData;
            h = fp.XData;
        end
        
        function [h, a] = calculateAcceleration(self)
            %calculateAcceleration get acceleration values vs altitude
            fp = fplot(@(x) self.a(x), [self.termDescentAlt self.h_e]);
            a = fp.YData;
            h = fp.XData;
        end
        
        function [h, V, a] = calculateTerminalDescent(self)
            %calculateTerminalDescent get acceleration and velocity
            %values vs altitude during terminal descent
            
            %Initialize variables
            h = zeros(1000, 1);
            V = zeros(1000, 1);
            a = zeros(1000, 1);
            
            % Initial deployment area
            S = pi * (self.vehicle.D_flotation_initial/2)^2;
            
            % Initial results
            V_termDescent = v(self, self.termDescentAlt);
            Drag_terminal = D_terminal(self, self.planet.rho(self.termDescentAlt), V_termDescent, self.vehicle.Cd_flotation, S);
            Drag_opening = D_initial(self, self.vehicle.Cx, Drag_terminal);
            
            a(1) = a_terminal(self, self.planet.g_o, Drag_opening, self.vehicle.m);
            h(1) = self.termDescentAlt;
            V(1) = V_termDescent;
            
            % Loop through rest of descent
            i = 1;
            deployedHalf = false;
            deployedFull = false;
            while h(i) > self.landAlt
                i = i+1;
                V(i) = v_termianl(self, V(i-1), a(i-1), self.TIME_STEP);
                h(i) = h_terminal(self, h(i-1), V(i), self.TIME_STEP);
                Drag = D_terminal(self, self.planet.rho(h(i)), V(i), self.vehicle.Cd_flotation, S);
                a(i) = a_terminal(self, self.planet.g_o, Drag, self.vehicle.m);
                
                % Deploy half of flotation device
                if h(i) < self.floatDeployAlt
                    Diameter = self.vehicle.D_flotation_initial + ...
                        (self.vehicle.D_flotation - self.vehicle.D_flotation_initial)/2;
                    S = pi * (Diameter/2)^2;
                    if ~deployedHalf
                        S = self.vehicle.Cx * S;
                        deployedHalf = true;
                    end
                    
                    % Deply rest of flotation device
                    if h(i) < self.landAlt + (self.floatDeployAlt - self.landAlt)/2
                        Diameter = self.vehicle.D_flotation;
                        S = pi * (Diameter/2)^2;
                        if ~deployedFull
                            S = self.vehicle.Cx * S;
                            deployedFull = true;
                        end
                    end
                end
            end
            
            %Trim off unused indices
            h = h(1:i-1);
            V = V(1:i-1);
            a = a(1:i-1);
        end
        
        function [h, drd] = calculateDRD(self)
            %calculateDRD get DRD values vs altitude
            fp = fplot(@(x) self.drd(x), [self.termDescentAlt self.h_e]);
            drd = fp.YData;
            drd = max(drd) - drd;
            h = fp.XData./1000;
        end
        
        function [h, q] = calculateq(self)
            %calculateq get heat flux values vs altitude
            fp = fplot(@(x) self.q(x), [self.termDescentAlt self.h_e]);
            q = fp.YData;
            h = fp.XData./1000;
        end
        
        function acceleration = a(self, h)
            %ACCELERATION acceleration of vehicle in g's
            acceleration = - (1 - v(self, h).^2./ ...
                       (self.planet.r(h).*self.planet.g_o))./ ...
                                  self.LtD;
        end
        
        function velocity = v(self, h)
            %VELOCITY velocity of vehicle (m/s)
            velocity = sqrt(self.planet.g_o .* self.planet.r(h)./ ...
                       (1 + self.planet.rho(h) .* self.LtD .* self.planet.r(h)./(2 .* self.vehicle.BC)));
        end
        
        function drd = drd(self, h)
            %DRD Down range distance of vehicle (m)
            drd = -self.planet.r(h)./2 .* self.LtD .* log(1 - v(self, h).^2./ ...
                                                 (self.planet.r(h).*self.planet.g_o));
        end
        
        function q = q(self, h)
            %q heat rate experienced by vehicle (W/cm^2)
            q = (self.planet.c .* sqrt(self.planet.rho(h) ./ (self.vehicle.D_Vehicle ./6)) .* v(self,h).^3) ./ 100.^2;
        end

        function Drag = D_terminal(self, rho, V, Cd, A)
            %D_terminal Drag during terminal descent
            Drag = Cd*A*rho*V^2 / 2;
        end

        function a = a_terminal(self, g, D_i, m)
            %a_terminal Acceleration during terminal descent
            a = g - D_i/m;
        end

        function velocity = v_termianl(self, v_i, a, dt)
            %v_terminal Velocity during terminal descent
            velocity = v_i + a*dt;
        end
        
        function altitude = h_terminal(self, h_i, V_i, dt)
            %h_terminal Altitude during terminal descent
            altitude = h_i - V_i*dt*cos(self.trim);
        end
        
        function Drag_initial = D_initial(self, Cx, D)
            %D_initial Initial drag when first opening flotation device
            Drag_initial = Cx*D;
        end
    end
end

