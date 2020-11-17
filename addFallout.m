function [mass2add] = addFallout(u_outer,diameter_outer,g,rho_B_outer,viscosity,const,fallout_inner)

    %set variables
    u_th = u_outer;                                                        % Top hat velocity                                                      %top hat radius
    
    mass2add = zeros(const.no_class_inner,1);                              % Create output array for fallout of each size
    
    for i = 1:const.no_class_inner 
        
            % Calculate settling velocity
                % Calculate potential fallout velocity in each RE regieme
                lowREv   = ((diameter_outer(i)^2)*g*(2000 - rho_B_outer))/(18*viscosity);
                interREv = diameter_outer(i)*((4*(g^2)*((2000 - rho_B_outer)^2)/(225*viscosity*rho_B_outer))^(1/3));
                highREv  = ((3.1*diameter_outer(i)*g*(2000 - rho_B_outer)/rho_B_outer))^0.5;

                % Calculate Re using tyhe velocities calculated for each RE regieme
                lowRE   = (diameter_outer(i)*lowREv*rho_B_outer) / viscosity;
                interRE = (diameter_outer(i)*interREv*rho_B_outer) / viscosity;
                highRE  = (diameter_outer(i)*highREv*rho_B_outer) / viscosity;

                % Determine which velocity should be used based on the RE number
                if  lowRE > 0 && lowRE <= 6
                    u_ti = lowREv; 
                elseif highRE >= 500
                    u_ti = highREv;
                else
                    u_ti = interREv;
                end

               
            if (2*u_th) > (u_ti)                                           % If we have an envolope folloew through with this method
                % Supported so added to the outer
                mass2add(i)    = fallout_inner(i); 
            else
                % Not supported so not added to outer
                mass2add(i)    = 0;                                        %If we dont have an envolope becuse the particle is settling faster at all points, loose all mass
            end

    end
       
end

