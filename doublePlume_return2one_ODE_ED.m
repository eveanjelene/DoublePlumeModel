function[dy_ds] = doublePlume_return2one_ODE_ED(z, y, const, atmos)


%% Set variables

    % Inner plume
    m_d_inner     = y(1);                                                  % Value of MFR dry gas phase from last step inner
    m_v_inner     = y(2);                                                  % Value of MFR wet vapor phase from last step inner
    psi_inner     = y(3);                                                  % Value of momentum from last step inner
    Q_inner       = y(4);                                                  % Value of energy from last step inner
    P             = y(5);                                                  % Value of pressure from last step inner
    m_class_inner = y(6:6 + const.no_class_inner - 1);                     % Value of MFR solid phase of each grainsize from last step inner
    m_s_inner     = sum(m_class_inner);                                    % Value of MFR solid phase of all grainsize from last step inner
    
    % Outer plume
    m_d_outer     = y(6 + const.no_class_inner);                           % Value of MFR dry gas phase from last step outer
    m_v_outer     = y(7 + const.no_class_inner);                           % Value of MFR wet vapor phase from last step outer
    psi_outer     = y(8 + const.no_class_inner);                           % Value of momentum from last step outer
    Q_outer       = y(9 + const.no_class_inner);                           % Value of energy from last step inner
    m_class_outer = y(10 + const.no_class_inner:(10 + const.no_class_inner) + (const.no_class_outer - 1)); %Value of mer solid phase of each grainsize from last step outer
    m_s_outer     = sum(m_class_outer);                                    % Value of MFR solid phase of all grainsize from last step outer


%% define constants
    theta_a0       = atmos.theta_a0;                                       % Initial temperature of atmospher at vent (k)
    H1             = atmos.H1;                                             % Height of the tropopause (m)
    H2             = atmos.H2;                                             % Height of the top of the stratospher (m)
    tempGrad_t     = atmos.tempGrad_t;                                     % Temperature gradient in the troposphere (k/m)
    tempGrad_s     = atmos.tempGrad_s;                                     % Temperature gradient in the startosphere (k/m)
    
    R_v            = const.R_v;                                            % Specifc gas constant of water vapor
    R_d            = const.R_d;                                            % Specifc gas constant of dry air
    C_s            = const.C_s;                                            % Specifc heat capacity (at constant pressure) solid
    C_d            = const.C_d;                                            % Specifc heat capacity (at constant pressure) dry air
    C_v            = const.C_v;                                            % Specifc heat capacity (at constant pressure) water vapour
    g              = const.g;                                              % Gravitational accelertion
    viscosity      = const.viscosity;                                      % Viscosity of air
    diameter_inner = const.diameter_inner;                                 % Vector of grain-sizes of inner plume
    no_class_inner = const.no_class_inner;                                 % No. grainsize classes of inner
    diameter_outer = const.diameter_outer;                                 % Vector of grain-sizes of outer plume
    no_class_outer = const.no_class_outer;                                 % No. grainsize classes of inner
    alpha          = const.alpha;                                          % Entrainment coefficant from outer to inner
    beta           = const.beta;                                           % Entrainment coefficant from inner to outer
    landa          = const.landa;                                          % Entrainment coefficant from ambient to outer

%% define atmospheric properties
    
% Temperature [K] from Wood et al., 1988
    if z <= H1
        theta_a = theta_a0 - tempGrad_t * z;
    elseif (z> H1) && (z<=H2)
        theta_a = theta_a0 - tempGrad_t * H1;
    else
        theta_a = theta_a0 - tempGrad_t * H1 + tempGrad_s * (z-H2);
    end
    
    eps    = R_d/R_v;
    w_a    = 0;
    rho_aB = P/(R_v*theta_a)*(1+w_a)/(w_a+eps);                            % Bulk density of the atmosphere
  
    rhophi_av = 0;                                                         % Density of water vapour phase in atmosphere set to zeros = no humidity so wet vapor phase of atmos P/(R_v*theta_a)*w_a/(w_a+eps);
    rhophi_ad = P/(R_v*theta_a)*1/(w_a+eps);                               % Density of dry gas phase in atmosphere 
    C_aB      = (C_d + w_a*C_v)/(1 + w_a);                                 % Specific heat capacity of the atmosphere

%% Determine consitutive equations of inner and outer plumes
    
    % Consitutive inner plume
    m_inner     = m_d_inner + m_v_inner + m_s_inner;                       % Total MFR in the inner region (all the phases)
    u_inner     = psi_inner/m_inner;                                       % Velocity of inner plume region (m/s)
    C_B_inner   = (m_d_inner*C_d + m_v_inner*C_v + m_s_inner*C_s)/m_inner; % Specific heat capcity of inner plume region
    theta_inner = (1/C_B_inner/m_inner)*Q_inner;                           % Temperature of inner plume
    n_d_inner   = m_d_inner/m_inner;                                       % Dry air MFR
    n_v_inner   = m_v_inner/m_inner;                                       % Water vapor MFR
    n_s_inner   = m_s_inner/m_inner;                                       % Solid MFR
    R_g_inner   = (n_v_inner/(n_v_inner + n_d_inner))*R_v + ...
                  (1 - n_v_inner/(n_v_inner + n_d_inner))*R_d;             % Specific gas constant in inner plume region
    rho_g_inner = P/R_g_inner/theta_inner;                                 % Density of gas phases in inner plume region
        
    if (m_v_inner <=0 )
        rho_B_inner = const.rho_s;     
        rhophi_d_inner = 0;                                                % Volume fraction dry gas
        rhophi_v_inner = 0;                                                % Volume fraction water vapour
    else    
        rho_B_inner = ((n_v_inner+n_d_inner)/rho_g_inner + (1-n_v_inner-n_d_inner)/const.rho_s)^(-1);     
        rhophi_d_inner = (P/R_v/theta_inner)*m_d_inner/(m_v_inner+eps*m_d_inner); % Volume fraction dry gas
        rhophi_v_inner = (P/R_v/theta_inner)*m_v_inner/(m_v_inner+eps*m_d_inner); % Volume fraction water vapour
    end
      
    r_inner = (m_inner/(rho_B_inner*u_inner))^(1/2);                       % Radius of the inner plume region (m)

    % Consitutive outer plume
    m_outer = m_d_outer + m_v_outer + m_s_outer;                           % Total MFR in the outer region (all the phases)
    u_outer = psi_outer/m_outer;                                           % Velocity of outer plume region
    C_B_outer = (m_d_outer*C_d + m_v_outer*C_v + m_s_outer*C_s)/m_outer;   % Specific heat capcity of outer plume region
    theta_outer = (1/C_B_outer/m_outer)*Q_outer;                           % Temperature of outer plume
    n_d_outer = m_d_outer/m_outer;                                         % Dry air MFR
    n_v_outer = m_v_outer/m_outer;                                         % Water vapor MFR
    n_s_outer = m_s_outer/m_outer;                                         % Solid MFR
    R_g_outer = (n_v_outer/(n_v_outer + n_d_outer))*R_v + ...
                (1 - n_v_outer/(n_v_outer + n_d_outer))*R_d;               % Specific gas constant in outer plume region
    rho_g_outer = P/R_g_outer/theta_outer;                                 % Density of gas phases in the outer plume region
    rho_B_outer = ((n_v_outer + n_d_outer)/rho_g_outer + (1 - n_v_outer-n_d_outer)/const.rho_s)^(-1); %Bulk density of the outer plume region
    
    rhophi_d_outer = (P/R_v/theta_outer)*m_d_outer/(m_v_outer + eps*m_d_outer); % Volume fraction dry gas
    rhophi_v_outer = (P/R_v/theta_outer)*m_v_outer/(m_v_outer + eps*m_d_outer); % Volume fraction water vapour
    r_outer = ((r_inner^2) + m_outer/(rho_B_outer*u_outer))^(1/2);         % Radius of the outer plume region (m)
        
    % Caluclate particle fallout inner
    [fallout_inner] = sedInner(u_inner,r_inner,diameter_inner,g,rho_B_inner,viscosity,const,m_class_inner);

    % Particle fallout outer
    i = 1;
    fallout_outer = zeros(no_class_outer,1);                               % Array to store what is lost from each grainsize in the plume
    while i <= no_class_outer
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
            falloutvelocity = lowREv; 
        elseif highRE >= 500
            falloutvelocity = highREv;
        else
            falloutvelocity = interREv;
        end

        dparticles_ds_outer = -const.prob*(m_class_outer(i)/r_outer)*(falloutvelocity/u_outer); %change in mass of particles of set size
        fallout_outer(i)    = dparticles_ds_outer;
        i = i +1;
    end
        
        
%% Entrainment velocities
    u_alpha = alpha * (u_inner-u_outer) * (rho_B_inner/rho_B_outer)^(0.5); % Velocity of the mass flux going from the outer to the inner
    u_beta = beta * u_outer * (rho_B_outer/rho_B_inner)^(0.5);             % Velocity of the mass flux going from the inner to the outer
    u_landa = landa * u_outer * (rho_B_outer/rho_aB)^(0.5);                % Atmospheric entrainment
      
%% Detmine how much is added to the outer plume (i.e. doesn't deposit straight out)
    [mass2add] = addFallout(u_outer,diameter_outer,g,rho_B_outer,viscosity,const,fallout_inner);
    
%% Governing equations
    
    %%%%% INNER EQUATIONS %%%%%
    %%%% Mass %%%%
    % Mass flux change of the inner region (dry air) 
    dm_d_ds_inner = 2*u_alpha*r_inner*rhophi_d_outer - 2*u_beta*r_inner*rhophi_d_inner; 
    
    % Mass flux change of the inner region (vapour air) 
    dm_v_ds_inner = 2*u_alpha*r_inner*rhophi_v_outer - 2*u_beta*r_inner*rhophi_v_inner; 
    
    % Mass flux change of the inner region (solid) 
    dm_s_ds_inner = sum(fallout_inner);                                    % Change to fallout
    
    % Total mass flux change of the inner region (all phases) 
    dm_ds_inner = dm_d_ds_inner + dm_v_ds_inner + dm_s_ds_inner;
    
    %%%% Momentum %%%%
    momentum_change_due_entrainment_inner = 2*u_alpha * r_inner * rhophi_d_outer * u_outer + ...
                                            2*u_alpha * r_inner * rhophi_v_outer * u_outer - ...
                                            2*u_beta * r_inner * rhophi_d_inner * u_inner - ...
                                            2*u_beta * r_inner * rhophi_v_inner * u_inner;
    
    % Momentum flux change of the inner region                                     
    dpsi_ds_inner = g*(rho_B_outer - rho_B_inner)*r_inner^2 + ...
                    momentum_change_due_entrainment_inner + ...
                    (u_inner*dm_s_ds_inner); 

    %%%% Energy %%%%    
    energy_change_due_entrainment_d_inner = 2*u_alpha*r_inner*rhophi_d_outer*C_d*theta_outer - ...
                                            2*u_beta*r_inner*rhophi_d_inner*C_d*theta_inner;
                                        
    energy_change_due_entrainment_v_inner = 2*u_alpha*r_inner*rhophi_v_outer*C_v*theta_outer - ...
                                            2*u_beta*r_inner*rhophi_v_inner*C_v*theta_inner;
                                        
    energy_change_due_entrainment_TOTAL_inner = energy_change_due_entrainment_d_inner + energy_change_due_entrainment_v_inner;
    
    % Energy flux change of the inner region
    dQ_ds_inner = energy_change_due_entrainment_TOTAL_inner - ...
                  rho_B_inner * u_inner * r_inner^2 *g * (rho_aB/rho_B_inner) + ...
                  C_s * theta_inner * dm_s_ds_inner; 
        
              
    %%%%% OUTER EQUATIONS %%%%%  
    %%%% Mass %%%%
    % Mass flux change of the outer region (dry air)
    dm_d_ds_outer = - 2*u_alpha*r_inner*rhophi_d_outer + ...
                    2*u_beta*r_inner*rhophi_d_inner + ...
                    2*u_landa*r_outer*rhophi_ad; 
    % Mass flux change of the outer region (vapour air)           
    dm_v_ds_outer = - 2*u_alpha*r_inner*rhophi_v_outer + ...
                    2*u_beta*r_inner*rhophi_v_inner + ...
                    2*u_landa*r_outer*rhophi_av;
                
    % Mass flux change of the outer region (solid)
    dm_s_ds_outer = sum(fallout_outer) - sum(mass2add);                    %change to fallout % minus inner as it is already a negative so need a negative to make it a positive
    
    % Total mass flux change of the outer region (all pahses)
    dm_ds_outer = dm_d_ds_outer + dm_v_ds_outer + dm_s_ds_outer;
    
    %%%% Momentum %%%%
    momentum_change_due_entrainment_outer = (-2*u_alpha*r_inner*rhophi_d_outer*u_outer) - ...
                                            (2*u_alpha*r_inner*rhophi_d_outer*u_outer) + ...
                                            (2*u_beta*r_inner*rhophi_d_inner*u_inner) + ...
                                            (2*u_beta*r_inner*rhophi_d_inner*u_inner);
    % Momentum flux change of the outer region                                      
    dpsi_ds_outer = g*(rho_aB - rho_B_outer)*(r_outer^2 - r_inner^2) + ...
                    momentum_change_due_entrainment_outer - ...
                    (u_inner*sum(mass2add)) + ...
                    (u_outer*sum(fallout_outer)); 
                
    %%%% Energy %%%%
    energy_change_due_entrainment_d_outer = - 2*u_alpha*r_inner*rhophi_d_outer*C_d*theta_outer + ...
                                              2*u_beta*r_inner*rhophi_d_inner*C_d*theta_inner + ...
                                              2*u_landa*r_outer*rhophi_ad*C_d*theta_a;
                                          
    energy_change_due_entrainment_v_outer = - 2*u_alpha*r_inner*rhophi_v_outer*C_v*theta_outer + ...
                                              2*u_beta*r_inner*rhophi_v_inner*C_v*theta_inner + ...
                                              2*u_landa*r_outer*rhophi_av*C_v*theta_a;
                                          
    energy_change_due_entrainment_TOTAL_outer = energy_change_due_entrainment_d_outer + energy_change_due_entrainment_v_outer;
    
    % Energy flux change of the outer region  
    dQ_ds_outer = energy_change_due_entrainment_TOTAL_outer - ...
                  rho_B_outer*u_outer*(r_outer^2-r_inner^2)*g*(rho_aB/rho_B_outer) + ...
                  (C_s*theta_outer*sum(fallout_outer)) - ...
                  C_s*theta_inner*sum(mass2add);
              
              
%               figure(1)
%               subplot(1,2,1)
%               scatter(dm_d_ds_outer,z,'b.')
%               hold on
%               
%               subplot(1,2,2)
%               scatter(dm_v_ds_outer,z,'b.')
%               hold on
%               subplot(1,4,3)
%               scatter(-rho_B_outer*u_outer*(r_outer-r_inner)^2*g*(rho_aB/rho_B_outer),z,'b.')
%               hold on
%               
%               subplot(1,4,4)
%               scatter(theta_outer,z,'b.')
%               hold on
              
%% assign change
    dy_ds = zeros((10 + const.no_class_inner) + (const.no_class_outer-1),1); %Create array to store rate of change for each variable

    dP_ds = -rho_aB*g;                                                     %Calculate the rate in change in pressure
    
    dy_ds(1)                             = dm_d_ds_inner;
    dy_ds(2)                             = dm_v_ds_inner;
    dy_ds(3)                             = dpsi_ds_inner;
    dy_ds(4)                             = dQ_ds_inner;
    dy_ds(5)                             = dP_ds;
    dy_ds(6:6 + const.no_class_inner -1) = fallout_inner(1:no_class_inner);
    dy_ds(6 + const.no_class_inner)      = dm_d_ds_outer;
    dy_ds(7 + const.no_class_inner)      = dm_v_ds_outer;
    dy_ds(8 + const.no_class_inner)      = dpsi_ds_outer;
    dy_ds(9 + const.no_class_inner)      = dQ_ds_outer;
    dy_ds(10 + const.no_class_inner:(10 + const.no_class_inner) + (const.no_class_outer - 1)) = fallout_outer(1:no_class_inner);
    
end

function[fallout_inner] = sedInner(u_inner,r_inner,diameter_inner,g,rho_B_inner,viscosity,const,m_class_inner)
    
    %set variables
    u_th = u_inner;                                                        %top hat velocity
    R_th = r_inner;                                                        %top hat radius
    
    b = R_th/sqrt(2);
    
    fallout_inner = zeros(const.no_class_inner,1);                         %Create output array for fallout of each size
    
    for i = 1:const.no_class_inner 
        
        if m_class_inner(i) == 0
            fallout_inner(i) = 0;
        else 
            %calculate settling velocity
                % Calculate potential fallout velocity in each RE regieme
                lowREv   = ((diameter_inner(i)^2)*g*(2000 - rho_B_inner))/(18*viscosity);
                interREv = diameter_inner(i)*((4*(g^2)*((2000 - rho_B_inner)^2)/(225*viscosity*rho_B_inner))^(1/3));
                highREv  = ((3.1*diameter_inner(i)*g*(2000 - rho_B_inner)/rho_B_inner))^0.5;

                % Calculate Re using tyhe velocities calculated for each RE regieme
                lowRE   = (diameter_inner(i)*lowREv*rho_B_inner) / viscosity;
                interRE = (diameter_inner(i)*interREv*rho_B_inner) / viscosity;
                highRE  = (diameter_inner(i)*highREv*rho_B_inner) / viscosity;

                % Determine which velocity should be used based on the RE number
                if  lowRE > 0 && lowRE <= 6
                    u_ti = lowREv; 
                elseif highRE >= 500
                    u_ti = highREv;
                else
                    u_ti = interREv;
                end
%                 
%                 
%                 if i == 1
%                     figure(3)
%                     hold on
%                     scatter(2*u_th,z,'b.')
% 
%                     scatter(u_ti,z,'r.')
% 
% %                     figure(4)
% %                     hold on
% %                     scatter(Re,z,'r.')
%                  end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if 2*u_th > (u_ti)%+0.1)                                               %if we have an envolope folloew through with this methos
                %work out Renv
                R_env = b * (log(((u_th*2)/u_ti)))^(0.5);

                %Calulate areas under the curve
                y_r = R_env/R_th;
                
                %Work out the ratio between the are under the whole cure and area to left of envelope
                area_ratio = 1 - erf(sqrt(2)*y_r);                          

                %Work out the amount of mass lost
                dparticles_ds_inner = -(m_class_inner(i)/r_inner)*area_ratio; %change in mass of particles of set size
                fallout_inner(i)    = dparticles_ds_inner;
            else
                fallout_inner(i)    = 0;% -(m_class_inner(i));                 %If we dont have an envolope becuse the particle is settling faster at all points, loose all mass
                
%                 if i == 1
%                     disp(190)
%                 end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
        end
    end

end