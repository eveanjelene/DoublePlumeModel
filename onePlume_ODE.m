function[dy_ds] = onePlume_ODE(z,y,const,atmos)


%% set variables
    %inner plume
    m_d_one       = y(1);                                                  %Value of mer dry gas phase from last step
    m_v_one       = y(2);                                                  %Value of mer wet vapour phase from last step 
    psi_one       = y(3);                                                  %Value of mpmentume from last step 
    Q_one         = y(4);                                                  %Value of energy from last step 
    P             = y(5);                                                  %Value of pressure from last step inner
    m_class_one   = y(6:6 + const.no_class_inner - 1);                     %Value of mer solid phase of each grainsize from last step 
    m_s_one= sum(m_class_one);                                             %Value of mer solid phase of all grainsize from last step 
    
%     if m_class_one(1) <= 0
%         disp(m_class_one(1))
%     end
%     m_class_one(m_class_one <= 0) = 0;
% %     m_v_inner(m_v_inner <= 0) = 0;
% % %         m_class_inner(m_class_inner <= 0.001) = 0;
% %      if m_v_inner < 0.0
% %         disp('error')
% %      end
%% define constants
    theta_a0       = atmos.theta_a0;                                       %Initial temperature of atmospher at vent (k)
    H1             = atmos.H1;                                             %Height of the tropopause (m)
    H2             = atmos.H2;                                             %Height of the top of the stratospher (m)
    tempGrad_t     = atmos.tempGrad_t;                                     %Temperature gradient in the troposphere (k/m)
    tempGrad_s     = atmos.tempGrad_s;                                     %Temperature gradient in the startosphere (k/m)
    
    R_v            = const.R_v;                                            %Specifc gas constant of water vapor
    R_d            = const.R_d;                                            %Specifc gas constant of dry air
    C_s            = const.C_s;                                            %Specifc heat capacity (at constant pressure) solid
    C_d            = const.C_d;                                            %Specifc heat capacity (at constant pressure) dry air
    C_v            = const.C_v;                                            %Specifc heat capacity (at constant pressure) water vapour
    g              = const.g;                                              %Gravitational accelertion
    viscosity      = const.viscosity;                                      %Viscosity of air
    diameter_inner = const.diameter_inner;                                 %Power coeffiacnt of GSD for inner
    no_class_inner = const.no_class_inner;                                 %No. grainsize classes of inner
    diameter_outer = const.diameter_outer;                                 %Power coeffiacnt of GSD for outer
    no_class_outer = const.no_class_outer;                                 %No. grainsize classes of inner
    alpha          = const.alpha;                                          %Entrainment coefficant from outer to inner
    beta           = const.beta;                                           %Entrainment coefficant from inner to outer
    landa          = const.landa;                                          %Entrainment coefficant from ambient to outer
    
    u_ambient      = 0;                                                    %Velocity of the atmosphere
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
    rho_aB = P/(R_v*theta_a)*(1+w_a)/(w_a+eps);                            %bulk density of the atmosphere
  
    rhophi_av = 0;                                                         %Density of water vapour phase in atmosphere set to zeros = no humidity so wet vapor phase of atmos P/(R_v*theta_a)*w_a/(w_a+eps);
    rhophi_ad = P/(R_v*theta_a)*1/(w_a+eps);                               %Density of dry gas phase in atmosphere 
    C_aB      = (C_d + w_a*C_v)/(1 + w_a);                                 %Specific heat capacity of the atmosphere
  
%% determine consitutive equations if inner and outer plumes
    
    % Consitutive plume
    m_one         = m_d_one + m_v_one + m_s_one;                           %Total mass flux in the plume (all the phases)
    u_one         = psi_one/m_one;                                         %Velocity of plume (m/s)
    C_B_one       = (m_d_one*C_d + m_v_one*C_v + m_s_one*C_s)/m_one;       %Specific heat capcity of plume
    theta_one     = (1/C_B_one)*(Q_one/m_one);                             %Temperature of inner
    n_d_one       = m_d_one/m_one;                                         %Dry air mass flux
    n_v_one       = m_v_one/m_one;                                         %Water vapor mass flux
    n_s_one       = m_s_one/m_one;                                         %Solid mass flux
    R_g_one       = (n_v_one/(n_v_one + n_d_one))*R_v + (1 - n_v_one/ ...
                    (n_v_one + n_d_one))*R_d;                              %Specific gas constant in plume
    rho_g_one     = P/R_g_one/theta_one;                                   %Density of gas phases in plume
    rho_B_one     = ((n_v_one+n_d_one)/rho_g_one + (1-n_v_one-n_d_one)/const.rho_s)^(-1);
    %%%%%%CHECK%%%%%%%%%%
    rhophi_d_one  = (P/R_v/theta_one)*m_d_one/(m_v_one+eps*m_d_one);       %Volume fraction dry gas
    rhophi_v_one  = (P/R_v/theta_one)*m_v_one/(m_v_one+eps*m_d_one);       %Volume fraction water vapour
    %%%%%%%%%%%%%%%%%%%%%
    r_one         = (m_one/(rho_B_one*u_one))^(1/2);                       %Radius of the plume (m)

    % Particle fallout 
    i = 1;
    fallout_one = zeros(size(no_class_outer));                             %Array to store what is lost from each grainsize in the plume
    while i <= no_class_outer
        % Calculate potential fallout velocity in each RE regieme
        lowREv   = ((diameter_outer(i)^2)*g*(2000 - rho_B_one))/(18*viscosity);
        interREv = diameter_outer(i)*((4*(g^2)*((2000 - rho_B_one)^2)/(225*viscosity*rho_B_one))^(1/3));
        highREv  = ((3.1*diameter_outer(i)*g*(2000 - rho_B_one)/rho_B_one))^0.5;

        % Calculate Re using tyhe velocities calculated for each RE regieme
        lowRE   = (diameter_outer(i)*lowREv*rho_B_one) / viscosity;
        interRE = (diameter_outer(i)*interREv*rho_B_one) / viscosity;
        highRE  = (diameter_outer(i)*highREv*rho_B_one) / viscosity;

        % Determine which velocity should be used based on the RE number
        if  lowRE > 0 && lowRE <= 6
            falloutvelocity = lowREv; 
        elseif highRE >= 500
            falloutvelocity = highREv;
        else
            falloutvelocity = interREv;
        end
        
        dparticles_ds_one = -const.prob*(m_class_one(i)/r_one)*(falloutvelocity/u_one);  %change in mass of particles of givern size i.e. the amount that has fallen out
        fallout_one(i)    = dparticles_ds_one;
        i = i + 1;
    end
        
%% Entrainment velocities
        u_beta = 0;                                                        %Velocity of the mass flux going from the inner to the outer
        u_alpha = 0;                                                       %Velocity of the mass flux going from the outer to the inner
        u_landa = 0.1 * u_one * (rho_B_one/rho_aB)^(0.5);                  %Atmospheric entrainment for a one plume based on Morton et al., 1956 with additional bousinesq approximation
                
%% Govening equations
    %u_beta is set to 0 so all thoose terms will cancel out. Plume is no wind, momentum term also cancels out so no mentuum added to the plume from the entrainment of air
        
    %%%% Mass %%%%
    % Mass flux change of the plume (dry air) 
    dm_d_ds_one   = 2*u_landa*r_one*rhophi_ad; % - 2*u_beta*r_inner*rhophi_ad; 
    
    % Mass flux change of the plume (vapour air) 
    dm_v_ds_one   = 2*u_landa*r_one*rhophi_av; %- 2*u_beta*r_inner*rhophi_av; 
    
    % Mass flux change of the plume (solid) 
    dm_s_ds_one   = sum(fallout_one); %change to fallout
    
    % Total mass flux change of the plume (all phases) 
    dm_ds_one     = dm_d_ds_one + dm_v_ds_one + dm_s_ds_one;
    
    %%%% Momentum %%%%
    momentum_change_due_entrainment_one = 2*u_landa*r_one*rhophi_ad*u_ambient + ...
                                          2*u_landa*r_one*rhophi_av * u_ambient;% -2*u_beta*r_inner*rhophi_d_inner*u_inner - 2*u_beta*r_inner*rhophi_v_inner*u_inner; %should be 0 if no wind
    
    % Momentum flux change of the plume                               
    dpsi_ds_one   = g*(rho_aB - rho_B_one)*r_one^2 + ...
                    momentum_change_due_entrainment_one + ...
                    (u_one*dm_s_ds_one); 
    
    %%%% Energy %%%% 
    energy_change_due_entrainment_d_one = 2*u_landa*r_one*rhophi_ad*C_d*theta_a;% - 2*u_beta*r_inner*rhophi_d_inner*C_d*theta_inner;
    energy_change_due_entrainment_v_one = 2*u_landa*r_one*rhophi_av*C_v*theta_a;% - 2*u_beta*r_inner*rhophi_v_inner*C_v*theta_inner;
    energy_change_due_entrainment_TOTAL_one = energy_change_due_entrainment_d_one + energy_change_due_entrainment_v_one;
    
    % Energy flux change of the plume
    dQ_ds_one = energy_change_due_entrainment_TOTAL_one - ...
                rho_B_one*u_one*r_one^2*g*(rho_aB/rho_B_one) + ...
                C_s*theta_one*dm_s_ds_one; 

%             figure(2)
%             hold on
%             scatter((u_one*dm_s_ds_one),z,'k.')
%              scatter(momentum_change_due_entrainment_one,z,'r.')
%              
%              figure(6)
%              subplot(1,2,1)
%              hold on 
%              scatter(u_one,z,'k.')
%              subplot(1,2,2)
%              hold on
%              scatter(g*(rho_aB - rho_B_one)*r_one^2,z,'b.')
%              
%              figure(3)
%             hold on
%             scatter(dm_s_ds_one,z,'k.')
%              scatter(dm_d_ds_one,z,'r.')
             
            
%% assign change
    dy_ds    = zeros((5 + const.no_class_inner-1),1);                      %Create array to store rate of change for each variable
    
    dP_ds    = -rho_aB*g;                                                  %Calculate the rate in change in pressure
    
    dy_ds(1)                              = dm_d_ds_one;
    dy_ds(2)                              = dm_v_ds_one;
    dy_ds(3)                              = dpsi_ds_one;
    dy_ds(4)                              = dQ_ds_one;
    dy_ds(5)                              = dP_ds;
    dy_ds(6:6 + const.no_class_inner - 1) = fallout_one; 
        
end