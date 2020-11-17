function [value,isterminal,direction] = doublePlume_return2one_stopInner(z,y,const,atmos)

    m_d_inner     = y(1);                                                  %Value of mer dry gas phase from last step inner
    m_v_inner     = y(2);                                                  %Value of mer wet vapor phase from last step inner
    psi_inner     = y(3);                                                  %Value of momentum from last step inner
    Q_inner       = y(4);                                                  %Value of energy from last step inner
    P             = y(5);                                                  %Value of pressure from last step inner
    m_class_inner = y(6:6 + const.no_class_inner - 1);                     %Value of mer solid phase of each grainsize from last step inner
    m_s_inner     = sum(m_class_inner);                                    %Value of mer solid phase of all grainsize from last step inner

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
    
%% atmosphere
    % Temperature [K] from Wood et al., 1988
    if z <= H1
        theta_a=theta_a0-tempGrad_t*z;
    elseif (z> H1) && (z<=H2)
        theta_a=theta_a0-tempGrad_t*H1;
    else
        theta_a=theta_a0-tempGrad_t*H1+tempGrad_s*(z-H2);
    end
    
    eps    = R_d/R_v;
    w_a    = 0;
    rho_aB = P/(R_v*theta_a)*(1+w_a)/(w_a+eps);                            %bulk density of the atmosphere
  
%% consititive equations
    %outer
    m_inner = m_d_inner + m_v_inner + m_s_inner;                           %Total mass flux in the inner region (all the phases)
    u_inner = psi_inner/m_inner;                                           %Velocity of inner plume region

    u       = u_inner;
    
    C_B_inner   = (m_d_inner*C_d + m_v_inner*C_v + m_s_inner*C_s)/m_inner; %Specific heat capcity of inner plume region
    theta_inner = (1/C_B_inner/m_inner)*Q_inner;                           %Temperature of inner plume
    n_d_inner   = m_d_inner/m_inner;                                       %Dry air mass flux
    n_v_inner   = m_v_inner/m_inner;                                       %Water vapor mass flux
    n_s_inner   = m_s_inner/m_inner;                                       %Solid mass flux
    R_g_inner   = (n_v_inner/(n_v_inner + n_d_inner))*R_v + ...
                  (1 - n_v_inner/(n_v_inner + n_d_inner))*R_d;             %Specific gas constant in inner plume region
    rho_g_inner = P/R_g_inner/theta_inner;                                 %Density of gas phases in inner plume region
        
    % Chocolate (to be checked)
    if (m_v_inner <=0 )
        rho_B_inner = const.rho_s;                                                   %Volume fraction water vapour
    else    
        rho_B_inner = ((n_v_inner+n_d_inner)/rho_g_inner + (1-n_v_inner-n_d_inner)/const.rho_s)^(-1);     
    end
      
%     if m_s_inner <= 10
%         u = 0.1;
%          disp('stop as no mass')


%%
u_ti = zeros(size(const.diameter_inner));
    for i = 1:const.no_class_inner 
%    if 
        diameter_inner = const.diameter_inner;
        viscosity = const.viscosity;
        g = const.g;
        %i = const.run;
        if m_class_inner(i) == 0
            u_ti(i) = 0;
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
                    u_ti(i) = lowREv; 
                elseif highRE >= 500
                    u_ti(i) = highREv;
                else
                    u_ti(i) = interREv;
                end
        end


    end
%% Locate the height when velocity passes through zero
value = [u-0.1,m_s_inner-1,(u_inner*2)-u_ti];%-0.1];                                                          %Detect velocity = 0
                                                                           %the 0.001 helps for stability
isterminal = [1,1,zeros(size(const.diameter_inner))+1];                                                          %Stop the integration
direction  = [0,0,zeros(size(const.diameter_inner))];                                                          %Detect all zeros (default)
