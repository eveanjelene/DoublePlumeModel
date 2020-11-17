function [value,isterminal,direction] = onePlume_stop(z,y,const,atmos,outerPlume)

%% set variables
    %inner plume
    m_d_one       = y(1);                                                  %Value of mer dry gas phase from last step
    m_v_one       = y(2);                                                  %Value of mer wet vapour phase from last step 
    psi_one       = y(3);                                                  %Value of mpmentume from last step 
    Q_one         = y(4);                                                  %Value of energy from last step 
    P             = y(5);                                                  %Value of pressure from last step inner
    m_class_one   = y(6:6 + const.no_class_inner - 1);                     %Value of mer solid phase of each grainsize from last step 
    m_s_one= sum(m_class_one);  
    
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
    rho_aB = P/(R_v*theta_a)*(1 + w_a)/(w_a + eps);                        %Bulk density of the atmosphere
   
%% consititive equations
    %outer
    m_one = m_d_one + m_v_one + m_s_one;                                   %Total mass flux in the plume (all the phases)
    u_one = psi_one/m_one;                                                 %Velocity of plume                        

    u     = u_one;
    
%% Locate the height when velocity passes through zero
value = [u-0.01];                                                          %Detect velocity = 0
                                                                           %the 0.001 helps for stability
isterminal = [1];                                                          %Stop the integration
direction  = [0];                                                          %Detect all zeros (default)
