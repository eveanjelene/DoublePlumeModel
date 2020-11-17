function[output] = double_PM(MER,input_i,input_o,addMass,gs_size,ratio)
%%%%%%%%%% inner Double plume model %%%%%%%%%%%%%%%%%%%%%%
%   - no wind
%   - bousinesq approximation used from Devenish et al 2018
%   - use rho_outer - rho_inner foe buoyancy terms of the inner plume
%   - no entrainment from inner to outer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up

% Model constants
    const.g         = 9.81;                                                % Gravitational accelertion
    const.R_v       = 461;                                                 % Specifc gas constant of water vapor
    const.R_d       = 287;                                                 % Specifc gas constant of dry air
    const.C_d       = 998;                                                 % Specifc heat capacity (at constant pressure) dry air
    const.C_v       = 1952;                                                % specifc heat capacity (at constant pressure) water vapor
    const.C_s       = 1250;                                                % specifc heat capacity (at constant pressure) solid
    const.L         = 2.257e6;                                             % Latent heat of vaporization
    const.rho_l     = 1000;                                                % Density liquid water in plume
    const.rho_s     = 2000;                                                % Density solids in plume
    const.prob      = 0.27;                                                % Probability of fallout
    const.viscosity = 1.81*(10^-5);                                        % Viscosity of air
    const.alpha     = 0;                                                   % Entrainment coefficant from outer to inner 
    const.beta      = 0.147;                                               % Entrainment coefficant from inner to outer
    const.landa     = 0.147;                                               % Entrainment coefficant from ambient to outer
    const.r_h       = 0;                                                   % set to no humidity

% Atmospheric constants    
    atmos.H1           = 11000;                                            % Tropopause height (m)
    atmos.H2           = 20000;                                            % Stratosphere top height (m)
    atmos.tempGrad_t   = 6.5/1000;                                         % Temperature gradient in the troposphere (k/m)
    atmos.tempGrad_s   = 2/1000;                                           % Temperature gradient in the startosphere (k/m)
    atmos.theta_a0     = 280;                                              % Temperature at the vent (k)
 
%% Initial conditions
    
    P_0 = 101325;                                                          % Initial pressure at vent (Pa)
    z_0 = 0;                                                               % Vent height (m asl)    
    
    inital_MER0 = MER/pi;                                                  % Initial MFR of both inner and outer
    initial_MER0_inner = inital_MER0 * (ratio);                            % Initial MFR of inner
    initial_MER0_outer = inital_MER0 * (1-ratio);                          % Initial MFR of outer
    
%% Inner plume initial conditions
    t_0_inner = input_i(2);                                                % Initial plume temp (k)
    u_0_inner = input_i(1);                                                % Initial plume velocity (m/s)
    n_0_inner = input_i(3);                                                % Initial gas mass fraction

    rho_v0_inner = P_0/const.R_v/t_0_inner;                                % Initial density of the water vapour phase
    rho_b0_inner = (n_0_inner/rho_v0_inner + ...
                       (1 - n_0_inner)/const.rho_s)^(-1);                  % Initial bulk density
        
    % Set up mass distributions 
    m_0_inner     = initial_MER0_inner;                                    % Initial MFR of inner (kg/s)
    r_0_inner     = sqrt(m_0_inner/(rho_b0_inner*u_0_inner));              % Initial radius of the inner plume (m)
    m_s0_inner    = m_0_inner * (1-n_0_inner);                             % Initial inner MER of the solid fraction
    m_d0_inner    = 0;                                                     % Initial inner MER of the dry gas fraction
    m_v0_inner    = n_0_inner * m_0_inner;                                 % Initial inner MER of the wet vapour fraction

 
    % Create grainsize array and mass fration array which is the same size as the outer plume
    diameterphi_inner   = -9:0.5:10;
    mass_fraction_inner = zeros(size(diameterphi_inner));
            
    one_gs = 0;                                                            % Set to 0 if multiple grain-size used
    x_phi = -9:0.5:10;
    
    % Create gaussian distribution
    [mass_fraction_i] = gaussDistribution(gs_size,0.75,x_phi);            % Inputs - median grain-size, standard deviation and grain-sizes
    
    % Correct gaassuian of extends out of given limits
    if sum(mass_fraction_i) <1
        x_phi = -8.5:0.5:10;
        [mass_fraction_i] = gaussDistribution(gs_size,0.75,x_phi);
        mass_fraction_inner(2:end) = mass_fraction_i;
        mass_fraction_inner(1) = 1 - sum(mass_fraction_inner);
    else
        mass_fraction_inner = mass_fraction_i;
    end
    
    mass_fraction_inner = mass_fraction_inner./sum(mass_fraction_inner);   % Normalise 


    diameter_inner       = (2.^(-diameterphi_inner))./1000;                % Grain-size in (m)
    no_class_inner       = length(diameter_inner);                         % Number of grain-size
    m_class_0_inner      = mass_fraction_inner .* m_s0_inner;              % MFR of each grain-size  

    % assign to the const 
    const.no_class_inner = no_class_inner;
    const.diameter_inner = diameter_inner;
            
    C_B0_inner  = (m_d0_inner*const.C_d + m_v0_inner*const.C_v + m_s0_inner * const.C_s)...
                  /m_0_inner;                                              % Initial bulk heat capcity
    psi_0_inner = m_0_inner * u_0_inner;                                   % Initial momentum
    Q_0_inner   = m_0_inner * C_B0_inner*t_0_inner;                        % Initial energy

%% Outer plume initial conditions
    t_0_outer = input_o(2);                                                % Initial plume temp (k)
    u_0_outer = input_o(1);                                                % Initial plume velocity (m/s)
    n_0_outer = input_o(3);                                                % Initial gas mass fraction

    rho_v0_outer = P_0/const.R_v/t_0_outer;                                % Initial density of the gas phase outer plume
    rho_b0_outer = (n_0_outer/rho_v0_outer + (1-n_0_outer)/const.rho_s)^(-1); % Initial bulk density

    % Set up mass distributions of outer plume
    m_0_outer  = initial_MER0_outer;                                       % Initial MER (kg/s)
    r_0_outer  = sqrt((m_0_outer/(u_0_outer*rho_b0_outer)) + r_0_inner^2); % Initial radius of outer plume (m)
    m_s0_outer = m_0_outer * (1-n_0_outer);                                % Initial MFR of the solid fraction
    m_d0_outer = 0;                                                        % Initial MFR of the dry gas fraction
    m_v0_outer = n_0_outer * m_0_outer;                                    % Initial MFR of the wet gas fraction
    
    % Create outer plume grain-size distribution
    x_phi = -9:0.5:10;
    [mass_fraction_outer] = gaussDistribution(0.5,1.5,x_phi);              % calculate gaussian grain-size distrinution
    mass_fraction_outer = mass_fraction_outer./sum(mass_fraction_outer);   % Normalise gaussian
    diameterphi_outer    = -9:0.5:10;                                      % Diameter of grainsize in phi
    diameter_outer       = (2.^(-diameterphi_outer))./1000;                % Diameter of grainsize in (m)
    no_class_outer       = length(diameter_outer);                         % Number of GS
    m_class_0_outer      = mass_fraction_outer .* m_s0_outer;              % MFR in each grainsize  

    % Set to const
    const.no_class_outer = no_class_outer;
    const.diameter_outer = diameter_outer;

    C_B0_outer  = (m_d0_outer*const.C_d + m_v0_outer*const.C_v + m_s0_outer*const.C_s)...
                  /m_0_outer;                                              % Initial bulk heat capcity
    psi_0_outer = m_0_outer * u_0_outer;                                   % Initial momentum
    Q_0_outer   = m_0_outer * C_B0_outer * t_0_outer;                      % Initial energy


%% Model the double plume region
    lastwarn('');
    
    options = odeset('RelTol',1e-11,'AbsTol',1e-11,'Events',@(z,y) doublePlume_return2one_stopInner(z,y,const,atmos));                     
    Sspan   = [z_0 40000];                                                 % Solve from z=0m to z=40000m

    IC      = [m_d0_inner m_v0_inner psi_0_inner Q_0_inner P_0 ...         % Initial conditions for ode solver
               m_class_0_inner(1:const.no_class_inner) m_d0_outer ...
               m_v0_outer psi_0_outer Q_0_outer ... 
               m_class_0_outer(1:const.no_class_outer)];                   
    
    % Solve the double plume region
    odeSolution = ode15s(@(z,y) doublePlume_return2one_ODE_ED(z,y,const,atmos),Sspan,IC,options);

    
    % Model output
    Y        = odeSolution.y';                                             % Output matrix from ode solver
    
    % Inner ODEs output i.e. respective paremter values in the plume with z
    m_d_inner_1     = Y(:,1);                                              % Ode inner output MFR dry gas phase
    m_v_inner_1     = Y(:,2);                                              % Ode inner output MFR wet gas phase
    psi_inner_1     = Y(:,3);                                              % Ode inner output momentum
    Q_inner_1       = Y(:,4);                                              % Ode inner output energy
    P_double_1      = Y(:,5);                                              % Ode output pressure
    m_class_inner_1 = Y(:,6:6 + const.no_class_inner - 1);                 % Ode inner output MFR of each grainsize in the solid phase
    m_s_inner_1     = sum(m_class_inner_1,2);                              % Ode inner output MFR of solid phase (i.e. sum of all grainsize)
    m_inner_1       = m_d_inner_1 + m_v_inner_1 + m_s_inner_1;             % Ode inner output total MFR in plume of all phases
    
    % Outer ODEs output
    m_d_outer_1     = Y(:,6 + const.no_class_inner);                       % Ode outer output MFR dry gas phase
    m_v_outer_1     = Y(:,7 + const.no_class_inner);                       % Ode outer output MFR wet vapour phase
    psi_outer_1     = Y(:,8 + const.no_class_inner);                       % Ode outer output momentum
    Q_outer_1       = Y(:,9 + const.no_class_inner);                       % Ode outer output energy
    m_class_outer_1 = Y(:,10 + const.no_class_inner:(10 + const.no_class_inner) + (const.no_class_outer - 1)); % Ode outer output MFR of each grainsize in the solid phase
    m_s_outer_1     = sum(m_class_outer_1,2);                              % Ode outer output mer of solid phase (i.e. sum of all grainsize)
    m_outer_1       = m_d_outer_1 + m_v_outer_1 + m_s_outer_1;             % Ode outer output total MFR in plume of all phases
    
    
    % Check if given warning
    [warnmsg, msgid] = lastwarn;
        
    % If model fails, return to main script 
    if ~isempty(warnmsg) 
       disp('fail')                                            
       output.success = 0;
       return
    end
     
    % If model doesn't fail set output.success to 1
    output.success = 1;
    
    % If ODE solver stops because of inner plume velocity reaches close to
    % 0, or total MFR is close to 0, confirm output variables. Else remove MFR of grain-size that
    % has fallen out and then re-run model from last height poiunt.
    if odeSolution.ie == 1 || odeSolution.ie == 2 || one_gs == 1
            m_d_inner     = m_d_inner_1;                                   % ODE inner output MFR dry gas phase
            m_v_inner     = m_v_inner_1;                                   % ODE inner output MRF wet gas phase
            psi_inner     = psi_inner_1;                                   % ODE inner output momentum
            Q_inner       = Q_inner_1;                                     % ODE inner output energy
            P_double      = P_double_1;                                    % ODE output pressure
            m_class_inner = m_class_inner_1;                               % ODE inner output MFR of each grainsize in the solid phase
            m_s_inner     = sum(m_class_inner,2);                          % ODE inner output MFR of solid phase (i.e. sum of all grainsize)
            m_inner       = m_d_inner + m_v_inner + m_s_inner;             % ODE inner output total MFR in plume of all phases

            % Outer ODEs output
            m_d_outer     = m_d_outer_1;                                   % ODE outer output mer dry gas phase
            m_v_outer     = m_v_outer_1;                                   % ODE outer output mer wet vapour phase
            psi_outer     = psi_outer_1;                                   % ODE outer output momentum
            Q_outer       = Q_outer_1;                                     % ODE outer output energy
            m_class_outer = m_class_outer_1;                               % ODE outer output mer of each grainsize in the solid phase
            m_s_outer     = sum(m_class_outer,2);                          % ODE outer output MFR of solid phase (i.e. sum of all grainsize)
            m_outer       = m_d_outer + m_v_outer + m_s_outer;             % ODE outer output total MFR in plume of all phases
            
            z_double      = odeSolution.x';                                % z array used in the ODE solver for the above matrix
            
    else %Else remove MFR of grain-size that has fallen out and then re-run model from last height point.
        
        % Create storage cell to store model output in each time model is restarted
        storage = cell(1,const.no_class_inner);           
        
        counter = 1;
        
        % Store model output from 1st run of the ODE solver
        storage{1} = odeSolution;
        
        % Cycle through 1 to number of grain-size as model can only stop
        % this way as the same number of grain-sizes used
        for jj = 1:const.no_class_inner
            
            % f=For first cucle assign ODE output deom the storage cell
            if jj ~= 1
                odeSolution = storage{counter};
            end
            
            
            % If stosp because no more envelope, remove MFR related to that grain-size and then re-run from last point
            if odeSolution.ie >= 3                                      
                
                options = odeset('RelTol',1e-11,'AbsTol',1e-11,'Events',@(z,y) doublePlume_return2one_stopInner(z,y,const,atmos));                     
                Sspan   = [odeSolution.x(end) 40000];                      % Solve from z=0m to z=40000m
                
                % Remove MFR related to grain-size that is no longer supported
                [newSource] = add2Surrounding(odeSolution,const);

                % Assign new start of inner plume
                psi_0_inner     = newSource.psi_0_inner;
                Q_0_inner       = newSource.Q_0_inner;
                m_class_0_inner = newSource.m_class_0_inner;
                m_v0_inner      = newSource.m_v0_inner;
                m_d0_inner      = newSource.m_d0_inner;

                % Assign new start of outer plume
                psi_0_outer     =  newSource.psi_0_outer;
                Q_0_outer       =  newSource.Q_0_outer;
                m_class_0_outer = newSource.m_class_0_outer; 
                m_v0_outer      = newSource.m_v0_outer;
                m_d0_outer      = newSource.m_d0_outer;

                P_0      = newSource.P_0;

                
                % Check to see if still solid MFR in inner plume, if not leave loop and transition to single plume
                if sum(m_class_0_inner) == 0
                    disp('stop running model as nothing left in inner')
                    break
                end

                counter = counter + 1;

                IC      = [m_d0_inner m_v0_inner psi_0_inner Q_0_inner P_0 ...
                           m_class_0_inner(1:const.no_class_inner) m_d0_outer ...
                           m_v0_outer psi_0_outer Q_0_outer ... 
                           m_class_0_outer(1:const.no_class_outer)];                   %initial conditions for ode solver
                
                % Re-run model
                odeSolution = ode15s(@(z,y) doublePlume_return2one_ODE_ED(z,y,const,atmos),Sspan,IC,options);
                
                % Store model output to storage cell
                storage{counter} = odeSolution; 


            else
                break
            end
        end
        
        % Put restarts of runs together
        [Q_outer,psi_outer,m_d_outer,m_v_outer,m_class_outer,m_s_outer,m_outer,Q_inner,psi_inner,m_d_inner,m_v_inner,m_class_inner,m_s_inner,m_inner,P_double,z_double] = putTogether(storage,counter,const);
        
        
    end
    
%% Determine double plume quanities
    
    % Inner plume quantities
    u_inner     = psi_inner./m_inner;                                      % Velocity of the inner plume region (m/s)
    C_B_inner   = (m_d_inner.*const.C_d + m_v_inner.*const.C_v + ...
                  sum(m_s_inner,2).*const.C_s)./m_inner;                   % Heat capacity of the inner plume
    theta_inner = (1./C_B_inner).*(Q_inner./m_inner);                      % Temperature of the inner plume (k)
    n_d_inner   = m_d_inner./m_inner;                                      % Dry dry gas MFR
    n_v_inner   = m_v_inner./m_inner;                                      % Water vapor MFR
    n_s_inner   = m_s_inner./m_inner;                                      % Solid MFR
    R_g_inner   = (n_v_inner./(n_v_inner+n_d_inner)).*const.R_v + ...
                  (1 - n_v_inner./(n_v_inner+n_d_inner)).*const.R_d;       % Specific gas constant of the inner plume      
    rho_g_inner = P_double./R_g_inner./theta_inner;                        % Density of the gas phases of the inner plume
    rho_b_inner = ((n_v_inner+n_d_inner)./rho_g_inner + (1-n_v_inner-n_d_inner)./const.rho_s).^(-1); % Bulk inner plume denisty
    r_inner     = (m_inner./(rho_b_inner.*u_inner)).^(1/2);                % Radius (m)
    
    % Outer plume quantities
    u_outer     = psi_outer./m_outer;                                      % Velocity of the outer plume region (m/s)
    C_B_outer   = (m_d_outer.*const.C_d + m_v_outer.*const.C_v + ...
                  sum(m_s_outer,2).*const.C_s)./m_outer;                   % Heat capacity of the outer plume
    theta_outer = (1./C_B_outer).*(Q_outer./m_outer);                      % Temperature of the outer plume (k)
    n_d_outer   = m_d_outer./m_outer;                                      % Dry gas MFR
    n_v_outer   = m_v_outer./m_outer;                                      % Water vapor MFR
    n_s_outer   = m_s_outer./m_outer;                                      % Solid mass MFR
    R_g_outer   = (n_v_outer./(n_v_outer+n_d_outer)).*const.R_v + ... 
                  (1 - n_v_outer./(n_v_outer+n_d_outer)).*const.R_d;       % Specific gas constant of the outer plume      
    rho_g_outer = P_double./R_g_outer./theta_outer;                        % Density of the gas phases of the outer plume
    rho_b_outer = ((n_v_outer+n_d_outer)./rho_g_outer + (1-n_v_outer-n_d_outer)./const.rho_s).^(-1); % Bulk outer plume denisty
    r_outer     = (r_inner.^2 + m_outer./(rho_b_outer.*u_outer)).^(1/2);   % Radius (m)
        
%% Model in the single plume region

    % Set esps for the single plume zone
    % 1 = add all left in inner to the single plume
    % 0 = add nothing left in inner to the single plume 
    % 2 = add left in inner to the single plume if supported
    if addMass == 1
        disp('addMass is 1 so all is added')
        m_d0_one                              = m_d_inner(end) + m_d_outer(end); % Initial MFR dry air of single plume
        m_v0_one                              = m_v_inner(end) + m_v_outer(end); % Initial MFR water vapor of single plume
        psi_0_one                             = psi_inner(end) + psi_outer(end); % Initial momentum flow rate of single plume
        Q_0_one                               = Q_inner(end) + Q_outer(end); % Initial heat flow rate of single plume
        m_class_0_one(1:const.no_class_inner) = m_class_outer(end,:) + m_class_inner(end,:); % Initial MFR solid of single plume
	elseif addMass == 0
        disp('addMass is 0 so nothing is added')
        m_d0_one                              = m_d_outer(end);            % Initial MFR dry air of single plume
        m_v0_one                              = m_v_outer(end);            % Initial MFR water vapor of single plume
        psi_0_one                             = psi_outer(end);            % Initial momentum flow rate of single plume
        Q_0_one                               = Q_outer(end);              % Initial heat flow rate of single plume
        m_class_0_one(1:const.no_class_inner) = m_class_outer(end,:);      % Initial MFR solid of single plume
	elseif addMass == 2
        disp('addMass is 2 so added depending on u')

        u_outer_lastStep   = u_outer(end);                                 % Set velocity from last modelled height
        rho_B_lastStep     = rho_b_outer(end);                             % Set density from last modelled height

        % Determine what MFR actually goes into outer (i.e. into the one)
        [mass2add] = addFallout(u_outer_lastStep,diameter_outer,const.g,rho_B_lastStep,const.viscosity,const,m_class_inner(end,:));

        % Work out value of velocity and bulk heat capacity from last step of the inner plume
        u_inner_lastStep    = u_inner(end);
        C_B_inner_lastStep  = C_B_inner(end);                

        % Calculate what momentum and energy flow rate should be added to the outer from the falling inner
        psi_added = u_inner_lastStep * (sum(mass2add) + m_d_inner(end) + m_v_inner(end));
        Q_added   = C_B_inner_lastStep * theta_inner(end) * (sum(mass2add) + m_d_inner(end) + m_v_inner(end));

        % set onePM start conditions
        m_class_0_one = m_class_outer(end,:) + mass2add';                  % Initial MFR solid of single plume               
        Q_0_one       = Q_outer(end)   + Q_added;                          % Initial heat flow rate of single plume
        psi_0_one     = psi_outer(end) + psi_added;                        % Initial momentum flow rate of single plume
        m_d0_one      = m_d_inner(end) + m_d_outer(end);                   % Initial MFR dry air of single plume
        m_v0_one      = m_v_inner(end) + m_v_outer(end);                   % Initial MFR water vapor of single plume

	else
         disp('Error when chosing whether to add mass at the top of the inner plume or not')
	end
    
    % Initial pressure at start of the single plume
    P_0_one = P_double(end);

	options = odeset('RelTol',1e-11,'AbsTol',1e-11,'Events',@(z,y) onePlume_stop(z,y,const,atmos));                    
    Sspan = [z_double(end) 40000];                                         % Solve from where double PM finished to z=40000m
    IC    = [m_d0_one m_v0_one psi_0_one Q_0_one P_0_one m_class_0_one(1:const.no_class_inner)]; % Initial conditions
    
    % Run ODE solver for single plume
    odeSolution = ode45(@(z,y) onePlume_ODE(z,y,const,atmos),Sspan,IC,options);
    
    % Get output from ODE solver
    Y_one  = odeSolution.y';
    z_one  = odeSolution.x';
    
    m_d_one     = Y_one(:,1);                                              % ODE output MFR dry gas phase
    m_v_one     = Y_one(:,2);                                              % ODE output MFR wet vapour phase
    psi_one     = Y_one(:,3);                                              % ODE output momentum flow rate
    Q_one       = Y_one(:,4);                                              % ODE output energy flow rate
    P_one       = Y_one(:,5);                                              % ODE output pressure
    m_class_one = Y_one(:,6:6 + const.no_class_inner - 1);                 % ODE outer output MFR of each grainsize in the solid phase
    m_s_one     = sum(m_class_one,2);                                      % ODE outer output MFR of solid phase (i.e. sum of all grainsize)
    
    m_one       = m_s_one + m_d_one + m_v_one;                             % ODE output total MFR in plume of all phases

%% Define atmos properties
    
    % Vector to store atmospheric temperature
    ambT = zeros(1, length(z_double) + length(z_one));                     
    
    % Assign vectors that where stored in the class
    H1         = atmos.H1;                                                 % Height of troposphere (m)
    H2         = atmos.H2;                                                 % Height of stratosphere (m) 
    theta_a0   = atmos.theta_a0;                                           % Initial temperature of atmosphere (k)
    tempGrad_t = atmos.tempGrad_t;                                         % Temperature gradient in troposphere
    tempGrad_s = atmos.tempGrad_s;                                         % Temperature gradient in troposphere
    
    % Calculate temp for the double plume region
    for i=1:length(z_double)
        %ambient temperature [K]
        if z_double(i) <= H1
            ambT(i) = theta_a0 - tempGrad_t*z_double(i);
        elseif (z_double(i)> H1) && (z_double(i)<=H2)
            ambT(i) = theta_a0 - tempGrad_t*H1;
        else
            ambT(i) = theta_a0 - tempGrad_t*H1 + tempGrad_s*(z_double(i) - H2);
        end
    end
    
    % Calculate temp for the single plume region
    for i=1:length(z_one)
        %ambient temperature [K]
        if z_one(i) <= H1
            ambT(i + length(z_double)) = theta_a0 - tempGrad_t*z_one(i);
        elseif (z_one(i)> H1) && (z_one(i)<=H2)
            ambT(i + length(z_double)) = theta_a0-tempGrad_t*H1;
        else
            ambT(i + length(z_double)) = theta_a0 - tempGrad_t*H1+tempGrad_s*(z_one(i) - H2);
        end
    end
    

    w_a       = 0;
    eps       =  const.R_d/const.R_v;
    rho_ab    =  [P_double;P_one]./(const.R_v.*ambT').*(1 + w_a)./(w_a + eps);      % Environmental density
    rhophi_av = P_double./(const.R_v.*ambT).*w_a./(w_a + eps);             % Volume fraction dry     
    rhophi_ad = P_double./(const.R_d.*ambT).*1./(w_a + eps);               % Volume fraction wet
  
%% Define single plume region quanities

    u_one     = psi_one./m_one;                                            % Velocity of the single plume
    C_B_one   = (m_d_one.*const.C_d + m_v_one.*const.C_v + sum(m_s_one,2).*const.C_s)./m_one; % Heat capacity of the single plume
    theta_one = (1./C_B_one).*(Q_one./m_one);                              % Temperature of the single plume
    
    n_d_one   = m_d_one./m_one;                                            % Dry air mass fraction of the single plume
    n_v_one   = m_v_one./m_one;                                            % Water vapor mass fraction of the single plume
    n_s_one   = m_s_one./m_one;                                            % Solid mass fraction of the single plume
    R_g_one   = (n_v_one./(n_v_one + n_d_one)).*const.R_v + (1 - n_v_one./(n_v_one + n_d_one)).*const.R_d; 
    rho_g_one = P_one./R_g_one./theta_one;                                 % Gas density of the single plume
    rho_b_one = ((n_v_one + n_d_one)./rho_g_one + (1 - n_v_one - n_d_one)./const.rho_s).^(-1);  % Bulk density of the single plume
    r_one     = (m_one./(rho_b_one.*u_one)).^(1/2);                        % Radius of the single plume
    
%% Create output
        
    % Add initialMFRs to output
    output.m0_inner = m_0_inner * pi;
    output.m0_outer = m_0_outer * pi;
    
    % Add plume height to output
    output.z_double = z_double;
    output.z_one = z_one;
    
    % Add velocities to output
    output.u_inner = u_inner;
    output.u_outer = u_outer;
    output.u_one = u_one;
    
    % Add temperature to output
    output.temp_inner = theta_inner;
    output.temp_outer = theta_outer;
    output.temp_one = theta_one;
    
    % Add raduii to output
    output.r_inner = r_inner;
    output.r_outer = r_outer;
    output.r_one = r_one;
    
    % Add densities to output
    output.rho_b_inner = rho_b_inner;
    output.rho_b_outer = rho_b_outer;
    output.rho_b_one = rho_b_one;
    
%     % Plot dynamics
%     figure(1)
%         subplot(2,2,1)
%             hold on
%             plot(u_inner,z_double,'b--')
%             plot(u_outer,z_double,'b:')
%             plot(u_one,z_one,'b:')
%         subplot(2,2,2)
%             hold on
%             plot(theta_inner,z_double,'b--')
%             plot(theta_outer,z_double,'b:')
%             plot(theta_one,z_one,'b:')
%         subplot(2,2,3)
%             hold on
%             plot(r_inner,z_double,'b--')
%             plot(r_outer,z_double,'b:')
%             plot(r_one,z_one,'b:')
%         subplot(2,2,4)
%             hold on
%             plot(rho_b_inner,z_double,'b--')
%             plot(rho_b_outer,z_double,'b:')
%             plot(rho_b_one,z_one,'b:')



    end