function[newSource] = add2Surrounding(odeSolution,const)

    diameter_outer = const.diameter_outer;

%% Assign last values from last ODE solver to the parameters
    
    Y = odeSolution.y';
    
    % Inner
    m_d0_inner = Y(end,1); 
    m_v0_inner = Y(end,2); 
    psi_0_inner = Y(end,3); 
    Q_0_inner = Y(end,4); 
    mass_solid = Y(end,6:6 + const.no_class_inner - 1); 
    
    % Outer
    m_d0_outer = Y(end,6 + const.no_class_inner); 
    m_v0_outer = Y(end,7 + const.no_class_inner); 
    psi_0_outer = Y(end,8 + const.no_class_inner); 
    Q_0_outer = Y(end,9 + const.no_class_inner);  
    m_class_0_outer = Y(end,10 + const.no_class_inner:(10 + const.no_class_inner) + (const.no_class_outer - 1));
    
    P_0 = Y(end,5);      
 
%% Calculate and adjust the psi and Q lost from the inner      
    
    %Calculate
    u_innerLost = psi_0_inner / (sum(mass_solid) + m_d0_inner + m_v0_inner); %work out the last value of u of the inner
    C_B_innerLost   = (m_d0_inner.*const.C_d + m_v0_inner.*const.C_v + ...   %work out the last value of c_B of the inner
                   sum(mass_solid).*const.C_s)./(sum(mass_solid) + m_d0_inner + m_v0_inner);                   
    temp_innerLost = (1./C_B_innerLost).*(Q_0_inner./(sum(mass_solid) + m_d0_inner + m_v0_inner)); %work out the last value of temp of the inner                      

    psi_lost = sum((u_innerLost*mass_solid(odeSolution.ie-2)));     %work out the amount of psi lost for given sie lost for the next step (i.e. start of the ne model)
    Q_lost = sum(const.C_s*temp_innerLost*mass_solid(odeSolution.ie-2)); %work out the amount of energy lost for given sie lost for the next step (i.e. start of the ne model)

    % Adjust the psi,mass and energy of the starting position of the next model run

% % % % % %     Q_0_inner = Q_0_inner - Q_lost;
% % % % % %     psi_0_inner = psi_0_inner - psi_lost;
% % % % % %                 
                
%% Determine if what has fallout out of the inner goes into the outer
    
    % Mass of what has fallen out of inner
    massFallout_inner = mass_solid(odeSolution.ie-2);
                    
    %calclaute variables of the surrounding plume
    u_outer = psi_0_outer / (sum(m_class_0_outer) + m_d0_outer + m_v0_outer); %work out the last value of u of the inner
    C_B_outer   = (m_d0_outer.*const.C_d + m_v0_outer.*const.C_v + ...
                   (sum(m_class_0_outer)).*const.C_s)./(m_d0_outer + m_v0_outer + (sum(m_class_0_outer)));                   %Heat capacity of the outer plume
    theta_outer = (1./C_B_outer).*(Q_0_outer./(m_d0_outer + m_v0_outer + (sum(m_class_0_outer))));
    g = const.g;
    n_d_outer   = m_d0_outer./(m_d0_outer + m_v0_outer + (sum(m_class_0_outer)));                                      %Dry gas mass flux
    n_v_outer   = m_v0_outer./(m_d0_outer + m_v0_outer + (sum(m_class_0_outer)));                                      %Water vapor mass flux
    %                     n_s_outer   = (sum(m_class_0_outer)./(m_d0_outer + m_v0_outer + (sum(m_class_0_outer)));  
    R_g_outer   = (n_v_outer./(n_v_outer+n_d_outer)).*const.R_v + ... 
    (1 - n_v_outer./(n_v_outer+n_d_outer)).*const.R_d;       %Specific gas constant of the outer plume      
    rho_g_outer = P_0./R_g_outer./theta_outer;                        %Density of the gas phases of the outer plume
    rho_B_outer = ((n_v_outer+n_d_outer)./rho_g_outer + (1-n_v_outer-n_d_outer)./const.rho_s).^(-1); %Bulk outer plume denisty
    viscosity = const.viscosity;
    
    % Detemine if any mass is added to the surrounding
    [mass2add_2outer] = addFallout_Single(u_outer,diameter_outer(odeSolution.ie-2),g,rho_B_outer,viscosity,const,massFallout_inner);
    
    % If mass from inner is not supported by outer, adjust accordingly
    if mass2add_2outer == 0
        m_d0_outer = Y(end,6 + const.no_class_inner); 
        m_v0_outer = Y(end,7 + const.no_class_inner); 
        psi_0_outer = Y(end,8 + const.no_class_inner); 
        Q_0_outer = Y(end,9 + const.no_class_inner);  
        m_class_0_outer = Y(end,10 + const.no_class_inner:(10 + const.no_class_inner) + (const.no_class_outer - 1));
    else
        m_d0_outer = Y(end,6 + const.no_class_inner); 
        m_v0_outer = Y(end,7 + const.no_class_inner);  
        psi_0_outer = Y(end,8 + const.no_class_inner) + psi_lost;
        Q_0_outer = Y(end,9 + const.no_class_inner) + Q_lost;
        m_class_0_outer = Y(end,10 + const.no_class_inner:(10 + const.no_class_inner) + (const.no_class_outer - 1));
        m_class_0_outer(odeSolution.ie-2) = m_class_0_outer(odeSolution.ie-2) + massFallout_inner;
    end

%% Adjust inner

    % Adjust inner grain-size distribution
    mass_solid(odeSolution.ie-2) = 0;
    m_class_0_inner              =  mass_solid;   
    
    % Adjust inner energy
    Q_0_inner = Q_0_inner - Q_lost;
    
    % Adjust inner psi
    psi_0_inner = psi_0_inner - psi_lost;
    
%% Assign to output

    % Assign new start of inner plume
    newSource.psi_0_inner     =  psi_0_inner;
    newSource.Q_0_inner       =  Q_0_inner;
    newSource.m_class_0_inner = m_class_0_inner; 
    newSource.m_v0_inner      = m_v0_inner;
    newSource.m_d0_inner      = m_d0_inner;
    
    % Assign new start of outer plume
    newSource.psi_0_outer     =  psi_0_outer;
    newSource.Q_0_outer       =  Q_0_outer;
    newSource.m_class_0_outer = m_class_0_outer; 
    newSource.m_v0_outer      = m_v0_outer;
    newSource.m_d0_outer      = m_d0_outer;
    
    newSource.P_0 = P_0;
    
end