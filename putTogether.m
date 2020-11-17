function[Q_outer,psi_outer,m_d_outer,m_v_outer,m_class_outer,m_s_outer,m_outer,Q_inner,psi_inner,m_d_inner,m_v_inner,m_class_inner,m_s_inner,m_inner,P_double,z_double] = putTogether(storage,counter,const)
%put all the seperate ode solver results together

    % Determine how many time looped
    noLoops = counter;
    
    % Create model output
    z_double = zeros(5000,1);
    P_double = zeros(5000,1);
    
    m_v_inner = zeros(5000,1);
    m_d_inner = zeros(5000,1);
    m_class_inner = zeros(5000,const.no_class_inner);
    Q_inner = zeros(5000,1); 
    psi_inner = zeros(5000,1); 
    
    m_v_outer = zeros(5000,1);
    m_d_outer = zeros(5000,1);
    m_class_outer = zeros(5000,const.no_class_inner);
    Q_outer = zeros(5000,1); 
    psi_outer = zeros(5000,1); 
    
    tally = 0;
    counter = 1;
    
    % Loop througth cell anf extract
    for i = 1:noLoops
        tempODEsolution = storage{i};
        
        %% Get outputs from ode Solver
        Y        = tempODEsolution.y';                                             %output matrix from ode solver
        z_double_1 = tempODEsolution.x';                                             %z array used in the ode solver for the above matrix
        
        array_length = length(z_double_1);
        tally = tally + array_length;
        
        % Inner ODEs output i.e. respective paremter values in the plume with z
        m_d_inner_1     = Y(:,1);                                                %ode inner output mer dry gas phase
        m_v_inner_1     = Y(:,2);                                                %ode inner output mer wet gas phase
        psi_inner_1     = Y(:,3);                                                %ode inner output momentum
        Q_inner_1       = Y(:,4);                                                %ode inner output energy
        P_double_1      = Y(:,5);                                                %ode output pressure
        m_class_inner_1 = Y(:,6:6 + const.no_class_inner - 1);                   %ode inner output mer of each grainsize in the solid phase

        % Outer ODEs output
        m_d_outer_1     = Y(:,6 + const.no_class_inner);                         %ode outer output mer dry gas phase
        m_v_outer_1     = Y(:,7 + const.no_class_inner);                         %ode outer output mer wet vapour phase
        psi_outer_1     = Y(:,8 + const.no_class_inner);                         %ode outer output momentum
        Q_outer_1       = Y(:,9 + const.no_class_inner);                         %ode outer output energy
        m_class_outer_1 = Y(:,10 + const.no_class_inner:(10 + const.no_class_inner) + (const.no_class_outer - 1)); %ode outer output mer of each grainsize in the solid phase

%         figure(23)
%         hold on 
%         plot(P_double_1,z_double_1,'-o')
    
        %% Add to out put array
        
        z_double(counter:counter+array_length-1) = z_double_1;
        
        P_double(counter:counter+array_length-1) = P_double_1;
        
        m_v_inner(counter:counter+array_length-1) = m_v_inner_1;
        m_d_inner(counter:counter+array_length-1) = m_d_inner_1;
        m_class_inner(counter:counter+array_length-1,:) = m_class_inner_1;
        psi_inner(counter:counter+array_length-1) = psi_inner_1;
        Q_inner(counter:counter+array_length-1) = Q_inner_1;
        
        m_v_outer(counter:counter+array_length-1) = m_v_outer_1;
        m_d_outer(counter:counter+array_length-1) = m_d_outer_1;
        m_class_outer(counter:counter+array_length-1,:) = m_class_outer_1;
        psi_outer(counter:counter+array_length-1) = psi_outer_1;
        Q_outer(counter:counter+array_length-1) = Q_outer_1;
        
        counter = counter + array_length;
    end
    
    P_double(z_double == 0) = [];
    
    
    
    m_v_inner(z_double == 0) = [];
    m_d_inner(z_double == 0) = [];
    m_class_inner(z_double == 0,:) = [];
    
    Q_inner(z_double == 0) = [];
    psi_inner(z_double == 0) = [];
    m_s_inner = sum(m_class_inner,2);
    m_inner = m_s_inner + m_v_inner + m_d_inner;
    
    m_v_outer(z_double == 0) = [];
    m_d_outer(z_double == 0) = [];
    m_class_outer(z_double == 0,:) = [];
    Q_outer(z_double == 0) = [];
    psi_outer(z_double == 0) = [];
    m_s_outer = sum(m_class_outer,2);
    m_outer = m_s_outer + m_v_outer + m_d_outer;
    
    z_double(z_double == 0) = [];
    

%     figure(23)
%         hold on 
%         plot(Q_inner,z_double,'-*')
%         pause(0.1)

end