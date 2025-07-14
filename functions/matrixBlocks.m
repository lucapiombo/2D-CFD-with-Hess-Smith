function blocks = matrixBlocks(triggerFoil, velocities, U_infinity, alpha)
    % Computes influence coefficients for source and vortex panels
    
    [nPoints, nPanels] = size(velocities.U_s);
    n_hat = [-sin(triggerFoil.beta)', cos(triggerFoil.beta)'];
    t_hat = [cos(triggerFoil.beta)', sin(triggerFoil.beta)'];

    %Initialize sub-matrices
    A_s = zeros(nPoints, nPanels);
    a_v = zeros(nPoints,1);
    c_s = zeros(1,nPanels);
    c_v = 0;
    b_s = zeros(nPoints,1);
    
    for i = 1:nPoints % control points
            for j = 1:nPanels % panels
                
                % Costruct the sub-matrices (influence coefficients)
                A_s(i,j) = sum([velocities.U_s(i,j), velocities.V_s(i,j)].*n_hat(i,:));
                a_v(i) = a_v(i) + sum([velocities.U_v(i,j), velocities.V_v(i,j)].*n_hat(i,:));
    
            end
            
            b_s(i) = -sum(U_infinity*[cos(alpha), sin(alpha)].*n_hat(i,:));
    
    end
    
    b_v = -sum(U_infinity*[cos(alpha), sin(alpha)].*(t_hat(1,:)+t_hat(nPanels,:)));

    for j = 1:nPanels
        c_s(1,j) = sum([velocities.U_s(1,j), velocities.V_s(1,j)].*t_hat(1,:)) + ...
                   sum([velocities.U_s(nPanels,j), velocities.V_s(nPanels,j)].*t_hat(nPanels,:));
        c_v = c_v + sum([velocities.U_v(1,j), velocities.V_v(1,j)].*t_hat(1,:)) + ...
              sum([velocities.U_v(nPanels,j), velocities.V_v(nPanels,j)].*t_hat(nPanels,:));
    end

    blocks.A_s = A_s;
    blocks.a_v = a_v;
    blocks.c_s = c_s;
    blocks.c_v = c_v;
    blocks.b_s = b_s;
    blocks.b_v = b_v;
end
