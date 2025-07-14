function velocities = inducedSpeeds(targetFoil, triggerFoil)
    % Computes influence coefficients for source and vortex panels

    nPoints = length(targetFoil.x_c);
    nPanels = length(triggerFoil.x) - 1;

    % Initialize velocity sources and vortices
    U_s = zeros(nPoints, nPanels);
    V_s = zeros(nPoints, nPanels);
    U_v = zeros(nPoints, nPanels);
    V_v = zeros(nPoints, nPanels);
    
    for i = 1:nPoints % control points
            for j = 1:nPanels % panels
            
                %(1) Determine the LOCAL coordinate of the i-th control point w.r.t all the panels
                Q = [cos(triggerFoil.beta(j)), sin(triggerFoil.beta(j)); -sin(triggerFoil.beta(j)), cos(triggerFoil.beta(j))];
                local_c = Q * [targetFoil.x_c(i)-triggerFoil.x(j); targetFoil.y_c(i)-triggerFoil.y(j)];
                eta_c = local_c(1);
                ni_c = local_c(2);
        
                % Radii of j-th panel w.r.t i-th centroid 
                r1 = sqrt(eta_c^2 + ni_c^2);
                r2 = sqrt((triggerFoil.l_panel(j)-eta_c)^2 + ni_c^2);
                
                % Angles of the radius w.r.t. LOCAL reference system
                teta1 = atan2(ni_c ,eta_c);
                teta2 = atan2(ni_c, eta_c-triggerFoil.l_panel(j));     
                
                if abs(teta1)<10^(-12) && abs(teta2)>3 % conditions to atan2 to do not screw-up
                    teta1 = 0;
                    teta2 = pi;
                end
                if abs(teta2)<10^(-12) && abs(teta1)>3 % conditions to atan2 to do not screw-up
                    teta2 = 0;
                    teta1  = -pi;
                end
    
                %(2) Determine the induced velocity in LOCAL reference (for a UNITARY SOURCE)
                u_s_l = -0.5/pi*log(r2/r1);
                v_s_l = 0.5/pi*(teta2-teta1);
                
                %(3) Go back to GLOBAL
                Q = [cos(triggerFoil.beta(j)), sin(triggerFoil.beta(j)); -sin(triggerFoil.beta(j)), cos(triggerFoil.beta(j))];
                vel_global = Q' * [u_s_l; v_s_l];

                % Velocity due to each SOURCE
                U_s(i,j) = vel_global(1);
                V_s(i,j) = vel_global(2);
                
                 % Velocity due to each VORTEX
                U_v(i,j) = V_s(i,j);
                V_v(i,j) = -U_s(i,j);
                
            end
    end

    velocities.U_s = U_s;
    velocities.V_s = V_s;
    velocities.U_v = U_v;
    velocities.V_v = V_v;
end
