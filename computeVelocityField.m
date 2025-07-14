function [velSource, velVortex] = computeVelocityField(velocity, gamma, q)
    % Computes velocity at field points due to source and vortex panels.

    [nPoints, nPanels] = size(velocity.U_s);
    
    % TOTAL velocity due to SOURCES and VORTEX
    velSource = zeros(nPoints,2);
    velVortex = zeros(nPoints,2);
    
    for i =1:nPoints % points
        for j = 1:nPanels % panels
            
            velSource(i,1) = velSource(i,1) + velocity.U_s(i,j) * q(j);
            velSource(i,2) = velSource(i,2) + velocity.V_s(i,j) * q(j);
            velVortex(i,1) = velVortex(i,1) + velocity.U_v(i,j) * gamma;
            velVortex(i,2) = velVortex(i,2) + velocity.V_v(i,j) * gamma;

        end
    end
end
