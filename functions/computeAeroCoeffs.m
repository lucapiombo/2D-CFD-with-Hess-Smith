function [Cp, Cl, Cm_LE] = computeAeroCoeffs(airfoil, U_source,U_vortex, U_infinity, alpha, c)
    % Compute Cp, Cl, and Cm_LE for given panels using total velocity at control points.
    
    n_hat = [-sin(airfoil.beta)', cos(airfoil.beta)'];
    t_hat = [cos(airfoil.beta)', sin(airfoil.beta)'];
    n = length(airfoil.x_c);
    
    %Initialize variables:
    U_t = zeros(1,n); %Tangent projection of flow velocity
    
    for i =1:n
        U_t(i) = sum((U_infinity*([cos(alpha), sin(alpha)])+U_source(i,:)+U_vortex(i,:)).*t_hat(i,:)); %velocità U effettiva su ogni pannello
    end  
    
    %PRESSURE COEFFICIENT of each panel:
    Cp = 1 - (U_t/U_infinity).^2; 
    
    %LIFT COEFFICIENT:
    Cl=0;
    for i=1:n
        Cl = Cl + Cp(i)*airfoil.l_panel(i)/c*n_hat(i,2);
    end
     Cl = -Cl;
    
    %MOMENT COEFFICIENT:
    Cm_LE = 0;
    for i = 1:n
        Cm_LE = Cm_LE + Cp(i)*airfoil.l_panel(i)/c^2*(airfoil.x_c(i)*n_hat(i,2) ...
            -airfoil.y_c(i)*n_hat(i,1));
    end
    
    Cm_LE = -Cm_LE;

end
