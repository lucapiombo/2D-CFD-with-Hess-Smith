function [A_s,a_v,c_s,c_v,b_s,b_v] = Influence_coefficient(beta,x_c,z_c,x,z,l_panel,U_infinity,alpha)
%--------------------------------------------------------------------------
%This function is used to calculate the influence coefficients associated
%to the centroides you pass as input (x_c,z_c), due to the vortex and
%sources distribution of the geometry you pass as imput (beta, (x,z), l_panel).

%INPUTS: orientation of each panel with respect to a global reference 
%system (beta), coordinate of centroids (x_c,z_c), extrema (x,z), length of
%the panels (l_panel), stream flow velocity and inclination (U_infinity
%and alpha).

%OUTPUTS: Influence coefficients A_s, a_v, c_s, c_v, b_s, b_v
%--------------------------------------------------------------------------



%Initialize elements:
n = length(x_c);
a_v = zeros(n,1);
c_v = 0;
n_hat = [-sin(beta)', cos(beta)'];
t_hat = [cos(beta)', sin(beta)'];

for i = 1:n %control points
        for j = 1:n %panels
        
        %1) -------- Determine the LOCAL coordinate of the i-th control point w.r.t all the panels: --------
        Q = [cos(beta(j)), sin(beta(j)); -sin(beta(j)), cos(beta(j))];
        local_c = Q * [x_c(i)-x(j); z_c(i)-z(j)];
        eta_c = local_c(1);
        ni_c = local_c(2);

        
        %Radii of j-th panel w.r.t i-th centroid 
        r1 = sqrt(eta_c^2 + ni_c^2);
        r2 = sqrt((l_panel(j)-eta_c)^2 + ni_c^2);
        
        %Angles of the radius w.r.t. LOCAL reference system
        teta1 = atan2(ni_c ,eta_c);
        teta2 = atan2(ni_c, eta_c-l_panel(j));     
        %conditions to atan2 to do not screw-up
        if abs(teta1)<10^(-12) && abs(teta2)>3 
            teta1 = 0;
            teta2 = pi;
        end
        if abs(teta2)<10^(-12) && abs(teta1)>3
            teta2 = 0;
            teta1  = -pi;
        end

   
        %2) -------- Determine the induced velocity in LOCAL reference (for a UNITARY SOURCE): --------
        u_s_l(i,j) = -0.5/pi*log(r2/r1);
        v_s_l(i,j) = 0.5/pi*(teta2-teta1);
        
        
        %3) -------- Go back to GLOBAL: --------
        Q = [cos(beta(j)), sin(beta(j)); -sin(beta(j)), cos(beta(j))];
        vel_global = Q' * [u_s_l(i,j); v_s_l(i,j)];
        u_s(i,j) = vel_global(1);
        v_s(i,j) = vel_global(2);
        
        %Velocity (due to a UNITARY VORTEX) directly in global evaluated simply like this:
        u_v(i,j) = v_s(i,j);
        v_v(i,j) = -u_s(i,j);
        
        
        
        
        
%--------------- Costruct the sub-matrices (influence coefficients): ---------------
        A_s(i,j) = sum([u_s(i,j),v_s(i,j)].*n_hat(i,:));
        a_v(i) = a_v(i) + sum([u_v(i,j), v_v(i,j)].*n_hat(i,:));

        
        end
        
        b_s(i) = -sum(U_infinity*[cos(alpha), sin(alpha)].*n_hat(i,:));

end

b_v = -sum(U_infinity*[cos(alpha), sin(alpha)].*(t_hat(1,:)+t_hat(n,:)));

for j = 1:n
    c_s(j) = sum([u_s(1,j), v_s(1,j)].*t_hat(1,:)) + sum([u_s(n,j), v_s(n,j)].*t_hat(n,:));
    c_v = c_v + sum([u_v(1,j), v_v(1,j)].*t_hat(1,:)) + sum([u_v(n,j), v_v(n,j)].*t_hat(n,:));
end

end