function [u_s,v_s,u_v,v_v,U] = flow_velocity(x_field,z_field,x,z,beta,q,gamma,U_infinity, alpha,l_panel)
%--------------------------------------------------------------------------
%This function is used to calculate the real velociy perceived in a certain
%coordinate of your field (x_field,z_field) once you have got the correct
%value of the strength of the sources and vortices on the panels where the
%extrema have coordinates (x,z), length l_panel and orientation w.r.t. the
%gloabl system beta.

%INPUTS: orientation of each panel with respect to a global reference 
%system (beta), coordinate of the desired points of the field (x_field,z_field), 
%extrema (x,z), length of the panels (l_panel), stream flow velocity and 
%inclination (U_infinity and alpha).

%OUTPUTS: Velocity perceived in point (x_field,z_field) due to:
%1) The sources alone u_s, v_s: these are matrix of as rows as the desired
%number of points to evaluate and as number of columns as the number of
%panels considered. Each element define the value of the x or z component of
%velocity in the point defined by the number of the row due to the source 
%distribution defined by the number of the column.
%2) The vortices alone [u_v, v_v]: these are matrix of as rows as the desired
%number of points to evaluate and as number of columns as the number of
%panels considered. Each element define the value of the x or z component of
%velocity in the point defined by the number of the row due to the vortex 
%distribution defined by the number of the column.
%3) Sources, Vortices, Streamflow toghether U: thi is a matrix of two
%components where the columns identify the x or z component of the velocity
%related to the point of the i-th row.
%--------------------------------------------------------------------------



%Initialize elements:
n = length(x_field);

for i = 1:n %field points
    for j = 1:(length(x)-1) %panels
        
        %1) -------- Determine the LOCAL coordinate of the i-th field point w.r.t all the panels: --------
        Q = [cos(beta(j)), sin(beta(j)); -sin(beta(j)), cos(beta(j))];
        local_c = Q * [x_field(i)-x(j); z_field(i)-z(j)];
        eta_c = local_c(1);
        ni_c = local_c(2);


        %Radii of j-th panel w.r.t i-th field point 
        r1 = sqrt(eta_c^2 + ni_c^2);
        r2 = sqrt((l_panel(j)-eta_c)^2 + ni_c^2);
        
        %Angles of the radius w.r.t. LOCAL reference system
        teta1 = atan2(ni_c, eta_c);
        teta2 = atan2(ni_c, eta_c-l_panel(j));
        %conditions to atan2 to do not screw-up
        if abs(teta1)<10^(-12) && abs(teta2)>3
            teta1 = 0;
            teta2 = pi;
        end
        if abs(teta2)<10^(-12) && abs(teta1)>3
            teta2 = 0;
            teta1 = -pi;
        end
        
        
        %2) -------- Determine the induced velocity in LOCAL reference (for the SOURCE): --------
        u_s_l(i,j) = -0.5/pi*log(r2/r1);
        v_s_l(i,j) = 0.5/pi*(teta2-teta1);
        
        
        %3) -------- Go back to GLOBAL: --------
        Q = [cos(beta(j)), sin(beta(j)); -sin(beta(j)), cos(beta(j))];
        vel_global = Q' * [u_s_l(i,j); v_s_l(i,j)];
        
%--------------- Velocity (due to a SOURCE): ---------------
        u_s(i,j) = vel_global(1);
        v_s(i,j) = vel_global(2);
        
        
%--------------- Velocity due to a VORTEX: ---------------
        u_v(i,j) = v_s(i,j);
        v_v(i,j) = -u_s(i,j);
        
    end  
end




U_source = zeros(n,2);
U_vortex = zeros(n,2);

for i =1:n %field points
    for j = 1:1:(length(x)-1) %panels
        
        U_source(i,1) = U_source(i,1)+u_s(i,j)*q(j);
        U_source(i,2) = U_source(i,2)+v_s(i,j)*q(j);
        U_vortex(i,1) = U_vortex(i,1)+u_v(i,j)*gamma;
        U_vortex(i,2) = U_vortex(i,2)+v_v(i,j)*gamma;
        

    end
%--------------- TOTAL velocity: ---------------
    U(i,:) = U_infinity*[cos(alpha), sin(alpha)] + U_source(i,:) + U_vortex(i,:);
    
end

end
