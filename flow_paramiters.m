function [Cp,Cl,Cm_LE] = flow_paramiters(beta,n,l_panel,x_field,z_field,c,U_source,U_vortex,U_infinity,alpha)
%--------------------------------------------------------------------------
%This function is used to calculate the value of the flow field paramiters
%observed in the points of coordinates (x_field,z_field)

%INPUTS: orientation of each panel with respect to a global reference 
%system (beta), number of panels (n), length of the panels (l_panel)
%coordinate of the desired points of the field (x_field,z_field), chord (c)
%strength of the sources and vorices (q and gamma), u_s, v_s, u_v, v_v, 
%U_infinity, alpha.

%OUTPUTS: Pressure coefficient Cp of the points (x_field,z_field) and if they
%are the centroides of your body a then you also get the lift coefficient
%Cl and the moment coefficient with respect to the leading edge Cm_LE
%--------------------------------------------------------------------------



%Define the normal and tangent to each panel:
n_hat = [-sin(beta)', cos(beta)'];
t_hat = [cos(beta)', sin(beta)'];


%1)-------- Compute the TANGENT component of the VELOCITY U_t: --------
%Initialize variables:
U_t = zeros(1,n); %Tangent projection of flow velocity

for i =1:n
    U_t(i) = sum((U_infinity*([cos(alpha), sin(alpha)])+U_source(i,:)+U_vortex(i,:)).*t_hat(i,:)); %velocità U effettiva su ogni pannello
end 



%2)------------ COEFFICIENTS: ------------


%PRESSURE COEFFICIENT of each panel:
Cp = 1 - (U_t/U_infinity).^2; 



%LIFT COEFFICIENT:
Cl=0;
for i=1:n
    Cl = Cl + Cp(i)*l_panel(i)/c*n_hat(i,2);
end
 Cl = -Cl;

%CL_jouk = 2*GAMMA./(U_infinity*([cos(alpha), sin(alpha)]));


%MOMENT COEFFICIENT:
Cm_LE = 0;
for i = 1:n
    Cm_LE = Cm_LE + Cp(i)*l_panel(i)/c^2*(x_field(i)*n_hat(i,2)-z_field(i)*n_hat(i,1)); %non sono sicuro
end

Cm_LE = -Cm_LE;

end