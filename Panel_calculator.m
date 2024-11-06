function [x,z,x_c,z_c,l_panel,beta] = Panel_calculator(c,t,N_panel,m,p,AoA)
%--------------------------------------------------------------------------
%This function is used create the geometry of a 4 digit NACA airfoil by
%means of a certain number of panels using half-cosine distribution

%INPUTS: chord (c), thickness (t), half of the number of panels desired 
%(N_panel), maximum chamber (m), location of the maximum chamber (m), angle 
%of attack (AoA)

%OUTPUTS: coordinates of the extrema of the panels (x,z), the centroids
%(x_c,z_c), length of each panel (l_panel) and orientation of each panel
%with respect to a global reference system (beta)
%--------------------------------------------------------------------------



%--------------- Find the (x,y) coordinates of extrema: ---------------

%Create and anle-spacing:
delta_theta = pi/(N_panel); % Delta difference between each angles
theta = 0:delta_theta:2*pi; % Angles vector


%Corresponding x-coordinates (half-cosine method):
x = c/2*(1 - cos(theta));% Define the x-coordinates on the chordline that will
                         % correspond to the x-coordinates of your extreme
                         % points of each panel

x_t = x(1:N_panel+1); % For the top side and goes from the leading 
                      % edge to the trailing edge

%You have to distinguish between symmmetric and not:
if m == 0 && p == 0  %for symmetric airfoils

    x_b = x(N_panel+1:2*N_panel); %For the bottom and goes from the 
                                %trailing edge to the leading edge
    
    %Calculate the z-coordinates:
    z_t = 5*t*(0.2969*sqrt(x_t/c) + (-0.1260)*(x_t/c) + ...
        (-0.3516)*(x_t/c).^2 + 0.2843*(x_t/c).^3 + ...
        (-0.1036)*(x_t/c).^4); %thickness distribution for the airfoil (for
                               %the top part)
    z_b = -(z_t(end:-1:2));%thickness distribution for the airfoil (for
                           %the top part)
    
    %Calculate extreme coordinates of the panel:
    x = horzcat(x_b,x_t); %makes the concatenation of the two vectors in input
    z = horzcat(z_b,z_t); %makes the concatenation of the two vectors in input
    z(1) = 0;
    z(end) = 0;
    
else %for not symmetric airfoils
    
    %Divide into two parts: that where x is less than p*c and that where is
    %equal to or greater than p*c.
    %In that way I can evaluate the equation of the airfoil as the union of
    %these two parts.
    
    
    %Calculate the z-coordinates:
    z_t = 5*t*(0.2969*sqrt(x_t/c) + (-0.1260)*(x_t/c) + ...
        (-0.3516)*(x_t/c).^2 + 0.2843*(x_t/c).^3 + ...
        (-0.1036)*(x_t/c).^4); % Thickness distribution
       
    
    %Determine the transition, that is the point where p*c = x
    [~,x_transition] = min(abs(x_t - p*c)); %x where p*c = c

    %For x < p*c
    z_camb1 = m*x_t(1:x_transition-1)/p^2.*(2*p - ...
        x_t(1:x_transition-1)/c); %Chamber thickness distribution
    grad = 2*m/(p^2)*(p - x_t(1:x_transition-1)/c); %angle of inclination 1
    kappa1 = atan(grad); %slope coefficient 1
    
    %For x >= p*c
    z_camb2 = m*(c - x_t(x_transition:end))/((1 - p)^2).*(1 + ...
        x_t(x_transition:end)/c - 2*p); %Chamber thickness distribution
    grad = 2*m/((1 - p)^2)*(p - x_t(x_transition:end)/c); %angle of inclination 2
    kappa2 = atan(grad); %slope coefficient 2
    
    z_camb = horzcat(z_camb1,z_camb2); %thickness of the chamber vector
    kappa = horzcat(kappa1,kappa2); %slope coefficients vector
    
    
    x_u = x_t - z_t.*sin(kappa); %x-coordinates upper surface
    x_l = fliplr(x_t + z_t.*sin(kappa)); %x-coordinates lower surface
    
    z_u = z_camb + z_t.*cos(kappa); %z-coordinates upper surface
    z_l = fliplr(z_camb - z_t.*cos(kappa)); %z-coordinates lower surface
    
    %Calculate extreme coordinates of the panel:
    z = horzcat(z_l(1:end-1),z_u);
    x = horzcat(x_l(1:end-1),x_u);      
end

%--------------- Impose the angle of attack: ---------------
n = N_panel*2;
for i = 1:n+1
        x(i) = sum([cos(AoA), sin(AoA)].*[x(i), z(i)]);
        z(i) = sum([-sin(AoA), cos(AoA)].*[x(i), z(i)]);
end

%--------------- Define control points: ---------------
for i =1:n
    
    dx(i) = x(i+1)-x(i);
    dz(i) = z(i+1)-z(i); 
    
    %Control points:
    x_c(i) = dx(i)/2 + x(i); %x-coordinate
    z_c(i) = dz(i)/2 + z(i); %y-coordinate
    
    
    %--------------- Length of each panel ---------------
    l_panel(i) = sqrt(dx(i).^2 + dz(i).^2);
    
    
    %--------------- Define the angle beta ---------------
    beta(i) = atan2(dz(i),dx(i)); %the orientation of the panels with respect to the global
end

end