%% Set-up inputs
clc,clear,close
addpath('functions');

% --------- Flow paramiters Inputs: ---------
U_infinity = 1; %Velocity of the streamflow
alpha = 0;      %Angle of the stream flow

% --------- Airfoil Inputs: ---------
%AIRFOIL 1:
nPanels(1) = 100;        % Number of panels
NACA4(1) = 0012;        % NACA airfoil denomination ( 4 digit: XXXX )
c(1) = 1;               %chord of main element
AoA(1) = deg2rad(0);  %Angle of attack (deg)

%AIRFOIL 2:
nPanels(2) = 100;        % Number of panels
NACA4(2) = 3312;        % NACA airfoil denomination ( 4 digit: XXXX )
c(2) = 0.4;             %chord of main element
AoA(2) = deg2rad(2.5);    %Angle of attack (deg)

% --------- Airfoil position: ---------
h = 0.2;    %heigth of the first airfoil
x12 = 1.1;  %distance between two airfoils
y12 = 0.1;  %height between two airfoils

%% Build airfoil

% Generate the airfoil panels
figure()
for i = 1:2
    
    airfoils(i) = NACA4_generator(NACA4(i), c(i), AoA(i), nPanels(i));

    % Adjust airfoil 1 position
    if i == 1
        airfoils(i).y = airfoils(i).y + h;
        airfoils(i).y_c = airfoils(i).y_c + h;
    
    % Adjust airfoil 2 position
    elseif i == 2
        airfoils(i).y = airfoils(i).y + (y12 + h);
        airfoils(i).x = airfoils(i).x + x12;
        airfoils(i).y_c = airfoils(i).y_c + (y12 + h);
        airfoils(i).x_c = airfoils(i).x_c + x12; 
    end

    % show airfoil geometry:
    hold on
    plot(airfoils(i).x, airfoils(i).y)
    plot(airfoils(i).x_c, airfoils(i).y_c,'x')
end

% Add images foils:
for i = 1:2
    symm_foil(i).x = airfoils(i).x;
    symm_foil(i).y = -airfoils(i).y;
    symm_foil(i).x_c = airfoils(i).x_c;
    symm_foil(i).y_c = -airfoils(i).y_c;
    symm_foil(i).beta = airfoils(i).beta;
    symm_foil(i).l_panel = airfoils(i).l_panel;

    plot(symm_foil(i).x, symm_foil(i).y)
    plot(symm_foil(i).x_c, symm_foil(i).y_c,'x')
end
yline(0, 'k-', 'LineWidth', 2);
grid on
axis equal
legend('panels','control points')
hold off


%% Evaluation flow paramiters

% Compute coefficients due to panels induced velocities (source and vortex):
for i = 1:2 % targetFoil
    for j= 1:2 % triggerFoil
        velocities(i,j) = inducedSpeeds(airfoils(i), airfoils(j));
        real_coeffs(i,j) = matrixBlocks(airfoils(j), velocities(i,j), U_infinity, alpha);

        symm_velocities(i,j) = inducedSpeeds(airfoils(i), symm_foil(j));
        symm_coeffs(i,j) = matrixBlocks(symm_foil(j), symm_velocities(i,j), U_infinity, alpha);

        coeffs(i,j).A_s = real_coeffs(i,j).A_s + symm_coeffs(i,j).A_s;
        coeffs(i,j).a_v = real_coeffs(i,j).a_v - symm_coeffs(i,j).a_v;
        coeffs(i,j).c_s = real_coeffs(i,j).c_s + symm_coeffs(i,j).c_s;
        coeffs(i,j).c_v = real_coeffs(i,j).c_v - symm_coeffs(i,j).c_v;
    end
end

% Build and solve the linear sistem 
A = zeros(nPanels(1)+nPanels(2)+2,nPanels(1)+nPanels(2)+2);

%Build up the system:
A(1:nPanels(1), 1:nPanels(1)) = coeffs(1,1).A_s; 
A(1:nPanels(1),(nPanels(1)+1):(end-2)) = coeffs(1,2).A_s;
A((nPanels(1)+1):(end-2),1:nPanels(1)) = coeffs(2,1).A_s;
A((nPanels(1)+1):(end-2),(nPanels(1)+1):(end-2)) = coeffs(2,2).A_s; 
A(1:nPanels(1),end-1:end) = [coeffs(1,1).a_v, coeffs(1,2).a_v];
A((nPanels(1)+1):end-2,end-1:end) = [coeffs(2,1).a_v, coeffs(2,2).a_v];
A(end-1:end,1:end-2) = [coeffs(1,1).c_s, coeffs(1,2).c_s; coeffs(2,1).c_s, coeffs(2,2).c_s];
A(end-1:end,end-1:end) = [coeffs(1,1).c_v, coeffs(1,2).c_v; coeffs(2,1).c_v, coeffs(2,2).c_v];
b = [real_coeffs(1,1).b_s; real_coeffs(2,2).b_s; real_coeffs(1,1).b_v; real_coeffs(2,2).b_v];

%Define the solution:
solution = A\b;
q(:,1) = solution(1:nPanels(1));            %Strength front Airfoil
q(:,2) = solution((nPanels(1)+1):end-2);    %Strength rear Airfoil
gamma(1) = solution(end-1);                 %Front Airfoil
gamma(2) = solution(end);                   %Rear Airfoil


%% Visualize results on velocity

% ----- Plot settings: --------
xLim=[-0.5, 2.5];
yLim=[0,2];
mSize = 100;
nLines = 40;

% discretization
xm = linspace(xLim(1), xLim(2), mSize);
ym = linspace(yLim(1), yLim(2), mSize);
[X,Y] = meshgrid(xm, ym);

for i = 1:2
    airfoil_mask = ~inpolygon(X, Y, airfoils(i).x, airfoils(i).y);
    airfoil_mask(airfoil_mask==0) = 0;
    vel_U = zeros(length(xm),length(ym));
    vel_V = zeros(length(xm),length(ym));
end

for i = 1:length(xm)
    for j = 1:length(ym)

        target.x_c = X(i,j);
        target.y_c = Y(i,j);
        
        speeds1 = inducedSpeeds(target, airfoils(1));
        speeds2 = inducedSpeeds(target, airfoils(2));
        speeds_symm1 = inducedSpeeds(target, symm_foil(1));
        speeds_symm2 = inducedSpeeds(target, symm_foil(2));
        
        [velSource1, velVortex1] = computeVelocityField(speeds1, gamma(1), q(:,1));
        [velSource2, velVortex2] = computeVelocityField(speeds2, gamma(2), q(:,2));
        [velSource_symm1, velVortex_symm1] = computeVelocityField(speeds_symm1, gamma(1), q(:,1));
        [velSource_symm2, velVortex_symm2] = computeVelocityField(speeds_symm2, gamma(2), q(:,2));

        vel_U(i,j) = velSource1(1)+velVortex1(1)+ ...
                     velSource2(1)+velVortex2(1)+ ...
                     velSource_symm1(1) - velVortex_symm1(1)+ ...
                     velSource_symm2(1) - velVortex_symm2(1)+ ...
                     U_infinity*cos(alpha);
        vel_V(i,j) = velSource1(2)+velVortex1(2)+ ...
                     velSource2(2)+velVortex2(2)+ ...
                     velSource_symm1(2) - velVortex_symm1(2)+ ...
                     velSource_symm2(2) - velVortex_symm2(2)+ ...
                     U_infinity*sin(alpha);
    end
end
vel_U(airfoil_mask==0) = NaN;
vel_V(airfoil_mask==0) = NaN;

vel_magnitude = sqrt(vel_U.^2 + vel_V.^2);

% streamlines
figure()
hold on
streamline(X,Y,vel_U,vel_V,linspace(xLim(1),xLim(1),nLines), linspace(yLim(1), yLim(2),nLines));
for i = 1:2
    plot(airfoils(i).x, airfoils(i).y,'k','Linewidth',1.2)
    fill(airfoils(i).x, airfoils(i).y,'w')
end
xlim(xLim)
ylim(yLim)
hold off

% colormap magnitude
figure();
hold on;
pcolor(X, Y, vel_magnitude);
shading interp;
colormap(jet);
colorbar;
for i = 1:2
    plot(airfoils(i).x, airfoils(i).y,'k','Linewidth',1.2)
    fill(airfoils(i).x, airfoils(i).y,'w')
end
title('Colormap velocity magnitude');
xlabel('X');
ylabel('Y');
xlim(xLim)
ylim(yLim)
grid off;
hold off;

% colormap U-component
figure();
hold on;
pcolor(X, Y, vel_U);
shading interp;
colormap(jet);
colorbar;
for i = 1:2
    plot(airfoils(i).x, airfoils(i).y,'k','Linewidth',1.2)
    fill(airfoils(i).x, airfoils(i).y,'w')
end
title('Colormap velocity x-component');
xlabel('X');
ylabel('Y');
xlim(xLim)
ylim(yLim)
grid off;
hold off;

% colormap V-component
figure();
hold on;
pcolor(X, Y, vel_V);
shading interp;
colormap(jet);
colorbar;
for i = 1:2
    plot(airfoils(i).x, airfoils(i).y,'k','Linewidth',1.2)
    fill(airfoils(i).x, airfoils(i).y,'w')
end
title('Colormap velocity y-component');
xlabel('X');
ylabel('Y');
xlim(xLim)
ylim(yLim)
grid off;
hold off;


%% Visualize results on coefficients

Cp = zeros(2, nPanels(1));
Cl = zeros(1,2);
Cm_LE = zeros(1,2);

for i = 1:2
    [velSource_i1, velVortex_i1] = computeVelocityField(velocities(i,1), gamma(1), q(:,1));
    [velSource_i2, velVortex_i2] = computeVelocityField(velocities(i,2), gamma(2), q(:,2));
    [velSource_isymm1, velVortex_isymm1] = computeVelocityField(symm_velocities(i,1), gamma(1), q(:,1));
    [velSource_isymm2, velVortex_isymm2] = computeVelocityField(symm_velocities(i,2), gamma(2), q(:,2));
    
    velSource = velSource_i1 + velSource_i2 + velSource_isymm1 + velSource_isymm2;
    velVortex = velVortex_i1 + velVortex_i2 - velVortex_isymm1 - velVortex_isymm2;
    
    [Cp(i,:),Cl(i),Cm_LE(i)] = computeAeroCoeffs(airfoils(i), velSource, velVortex, U_infinity, alpha, c(i));
    fprintf("Lift coefficient of airfoil %d: %d\n",i, Cl(i))
    fprintf("Moment coefficient w.r.t. leading edge of airfoil %d:: %d\n",i, Cm_LE(i))
    
    figure()
    hold on
    plot(airfoils(i).x_c(nPanels/2:-1:1),Cp(i, nPanels/2:-1:1),'--r','Linewidth',1.5)
    plot(airfoils(i).x_c(nPanels/2+1:nPanels),Cp(i,nPanels/2+1:nPanels),'-r','Linewidth',1.5)
    set (gca,'Ydir','reverse')
    l = legend('Bottom','Top');
    set(l,'FontSize',15)
    xlabel('chord','FontSize', 10,'fontweight','bold')
    ylabel('C_p','FontSize', 10,'fontweight','bold')
    grid on
    axis on
    hold off
end