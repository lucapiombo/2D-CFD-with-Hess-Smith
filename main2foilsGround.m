%% Set-up inputs
clc,clear,close

% --------- Flow paramiters Inputs: ---------
U_infinity = 1; %Velocity of the streamflow
alpha = 0;      %Angle of the stream flow


% --------- Airfoil Inputs: ---------
%AIRFOIL 1:
nPanels(1) = 50;       % Number of panels
NACA4(1) = 3310;        % NACA airfoil denomination ( 4 digit: XXXX )
c(1) = 1;               %chord of main element
AoA(1) = deg2rad(0);    %Angle of attack (deg)

%AIRFOIL 2:
nPanels(2) = 40;       % Number of panels
NACA4(2) = 3312;        % NACA airfoil denomination ( 4 digit: XXXX )
c(2) = 0.4;             %chord of main element
AoA(2) = deg2rad(5);    %Angle of attack (deg)


% --------- Airfoil position: ---------
h = 0.8;    %heigth of the first airfoil
x12 = 1.1;  %distance between two airfoils
y12 = 0.5;  %height between two airfoils

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
        fprintf(i)
        airfoils(i).y = airfoils(i).y + (y12 + h);
        airfoils(i).x = airfoils(i).x + x12;
        airfoils(i).y_c = airfoils(i).y_c + (y12 + h);
        airfoils(i).x_c = airfoils(i).x_c + x12; 
    end

    %Images foils:
    airfoils(i).x_symm = airfoils(i).x;
    airfoils(i).y_symm = -airfoils(i).y;
    airfoils(i).x_c_symm = airfoils(i).x_c;
    airfoils(i).y_c_symm = -airfoils(i).y_c;

    % show airfoil geometry:
    hold on
    plot(airfoils(i).x, airfoils(i).y)
    plot(airfoils(i).x_c, airfoils(i).y_c,'x')
    plot(airfoils(i).x_symm, airfoils(i).y_symm)
    plot(airfoils(i).x_c_symm, airfoils(i).y_c_symm,'x')
end
grid on
axis equal
legend('panels','control points')
hold off


%% Evaluation flow paramiters

% Compute coefficients due to panels induced velocities (source and vortex):
for i = 1:2 % targetFoil
    for j= 1:2 % triggerFoil
        velocities(i,j) = inducedSpeeds(airfoils(i), airfoils(j));
        coeffs(i,j) = matrixBlocks(airfoils(j), velocities(i,j), U_infinity, alpha);
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
b = [coeffs(1,1).b_s; coeffs(2,2).b_s; coeffs(1,1).b_v; coeffs(2,2).b_v];

%Define the solution:
solution = A\b;
q(:,1) = solution(1:nPanels(1));            %Strength front Airfoil
q(:,2) = solution((nPanels(1)+1):end-2);    %Strength rear Airfoil
gamma(1) = solution(end-1);                 %Front Airfoil
gamma(2) = solution(end);                   %Rear Airfoil


%% Visualize results on velocity
% ----- Plot settings: --------
xLim=[-0.5, 2.5];
yLim=[-1,1];
mSize = 500;
nLines = 20;

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
        [velSource1, velVortex1] = computeVelocityField(speeds1, gamma(1), q(:,1));
        [velSource2, velVortex2] = computeVelocityField(speeds2, gamma(2), q(:,2));

        vel_U(i,j) = velSource1(1)+velVortex1(1)+velSource2(1)+velVortex2(1)+U_infinity*cos(alpha);
        vel_V(i,j) = velSource1(2)+velVortex1(2)+velSource2(2)+velVortex2(2)+U_infinity*sin(alpha);

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

[velSource11, velVortex11] = computeVelocityField(velocities(1,1), gamma(1), q(:,1));
[velSource12, velVortex12] = computeVelocityField(velocities(1,2), gamma(2), q(:,2));
[velSource21, velVortex21] = computeVelocityField(velocities(2,1), gamma(1), q(:,1));
[velSource22, velVortex22] = computeVelocityField(velocities(2,2), gamma(2), q(:,2));
velSource = velSource11 + velSource12 + velSource21 + velSource22;
velVortex = velVortex11 + velVortex12 + velVortex21 + velVortex22;
for i = 1:2
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