%% Set-up inputs
clc,clear,close

% --------- Airfoil Inputs: ---------
%AIRFOIL:
nPanels = 50;      % Number of panels
NACA4 = 0012;       % NACA airfoil denomination ( 4 digit: XXXX )
c = 1;              %chord of main element
AoA = deg2rad(2.5);   %Angle of attack (deg)

% --------- Flow paramiters Inputs: ---------
U_infinity = 1; %Velocity of the streamflow
alpha = 0;      %Angle of the stream flow


%% Build airfoil

% Generate the airfoil panels
airfoil = NACA4_generator(NACA4, c, AoA, nPanels);

% show airfoil geometry:
figure()
hold on
plot(airfoil.x, airfoil.y)
plot(airfoil.x_c, airfoil.y_c,'x')
grid on
axis equal
legend('panels','control points')
hold off

%% Evaluation flow paramiters

% Compute coefficients due to panels induced velocities (source and vortex):
velocities = inducedSpeeds(airfoil, airfoil);
coeffs = matrixBlocks(airfoil, velocities, U_infinity, alpha);

% Build and solve the linear sistem 
A = zeros(nPanels+1,nPanels+1);

A(1:nPanels,1:nPanels) = coeffs.A_s;
A(1:nPanels,end) = coeffs.a_v;
A(end,1:nPanels) = coeffs.c_s;
A(end,end) = coeffs.c_v;
b = [coeffs.b_s; coeffs.b_v];

solution = A\b;

q = solution(1:nPanels);
gamma = solution(end);

%% Visualize results on velocity
% ----- Plot settings: --------
xLim=[-0.5, 2.5];
yLim=[-1,1];
mSize = 100;
nLines = 20;

% discretization
xm = linspace(xLim(1), xLim(2), mSize);
ym = linspace(yLim(1), yLim(2), mSize);
[X,Y] = meshgrid(xm, ym);

airfoil_mask = ~inpolygon(X, Y, airfoil.x, airfoil.y);
airfoil_mask(airfoil_mask==0) = 0;
vel_U = zeros(length(xm),length(ym));
vel_V = zeros(length(xm),length(ym));

for i = 1:length(xm)
    for j = 1:length(ym)

        target.x_c = X(i,j);
        target.y_c = Y(i,j);
        speeds = inducedSpeeds(target, airfoil);
        [velSource, velVortex] = computeVelocityField(speeds, gamma, q);
        vel_U(i,j) = velSource(1)+velVortex(1)+U_infinity*cos(alpha);
        vel_V(i,j) = velSource(2)+velVortex(2)+U_infinity*sin(alpha);

    end
end
vel_U(airfoil_mask==0) = NaN;
vel_V(airfoil_mask==0) = NaN;

vel_magnitude = sqrt(vel_U.^2 + vel_V.^2);

% streamlines
figure()
hold on
streamline(X,Y,vel_U,vel_V,linspace(xLim(1),xLim(1),nLines), linspace(yLim(1), yLim(2),nLines));
plot(airfoil.x, airfoil.y,'k','Linewidth',1.2)
fill(airfoil.x,airfoil.y,'w')
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
plot(airfoil.x, airfoil.y,'k','Linewidth',1.2)
fill(airfoil.x,airfoil.y,'w')
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
plot(airfoil.x, airfoil.y,'k','Linewidth',1.2)
fill(airfoil.x,airfoil.y,'w')
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
plot(airfoil.x, airfoil.y,'k','Linewidth',1.2)
fill(airfoil.x,airfoil.y,'w')
title('Colormap velocity y-component');
xlabel('X');
ylabel('Y');
xlim(xLim)
ylim(yLim)
grid off;
hold off;

%% Visualize results on coefficients
[velSource, velVortex] = computeVelocityField(velocities, gamma, q);
[Cp,Cl,Cm_LE] = computeAeroCoeffs(airfoil, velSource, velVortex, U_infinity, alpha, c);
fprintf("Lift coefficient: %d\n",Cl)
fprintf("Moment coefficient w.r.t. leading edge: %d\n",Cm_LE)

if NACA4 == 0012 & (AoA == deg2rad(2.5) || AoA == deg2rad(0))
    fprintf('\nValidation data available, check the plots...\n')

    if AoA == deg2rad(2.5)
        XFLR5 = load('0012_2.5deg.txt');
    elseif AoA == deg2rad(0)
        XFLR5 = load('0012_0deg.txt');
    else
        warning('\n\nNo available data for validation!!!\n\n')
    end

    Cp_xef = XFLR5(:,2);
    x_xef = XFLR5(:,1);
    
    n = nPanels;
    figure()
    hold on
    plot(airfoil.x_c(nPanels/2:-1:1),Cp(nPanels/2:-1:1),'--b','linewidth', 1.5)
    plot(airfoil.x_c(nPanels/2+1:n),Cp(nPanels/2+1:n),'-b','linewidth', 1.5)
    plot(x_xef(50:-1:1,1),Cp_xef(50:99,1),'*k',x_xef(50:-1:1,1),Cp_xef(50:-1:1,1),'*k')
    plot(linspace(0,1),zeros(100,1),'k')
    set (gca,'Ydir','reverse')
    xlim([0,c])
    title('Pressure coefficient comparison','FontSize', 15)
    l = legend('Bottom','Top', 'XFLR 5','Location','southeast');
    set(l,'FontSize',15)
    xlabel('chord')
    ylabel('C_p')
    grid on
    hold off
else
    fprintf('\nValidation data NOT available for this arifoil and angle of attack.\n')

    figure()
    hold on
    plot(airfoil.x_c(nPanels/2:-1:1),Cp(nPanels/2:-1:1),'--r','Linewidth',1.5)
    plot(airfoil.x_c(nPanels/2+1:nPanels),Cp(nPanels/2+1:nPanels),'-r','Linewidth',1.5)
    set (gca,'Ydir','reverse')
    l = legend('Bottom','Top');
    set(l,'FontSize',15)
    xlabel('chord','FontSize', 10,'fontweight','bold')
    ylabel('C_p','FontSize', 10,'fontweight','bold')
    grid on
    axis on
    hold off
end