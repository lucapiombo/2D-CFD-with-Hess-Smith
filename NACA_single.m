%Script for evaluation of Cp, Cl and Cm_LE for a single airfoil. 

%INPUTS:
%Number of panels, Chord, Angle of attack (also possible to set a range of
%Aoa for Cl-alpha curve), U_infinity, Alpha

%OUTPUTS:
%Cp curve of the single airfoil, Cl (assuming a set of AoA also the
%Cl-alpha curve), Cm_LE
%N.B. --> If you want to save the plot de-indent the save command

%% Set-up inputs
clc,clear,close


% --------- Airfoil Inputs: ---------
%AIRFOIL:
N_panel = 50; % Number of point you divide the chord (the number of panels is 2*N_panel)
airfoil = 0012; % NACA airfoil denomination ( 4 digit: XXXX )
c = 1; %chord of main element
AoA = 2.5; %Range of your angle of attack
AoA = deg2rad(AoA);


% --------- Flow paramiters Inputs: ---------
U_infinity = 1; %Velocity of the streamflow
alpha = 0; %Angle of the stream flow




%% Evaluation flow paramiters with a cycle

% --------- Airfoil paramiters: ---------
%AIRFOIL 1:
m = 1/100*floor(airfoil/1000); %maximum camber of main element:
p = floor((airfoil-100*m*1000)/100)/10; %location of the maximum camber of main element:
t = c*(airfoil - p*10*100 - m*100*1000)/100; %associated thickness:
n = N_panel*2;

%Initialize vectors:
Cp = zeros(length(AoA),n);
Cl = zeros(length(AoA),1);
Cm_LE = zeros(length(AoA),1);


for i = 1:length(AoA)  
    
    %--------- Evaluate extrema and centroides of panels: ---------
    [x,z,x_c(i,:),z_c,l_panel,beta] = Panel_calculator(c,t,N_panel,m,p,AoA(i));
    t_hat = [cos(beta)', sin(beta)'];
    n_hat = [-sin(beta)',cos(beta)'];
    
    
    %--------- Compute the induced velocity of source and vortex by
    %panels-j w.r.t centroid-i: ---------
    
    [A_s,a_v,c_s,c_v,b_s,b_v] = Influence_coefficient(beta,x_c(i,:),z_c,x,z,l_panel,U_infinity,alpha); %velocity of airfoil 1 induced by airfoil 1 

    
    
    %--------- Build and solve the linear sistem: ---------    
    A = zeros(n+1,n+1);
    
    A(1:n,1:n) = A_s;
    A(1:n,end) = a_v;
    A(end,1:n) = c_s';
    A(end,end) = c_v;
    b = [b_s'; b_v];
    
    solution = A\b;
    
    q = solution(1:n);
    gamma = solution(end);
    
    
    %--------- Compute the flow field quantity: ---------
    [U_source,U_vortex] = source_vortex_velocity(x_c(i,:),z_c,x,z,beta,q,gamma,l_panel);
    [Cp(i,:),Cl(i),Cm_LE(i)] = flow_paramiters(beta,n,l_panel,x_c(i,:),z_c,c,U_source,U_vortex,U_infinity,alpha);
        
    % --------- Show the aerofoil SINGULARLY: ---------
    figure()
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]); % ingrandisci a tutto schermo
    hold on
    plot(x,z,'k','Linewidth',1.5)
%     plot(x_c,z_c,'x')
    title(strcat('NACA',{' '},num2str(airfoil), {' '},'for', {' '},'α =',{' '},num2str(rad2deg(AoA(i)))),'FontSize', 15)
    grid on
    xlim([-0.2,1.2])
    ylim([-0.3,0.3])
%     legend('panels','control points')
    hold off
%     saveas(gcf, strcat('Airfoil', num2str(airfoil)),'png')
end



%% Plot results:


% --------- Show Cp: ---------
for i = 1:length(AoA)
    figure()
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]); % ingrandisci a tutto schermo
    hold on
    plot(x_c(i,n/2:-1:1),Cp(i,n/2:-1:1),'--r','Linewidth',1.5)
    plot(x_c(i,n/2+1:n),Cp(i,n/2+1:n),'-r','Linewidth',1.5)
    set (gca,'Ydir','reverse')
    l = legend('Bottom','Top');
    set(l,'FontSize',15)
    title(strcat('Pressure coefficient NACA', {' '},num2str(airfoil), {' '},'at', {' '},'α =',{' '},num2str(rad2deg(AoA(i)))),'FontSize', 15)
    xlabel('chord','FontSize', 10,'fontweight','bold')
    ylabel('C_p','FontSize', 10,'fontweight','bold')
    grid on
    axis on
    hold off
end
% saveas(gcf, strcat('Cp Airfoil',num2str(airfoil)),'png')




% --------- Show Cl (only for more than 1 alpha): ---------
if length(AoA)>1
    figure()
    hold on
    
    plot(rad2deg(AoA),Cl,'-k')
    xlabel('alpha','FontSize', 10,'fontweight','bold')
    ylabel('C_l','FontSize', 10,'fontweight','bold')
    title(strcat('Lift coefficient NACA', {' 00'},num2str(airfoil)),'FontSize', 15)
    axis on
    grid on
    
    hold off
%     saveas(gcf, strcat('Cl Airfoil 00',num2str(airfoil)),'png')
end


%% Step(6): Post processing plotting the results:

XFLR5 = load('0012_2.5deg.txt');
Cp_xef = XFLR5(:,2);
x_xef = XFLR5(:,1);

figure()
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]); % ingrandisci a tutto schermo
hold on
plot(x_c(n/2:-1:1),Cp(n/2:-1:1),'--b','linewidth', 1.5)
plot(x_c(n/2+1:n),Cp(n/2+1:n),'-b','linewidth', 1.5)
plot(x_xef(n/2:-1:1),Cp_xef(50:99,1),'*k',x_xef(n/2:n-1),Cp_xef(50:-1:1,1),'*k')
plot(linspace(0,1),zeros(100,1),'k')
set (gca,'Ydir','reverse')
xlim([0,1])
title('Pressure coefficient comparison','FontSize', 15)
l = legend('Bottom','Top', 'XFLR 5','Location','southeast');
set(l,'FontSize',15)
xlabel('chord')
ylabel('C_p')
grid on
hold off
% saveas(gcf, strcat('Cp convalida Airfoil 25 and ',num2str(airfoil)),'png')


%% Step(7): Post processing plotting the results:

%Mesh plot:
d = 100; %finezza mesh
xm = linspace(-0.5,2.5,d); %x-dimension
zm = linspace(-1,1,d); %y-dmension
[X,Y] = meshgrid(xm,zm);

airfoil_mask = ~inpolygon(X, Y, x(1,:), z(1,:));

for i = 1:d
    for j = 1:d
        if airfoil_mask(i,j) == 0
            airfoil_mask(i,j) = 0;
        end
    end
end


for i = 1:length(xm)
    for j = 1:length(zm)
        if airfoil_mask(i,j) == 1
            [U_source1,U_vortex1] = source_vortex_velocity(X(i,j),Y(i,j),x(1,:),z(1,:),beta(1,:),q(:,1),gamma(1),l_panel(1,:));
            vel_U(i,j) = U_source1(1)+U_vortex1(1)+U_infinity*cos(alpha);
            vel_V(i,j) = U_source1(2)+U_vortex1(2)+U_infinity*sin(alpha);
                      
        else
            vel_U(i,j) = 0;
            vel_V(i,j) = 0;
            vel(i,j) = 0;        
        end
    end
end
vel_magnitude = sqrt(vel_U.^2 + vel_V.^2);

% Create a figure
figure()
hold on
streamline(X,Y,vel_U,vel_V,linspace(-0.2,-0.2,d), zm);
plot(x(1,1:(n(1)+1)),z(1,1:(n(1)+1)),'k','Linewidth',1.2)
fill(x(1,1:(n(1)+1)),z(1,1:(n(1)+1)),'w')
xlim([xm(1),xm(end)])
ylim([-0.3,0.3])
axis equal
hold off
%saveas(gcf, strcat('Streamlines Airfoil',num2str(airfoil(1)),num2str(airfoil(2))),'png')

figure();
hold on;
pcolor(X, Y, vel_magnitude);
shading interp;
colormap(jet);
colorbar;
plot(x(1,1:(n(1)+1)),z(1,1:(n(1)+1)),'k','Linewidth',1.2)
fill(x(1,1:(n(1)+1)),z(1,1:(n(1)+1)),'w')
axis equal;
title('Mappa colorata della magnitudine della velocità');
xlabel('X');
ylabel('Y');
grid off;
hold off;
%saveas(gcf, strcat('Color map Airfoil',num2str(airfoil(1)),num2str(airfoil(2))),'png')



