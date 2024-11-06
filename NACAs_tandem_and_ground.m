%Script for evaluation of Cp, Cl and Cm_LE for a tandem configuration
%with the ground effect.

%INPUTS:
%Number of panels, Chord, Angle of attack (also possible to set a range of
%Aoa for Cl-alpha curve), U_infinity, Alpha, Distances between airfoils and
%from the ground

%OUTPUTS:
%Cp curve of each airfoil, Cl (assuming a set of AoA for one of the two 
%airfoils will provide also the Cl-alpha curve), Cm_LE

%N.B. --> If you want to save the plot de-indent the save command


%% Set-up inputs
clc,clear,close


% --------- Flow paramiters Inputs: ---------
U_infinity = 1; %Velocity of the streamflow
alpha = 0; %Angle of the stream flow
alpha = deg2rad(alpha);



% --------- Airfoil Inputs: ---------
%AIRFOIL 1:
N_panel(1) = 50; % Number of point you divide the chord (the number of panels is 2*N_panel)
n(1) = N_panel(1)*2;
airfoil(1) = 3310; % NACA airfoil denomination ( 4 digit: XXXX )
c(1) = 1; %chord of main element
AoA(1) = 2.5; %Angle of attck of the airfoil
AoA(1) = deg2rad(AoA(1));

%AIRFOIL 2:
N_panel(2) = 50; 
n(2) = N_panel(2)*2;
airfoil(2) = 0012; % NACA airfoil denomination ( 4 digit: XXXX )
c(2) = 0.4;
AoA(2) = 5; %Angle of attck of the airfoil
AoA(2) = deg2rad(AoA(2));



% --------- Airfoil position: ---------
h = 0.8; %heigth of the first airfoil
x12 = 1.2; %distance between two airfoils
z12 = -0.5; %height between two airfoils




%% STEP (1): Buld the geometry:
% --------- Airfoil paramiters: ---------
%AIRFOIL 1:
m(1) = 1/100*floor(airfoil(1)/1000); %maximum camber of main element:
p(1) = floor((airfoil(1)-100*m(1)*1000)/100)/10; %location of the maximum camber of main element:
t(1) = c(1)*(airfoil(1) - p(1)*10*100 - m(1)*100*1000)/100; %associated thickness:

%AIRFOIL 2:
m(2) = 1/100*floor(airfoil(2)/1000); %maximum camber of main element:
p(2) = floor((airfoil(2)-100*m(2)*1000)/100)/10; %location of the maximum camber of main element:
t(2) = c(2)*(airfoil(2) - p(2)*10*100 - m(2)*100*1000)/100; %associated thickness:


%Initialization vectors
N = max(n);
x = zeros(2,N+1);
z = zeros(2,N+1);
x_c = zeros(2,N);
z_c = zeros(2,N);
l_panel_ = zeros(2,N);
beta = zeros(2,N);
x_symm = zeros(2,N+1);
z_symm = zeros(2,N+1);
x_c_symm = zeros(2,N);
z_c_symm = zeros(2,N);

for i = 1:2
    
    [x(i,1:(n(i)+1)),z(i,1:(n(i)+1)),x_c(i,1:n(i)),z_c(i,1:n(i)),l_panel(i,1:n(i)),beta(i,1:n(i))] = Panel_calculator(c(i),t(i),N_panel(i),m(i),p(i),AoA(i));
    
    if i == 1
        %AIRFOIL 1 position:
        z(i,:) = z(i,:) + h;
        z_c(i,:) = z_c(i,:) + h;
        %Images 1:
        x_symm(i,:) = x(i,:);
        z_symm(i,:) = -z(i,:);
        x_c_symm(i,:) = x_c(i,:);
        z_c_symm(i,:) = -z_c(i,:);
 
    elseif i == 2
        %AIRFOIL 2 position:
        z(i,:) = z(i,:) + (h + z12);
        x(i,:) = x(i,:) + x12;
        z_c(i,:) = z_c(i,:) + (h + z12);
        x_c(i,:) = x_c(i,:) + x12;
        %Images 2:
        x_symm(i,:) = x(i,:);
        z_symm(i,:) = -z(i,:);
        x_c_symm(i,:) = x_c(i,:);
        z_c_symm(i,:) = -z_c(i,:);
        
    end
    
end



% --------- Show the aerofoils SINGULARLY: ---------
figure()

subplot(211)
hold on
plot(x(1,1:(n(1)+1)),z(1,1:(n(1)+1)),'k','Linewidth',1.2)
plot(x_c(1,1:n(1)),z_c(1,1:n(1)),'x')
title('Airfoil 1 and extremes')
legend('panels','control points')
hold off

subplot(212)
hold on
plot(x(2,1:(n(2)+1)),z(2,1:(n(2)+1)),'k','Linewidth',1.2)
plot(x_c(2,1:n(2)),z_c(2,1:n(2)),'x')
title('Airfoil 2 and extremes')
legend('panels','control points')
hold off
% saveas(gcf, 'Tandem and ground geometry singularly','png')


% --------- Show the aerofoils and mirrored TOGHTTHER: ---------
figure()
hold on
plot(x(1,1:(n(1)+1)),z(1,1:(n(1)+1)),'k','Linewidth',1.2)
plot(x_symm(1,1:(n(1)+1)),z_symm(1,1:(n(1)+1)),'--k','Linewidth',0.6)
plot(x(2,1:(n(2)+1)),z(2,1:(n(2)+1)),'k','Linewidth',1.2)
plot(x_symm(2,1:(n(2)+1)),z_symm(2,1:(n(2)+1)),'--k','Linewidth',0.6)
plot(linspace(-1,2),zeros(1,100),'k')
title('Airfoil configuration with ground effect')
% legend('Real airfoils','Mirrored')
axis equal
grid off
hold off
% saveas(gcf, 'Tandem and ground geometry and mirrored','png')




%% Step (2): Evaluate the influence coefficients for the linear system
        
[A_s11,a_v11,c_s11,c_v11,b_s1,b_v1] = Influence_coefficient(beta(1,:),x_c(1,:),z_c(1,:),x(1,:),z(1,:),l_panel(1,:),U_infinity,alpha); %velocity of airfoil 1 induced by airfoil 1 
[A_s_symm,a_v_symm,c_s_symm,c_v_symm] = Influence_coefficient(beta(1,:),x_c(1,:),z_c(1,:),x_symm(1,:),z_symm(1,:),l_panel(1,:),U_infinity,alpha); %velocity of airfoil 1 induced by airfoil 1 
A_s11 = A_s11 + A_s_symm;
a_v11 = a_v11 - a_v_symm;
c_s11 = c_s11 + c_s_symm;
c_v11 = c_v11 - c_v_symm;

[A_s12,a_v12,c_s12,c_v12] = Influence_coefficient(beta(2,:),x_c(1,:),z_c(1,:),x(2,:),z(2,:),l_panel(2,:),U_infinity,alpha); %velocity of airfoil 1 induced by airfoil 2 
[A_s_symm,a_v_symm,c_s_symm,c_v_symm] = Influence_coefficient(beta(2,:),x_c(1,:),z_c(1,:),x_symm(2,:),z_symm(2,:),l_panel(2,:),U_infinity,alpha); %velocity of airfoil 1 induced by airfoil 2 
A_s12 = A_s12 + A_s_symm;
a_v12 = a_v12 - a_v_symm;
c_s12 = c_s12 + c_s_symm;
c_v12 = c_v12 - c_v_symm;

[A_s22,a_v22,c_s22,c_v22,b_s2,b_v2] = Influence_coefficient(beta(2,:),x_c(2,:),z_c(2,:),x(2,:),z(2,:),l_panel(2,:),U_infinity,alpha); %velocity of airfoil 2 induced by airfoil 2 
[A_s_symm,a_v_symm,c_s_symm,c_v_symm] = Influence_coefficient(beta(2,:),x_c(2,:),z_c(2,:),x_symm(2,:),z_symm(2,:),l_panel(2,:),U_infinity,alpha); %velocity of airfoil 2 induced by airfoil 2 
A_s22 = A_s22 + A_s_symm;
a_v22 = a_v22 - a_v_symm;
c_s22 = c_s22 + c_s_symm;
c_v22 = c_v22 - c_v_symm;

[A_s21,a_v21,c_s21,c_v21] = Influence_coefficient(beta(1,:),x_c(2,:),z_c(2,:),x(1,:),z(1,:),l_panel(1,:),U_infinity,alpha); %velocity of airfoil 2 induced by airfoil 1 
[A_s_symm,a_v_symm,c_s_symm,c_v_symm] = Influence_coefficient(beta(1,:),x_c(2,:),z_c(2,:),x_symm(1,:),z_symm(1,:),l_panel(1,:),U_infinity,alpha); %velocity of airfoil 2 induced by airfoil 1 
A_s21 = A_s21 + A_s_symm;
a_v21 = a_v21 - a_v_symm;
c_s21 = c_s21 + c_s_symm;
c_v21 = c_v21 - c_v_symm;


%% Step (3)&(4): Build and solve the linear sistem for the strength q_i and the vorticity gamma:
%Initialize matrix:
A = zeros(n(1)+n(2)+2,n(1)+n(2)+2);

%Build up the system:
A(1:n(1),1:n(1)) = A_s11; %ok
A(1:n(1),(n(1)+1):(end-2)) = A_s12; %ok
A((n(1)+1):(end-2),1:n(1)) = A_s21; %ok
A((n(1)+1):(end-2),(n(1)+1):(end-2)) = A_s22; %ok
A(1:n(1),end-1:end) = [a_v11, a_v12]; %ok
A((n(1)+1):end-2,end-1:end) = [a_v21, a_v22]; %ok
A(end-1:end,1:end-2) = [c_s11, c_s12; c_s21, c_s22]; %ok
A(end-1:end,end-1:end) = [c_v11, c_v12; c_v21, c_v22]; %ok
b = [b_s1'; b_s2'; b_v1; b_v2]; %ok

%Define the solution:
solution = A\b;
q(:,1) = solution(1:n(1)); %Strength front Airfoil
q(:,2) = solution((n(1)+1):end-2); %Strength rear Airfoil
gamma(1) = solution(end-1); %Front Airfoil
gamma(2) = solution(end); %Rear Airfoil
gamma(2) = solution(end);



%% Step(5): Compute the flow field quantity:

                      
for i = 1:2
    [U_source1,U_vortex1] = source_vortex_velocity(x_c(i,:),z_c(i,:),x(1,:),z(1,:),beta(1,:),q(:,1),gamma(1),l_panel(1,:));
    [U_sourcem1,U_vortexm1] = source_vortex_velocity(x_c(i,:),z_c(i,:),x_symm(1,:),z_symm(1,:),beta(1,:),q(:,1),gamma(1),l_panel(1,:));
    [U_source2,U_vortex2] = source_vortex_velocity(x_c(i,:),z_c(i,:),x(2,:),z(2,:),beta(2,:),q(:,2),gamma(2),l_panel(2,:));
    [U_sourcem2,U_vortexm2] = source_vortex_velocity(x_c(i,:),z_c(i,:),x_symm(2,:),z_symm(2,:),beta(2,:),q(:,2),gamma(2),l_panel(2,:));
    U_source = U_source1+U_source2+U_sourcem1+U_sourcem2;
    U_vortex = U_vortex1+U_vortex2-U_vortexm1-U_vortexm2;

    [Cp(i,:),Cl(i),Cm_LE(i)] = flow_paramiters(beta(i,:),n(i),l_panel(i,:),x_c(i,:),z_c(i,:),c(i),U_source,U_vortex,U_infinity,alpha);

end




%% Plot results:


% --------- Show Cp: ---------
%Togheter
figure()
hold on
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]); % ingrandisci a tutto schermo
plot(x_c(1,n(1)/2:-1:1),Cp(1,n(1)/2:-1:1),'--b','linewidth',1.5)
plot(x_c(1,n(1)/2+1:n(1)),Cp(1,n(1)/2+1:n(1)),'-b','linewidth',1.5)
plot(x_c(2,n(2)/2:-1:1)-x12,Cp(2,n(2)/2:-1:1),'--r','linewidth',1.5)
plot(x_c(2,n(2)/2+1:n(2))-x12,Cp(2,n(2)/2+1:n(2)),'-r','linewidth',1.5)
set (gca,'Ydir','reverse')
title('C_p distribution of tandem configuration and ground','FontSize',15)
l = legend('Bottom wing','Top wing','Bottom tail','Top tail');
set(l,'FontSize',15)
xlabel('chord','FontSize', 10,'fontweight','bold')
ylabel('C_p','FontSize', 10,'fontweight','bold')
grid on
axis on
% saveas(gcf, 'Cp togheter tandem e ground','png')


%Separately
figure()
hold on
plot(x_c(1,n(1)/2:-1:1),Cp(1,n(1)/2:-1:1),'--b',x_c(1,n(1)/2+1:n(1)),Cp(1,n(1)/2+1:n(1)),'-b')
plot(x_c(2,n(2)/2:-1:1),Cp(2,n(2)/2:-1:1),'--r',x_c(2,n(2)/2+1:n(2)),Cp(2,n(2)/2+1:n(2)),'-r')
plot(x(1,1:(n(1)+1)),z(1,1:(n(1)+1)),'k','Linewidth',1.2)
plot(x(2,1:(n(2)+1)),z(2,1:(n(2)+1)),'k','Linewidth',1.2)
set (gca,'Ydir','reverse')
title('Pressure coefficient with ground','FontSize', 15)
legend('Bottom wing','Top wing','Bottom tail','Top tail')
xlabel('chord','FontSize', 10,'fontweight','bold')
ylabel('C_p','FontSize', 10,'fontweight','bold')
grid on
axis on
% saveas(gcf, 'Cp singularly tandem e ground','png')



%% Step(6): Post processing plotting the results:


%Mesh plot:
d = 100; %finezza mesh
xm = linspace(-0.2,2.5,d); %x-dimension
zm = linspace(0,1.2,d); %y-dmension
[X,Y] = meshgrid(xm,zm);

airfoil_mask = ~inpolygon(X, Y, x(1,:), z(1,:));
airfoil_mask1 = ~inpolygon(X, Y, x(2,:), z(2,:));

for i = 1:d
    for j = 1:d
        if airfoil_mask1(i,j) == 0
            airfoil_mask(i,j) = 0;
        end
    end
end


for i = 1:length(xm)
    for j = 1:length(zm)
        if airfoil_mask(i,j) == 1
            [U_source1,U_vortex1] = source_vortex_velocity(X(i,j),Y(i,j),x(1,:),z(1,:),beta(1,:),q(:,1),gamma(1),l_panel(1,:));
            [U_sourcem1,U_vortexm1] = source_vortex_velocity(X(i,j),Y(i,j),x_symm(1,:),z_symm(1,:),beta(1,:),q(:,1),gamma(1),l_panel(1,:));
            [U_source2,U_vortex2] = source_vortex_velocity(X(i,j),Y(i,j),x(2,:),z(2,:),beta(2,:),q(:,2),gamma(2),l_panel(2,:));
            [U_sourcem2,U_vortexm2] = source_vortex_velocity(X(i,j),Y(i,j),x_symm(2,:),z_symm(2,:),beta(2,:),q(:,2),gamma(2),l_panel(2,:));
            vel_U(i,j) = U_source1(1)+U_source2(1)+U_sourcem1(1)+U_sourcem2(1)+U_vortex1(1)+U_vortex2(1)-U_vortexm1(1)-U_vortexm2(1)+U_infinity*cos(alpha);
            vel_V(i,j) = U_source1(2)+U_source2(2)+U_sourcem1(2)+U_sourcem2(2)+U_vortex1(2)+U_vortex2(2)-U_vortexm1(2)-U_vortexm2(2)+U_infinity*sin(alpha);
                      
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
plot(x(2,1:(n(2)+1)),z(2,1:(n(2)+1)),'k','Linewidth',1.2)
fill(x(2,1:(n(2)+1)),z(2,1:(n(2)+1)),'w')
plot(linspace(xm(1),xm(end)),zeros(1,100),'k','Linewidth',3)
axis equal
grid off
hold off
%saveas(gcf, strcat('Streamlines Airfoil ground',num2str(airfoil(1)),num2str(airfoil(2))),'png')



figure();
hold on;
pcolor(X, Y, vel_magnitude);
shading interp;
colormap(jet);
colorbar;
plot(x(1,1:(n(1)+1)),z(1,1:(n(1)+1)),'k','Linewidth',1.2)
fill(x(1,1:(n(1)+1)),z(1,1:(n(1)+1)),'w')
plot(x(2,1:(n(2)+1)),z(2,1:(n(2)+1)),'k','Linewidth',1.2)
fill(x(2,1:(n(2)+1)),z(2,1:(n(2)+1)),'w')
axis equal;
title('Mappa colorata della magnitudine della velocità');
xlabel('X');
ylabel('Y');
grid off;
hold off;
%saveas(gcf, strcat('Color map Airfoil ground',num2str(airfoil(1)),num2str(airfoil(2))),'png')