%% STEP (1): Defining the inputs

clc,clear,close

counter = -5;

for k = 1:30
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
    c(2) = 1;
    AoA(2) = counter; %Angle of attck of the airfoil
    AoA(2) = deg2rad(AoA(2));
    
    % --------- Airfoil paramiters: ---------
    %AIRFOIL 1:
    m(1) = 1/100*floor(airfoil(1)/1000); %maximum camber of main element:
    p(1) = floor((airfoil(1)-100*m(1)*1000)/100)/10; %location of the maximum camber of main element:
    t(1) = c(1)*(airfoil(1) - p(1)*10*100 - m(1)*100*1000)/100; %associated thickness:
    
    %AIRFOIL 2:
    m(2) = 1/100*floor(airfoil(2)/1000); %maximum camber of main element:
    p(2) = floor((airfoil(2)-100*m(2)*1000)/100)/10; %location of the maximum camber of main element:
    t(2) = c(2)*(airfoil(2) - p(2)*10*100 - m(2)*100*1000)/100; %associated thickness:
    
    
    % --------- Airfoil position: ---------
    h = 0.25; %heigth of the first airfoil
    z12 = 0.11; %height between two airfoils
    x12 = 1.2; %distance between two airfoils
    
    
    %% STEP (2): Buld the geometry:
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
            %Images 1:
            x_symm(i,:) = x(i,:);
            z_symm(i,:) = -z(i,:);
            x_c_symm(i,:) = x_c(i,:);
            z_c_symm(i,:) = -z_c(i,:);
            
        end
        
    end
    
    %% Step (2): Compute the induced velocity of source and vortex by panels-j w.r.t centroid-i
    
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
    A = zeros(n(1)+n(2)+2,n(1)+n(2)+2);
    
    A(1:n(1),1:n(1)) = A_s11; %ok
    A(1:n(1),(n(1)+1):(end-2)) = A_s12; %ok
    A((n(1)+1):(end-2),1:n(1)) = A_s21; %ok
    A((n(1)+1):(end-2),(n(1)+1):(end-2)) = A_s22;
    A(1:n(1),end-1:end) = [a_v11, a_v12];
    A((n(1)+1):end-2,end-1:end) = [a_v21, a_v22];
    A(end-1:end,1:end-2) = [c_s11, c_s12; c_s21, c_s22];
    A(end-1:end,end-1:end) = [c_v11, c_v12; c_v21, c_v22]; %ok
    b = [b_s1'; b_s2'; b_v1; b_v2];
    
    solution = A\b;
    
    q(:,1) = solution(1:n(1));
    q(:,2) = solution((n(1)+1):end-2);
    gamma(1) = solution(end-1);
    gamma(2) = solution(end);
    
    
    %% Step(5): Compute the flow field quantity:
    
    for i = 1:2
        [U_source1,U_vortex1] = source_vortex_velocity(x_c(i,:),z_c(i,:),x(1,:),z(1,:),beta(1,:),q(:,1),gamma(1),l_panel(1,:));
        [U_sourcem1,U_vortexm1] = source_vortex_velocity(x_c(i,:),z_c(i,:),x_symm(1,:),z_symm(1,:),beta(1,:),q(:,1),gamma(1),l_panel(1,:));
        [U_source2,U_vortex2] = source_vortex_velocity(x_c(i,:),z_c(i,:),x(2,:),z(2,:),beta(2,:),q(:,2),gamma(2),l_panel(2,:));
        [U_sourcem2,U_vortexm2] = source_vortex_velocity(x_c(i,:),z_c(i,:),x_symm(2,:),z_symm(2,:),beta(2,:),q(:,2),gamma(2),l_panel(2,:));
        U_source = U_source1+U_source2+U_sourcem1+U_sourcem2;
        U_vortex = U_vortex1+U_vortex2-U_vortexm1-U_vortexm2;
        
        [~,Cl(k,i),~] = flow_paramiters(beta(i,:),n(i),l_panel(i,:),x_c(i,:),z_c(i,:),c(i),U_source,U_vortex,U_infinity,alpha);
        
    end
    
    counter = counter + 1;
end


figure()
hold on
plot(-5:1:24,Cl(:,1),'-*b','linewidth',1)
grid on
xlabel('α_{tail}','FontSize', 10,'fontweight','bold')
ylabel('C_{l,wilng}','FontSize', 10,'fontweight','bold')
title('C_{l,wing} variation with ground effect','FontSize', 15)
hold off
%saveas(gcf, 'Find angle tail ground effect','png')


