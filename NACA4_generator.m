function airfoil = NACA4_generator(NACA4, c, AoA, n_Panels)
    % Generates 4-digit NACA airfoil panels with half-cosine spacing

    nChordPoints = n_Panels / 2;

    % Airfoil paramiters
    m = 1/100*floor(NACA4/1000);                %maximum camber of main element:
    p = floor((NACA4-100*m*1000)/100)/10;       %location of the maximum camber of main element:
    t = c*(NACA4 - p*10*100 - m*100*1000)/100;  %associated thickness:

    % Define coordinates with half-cosine method
    delta_theta = pi/(nChordPoints); 
    theta = 0:delta_theta:2*pi;
    x = c/2*(1 - cos(theta));
    x_top = x(1:nChordPoints+1);
    
    if m == 0 && p == 0 % symmetric airfoils
    
        x_bottom = x(nChordPoints+1:n_Panels);
        y_top = 5*t*(0.2969*sqrt(x_top/c) + (-0.1260)*(x_top/c) + ...
            (-0.3516)*(x_top/c).^2 + 0.2843*(x_top/c).^3 + ...
            (-0.1036)*(x_top/c).^4);
        y_bottom = -(y_top(end:-1:2));
        
        % Calculate extreme coordinates of the panel:
        x = horzcat(x_bottom,x_top);
        y = horzcat(y_bottom,y_top);
        y(1) = 0;
        y(end) = 0;
        
    else % not symmetric airfoils        
        
        y_top = 5*t*(0.2969*sqrt(x_top/c) + (-0.1260)*(x_top/c) + ...
            (-0.3516)*(x_top/c).^2 + 0.2843*(x_top/c).^3 + ...
            (-0.1036)*(x_top/c).^4); % Thickness distribution

        [~,idx] = min(abs(x_top - p*c));
        
        % For x < p*c
        y_camb1 = m*x_top(1:idx-1)/p^2.*(2*p - ...
            x_top(1:idx-1)/c);
        grad = 2*m/(p^2)*(p - x_top(1:idx-1)/c);
        kappa1 = atan(grad);
        
        %For x >= p*c
        y_camb2 = m*(c - x_top(idx:end))/((1 - p)^2).*(1 + ...
            x_top(idx:end)/c - 2*p);
        grad = 2*m/((1 - p)^2)*(p - x_top(idx:end)/c);
        kappa2 = atan(grad);
       
        y_camb = horzcat(y_camb1,y_camb2);
        kappa = horzcat(kappa1,kappa2);
       
        x_up = x_top - y_top.*sin(kappa); 
        x_down = fliplr(x_top + y_top.*sin(kappa));
        y_up = y_camb + y_top.*cos(kappa);
        y_down = fliplr(y_camb - y_top.*cos(kappa));
        
        %Calculate extreme coordinates of the panel:
        y = horzcat(y_down(1:end-1),y_up);
        x = horzcat(x_down(1:end-1),x_up);      
    end
    
    % Add angle of attack
    for i = 1:n_Panels+1
            x(i) = sum([cos(AoA), sin(AoA)].*[x(i), y(i)]);
            y(i) = sum([-sin(AoA), cos(AoA)].*[x(i), y(i)]);
    end
    
    for i =1:n_Panels
        dx(i) = x(i+1)-x(i);
        dy(i) = y(i+1)-y(i); 
        
        % control points:
        x_c(i) = dx(i)/2 + x(i);
        y_c(i) = dy(i)/2 + y(i);
        
        
        % panel length
        l_panel(i) = sqrt(dx(i).^2 + dy(i).^2);
        
        
        % orientation of the panels with respect to the global
        beta(i) = atan2(dy(i),dx(i));
    end
    
    airfoil.m = m;
    airfoil.p = p;
    airfoil.t = t;
    airfoil.x = x;
    airfoil.y = y;
    airfoil.x_c = x_c;
    airfoil.y_c = y_c;
    airfoil.l_panel = l_panel;
    airfoil.beta = beta;
    
end
