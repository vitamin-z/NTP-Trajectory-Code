function [R,V,t] = prop_sun(r_o, v_o,theta)

  %Constants  
    G = 6.67408*10^-20; % Gravitational Constant [km^3/kgs^2]
    m_e = 5.9723*10^24; % Mass of Earth [kg]
    m_J = 1.898*10^27;  % Mass of Jupiter [kg]
    m_s = 1.989*10^30;  % Mass of Sun [kg]
    V_e = 30;           % Velocity of Earth [km/s]
    r_e = 6378;         % Radius of Earth [km]
    ri = r_e + 390;     % Radius of Earth parking orbit [km] 
    r_J = 71492;        % Radius of Jupiter [km]
    rf = r_J*10;        % Radius of Jupiter parking orbit [km]
    R_e = 149597870.7;  % Radius of Earth orbit around sun [km]
    R_J = R_e*5.2;      % Radius of Jupiter orbit around sun [km]
    u_s = G*m_s;        % mu Sun [km^3/s^2]
    u_e = G*m_e;        % mu Earth [km^3/s^2]
    u_J = G*m_J;        % mu Jupiter [km^3/s^2]
    AU = 1.496e+08;     %1 AU[km]
    
  %Calculations          
     v_r_o = dot(v_o,r_o)/norm(r_o);     %Initial radial velocity [km/s]
     coe = coe_from_sv(r_o,v_o,u_s);
     h = coe(1);
     e = coe(2);
     a = coe(7);
     T = (2*pi*a^(3/2))/u_s^.5;

        %Calculation of r,f,g,g',f' for range of theta
        r = (h^2/(u_s))*(1./(1+((h^2/((u_s)*norm(r_o)))-1)*cosd(theta)-(h*v_r_o/(u_s))*sind(theta)));
        f = 1-((u_s)*r/h^2).*(1-cosd(theta));
        g = (r*norm(r_o).*sind(theta))./h;
        g_p = 1-(((u_s)*norm(r_o)*(1-cosd(theta)))./(h^2));
        f_p = u_s/h*(v_r_o/h*(1-cosd(theta))-sind(theta)/r);
        
        %Resulting Radius and Velocity vectors for each change in theta
    for i=1:length(theta-1)
        for j=1:3
            R(i,j) = f(i)*r_o(j)+g(i)*v_o(j);
            V(i,j) = f_p(i)*r_o(j)+g_p(i)*v_o(j);    
        end
    end
    
    %TAKE THE THETA VALUE ALONG WITH RADIUS VALUE AND SOLVE FOR
    %CORRESPONDING TIMES SINCE PERIAPSIS FOR EACH THETA VALUE
    
    for i = 1:length(theta-1)
        thetaR(i) = theta(i)*pi/180;
        %if e<1
            Me(i) = 2*atan(sqrt((1-e)/(1+e))*tan(thetaR(i)/2))-e*sqrt(1-e^2)*sin(thetaR(i))/(1+e*cos(thetaR(i)));
            if Me(i)<0
                Me(i)=2*pi-abs(Me(i));
            end
            E(i) = kepler_E(e,Me(i));
            t(i) = (Me(i)*T/(2*pi))/86400;
        %else 
%             Mh(i) = e*sqrt(e^2-1)*sin(thetaR(i))/(1+e*cos(thetaR(i)))...
%                 -log(((sqrt(e+1))+(sqrt(e-1)*tan(thetaR(i)/2)))/((sqrt(e+1))-(sqrt(e-1)*tan(thetaR(i)/2))));
%             if Mh(i)<0
%                 Mh(i)=2*pi-abs(Mh(i));
%             end
%             F(i) = kepler_H(e,Mh(i));
%             t(i) = (Mh(i)*h^3/((u_s^2)*(e^2-1)^3/2))/86400;
            
    end
end
