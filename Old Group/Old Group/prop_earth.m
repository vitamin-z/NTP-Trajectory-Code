function [R,V] = prop_earth(r_o, v_o,theta)

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
     h = norm(cross(r_o,v_o));           %Angular Momentum [kg*km^2/s]
     v_r_o = dot(v_o,r_o)/norm(r_o);     %Initial radial velocity [km/s]

        %Calculation of r,f,g,g',f' for range of theta
        r = (h^2/u_e)*(1./(1+((h^2/(u_e*norm(r_o)))-1)*cosd(theta)-(h*v_r_o/u_e)*sind(theta)));
        f = 1-(u_e*r/h^2).*(1-cosd(theta));
        g = (r*norm(r_o).*sind(theta))./h;
        g_p = 1-((u_e*norm(r_o)*(1-cosd(theta)))./(h^2));
        f_p = (u_e.*(1-cosd(theta))./(h.*sind(theta))).*((u_e.*(1-cosd(theta))./h^.2)-(1/norm(r_o))-(1./r));

        %Resulting Radius and Velocity vectors for each change in theta
    for i=1:length(theta)
        for j=1:3
            R(i,j) = f(i)*r_o(j)+g(i)*v_o(j);
            V(i,j) = f_p(i)*r_o(j)+g_p(i)*v_o(j);
        end
    end
end

