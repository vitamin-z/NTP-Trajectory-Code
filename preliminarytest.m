%% Jesse Owens Nuclear Thermal Propulsion Spacecraft Trajectories
% Primary Author: Zach Zoloty
% Preliminary n-body Sim and Trajectory code
% Possible Departure Dates assuming a 2033 window:
% 4/6/33
% 4/28/33
% 4/20/33
% 1/26/33

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Currently not in running condition, %
% but this will be our main file.     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add a switch to run both the high energy trajectory [1] and Hohmann
% transfer [2]

clc;
clear;
%%

% Ephemerides
% Current ephemerides are from 4/20/33 to 4/20/35 at 6 hour increments
% Ephemerides for the celestial bodies listed below:
% Get them from 4/20/33 until about 3 years later. That should be a
% sufficient amount of time to have regardless of the wait time at Mars.
% As for the time step between data, we want it as small as possible
% obviously, but too small makes the files way too large and too large
% will make the spacecraft jump around through space in a very disjointed
% manner.

% Outsourcing direct ephemeride database integration to Kevin?
% Columns 2-4 are X,Y,Z and columns 5-7 are U,V,W
% We need X, Y, and Z to plot so those are the state arrays
Earth = table2array(readtable('Earth42033.txt'));
state_e(:,1) = Earth(:,1);
state_e(:,2) = Earth(:,2);
state_e(:,3) = Earth(:,3);
state_e(:,4) = Earth(:,4);
state_e(:,5) = Earth(:,5);
state_e(:,6) = Earth(:,6);

Mars = table2array(readtable('Mars42033.txt'));
state_m(:,1) = Mars(:,1);
state_m(:,2) = Mars(:,2);
state_m(:,3) = Mars(:,3);
state_m(:,4) = Mars(:,4);
state_m(:,5) = Mars(:,5);
state_m(:,6) = Mars(:,6);

Jupiter = table2array(readtable('Jupiter42033.txt'));
state_j(:,1) = Jupiter(:,1);
state_j(:,2) = Jupiter(:,2);
state_j(:,3) = Jupiter(:,3);
state_j(:,4) = Jupiter(:,4);
state_j(:,5) = Jupiter(:,5);
state_j(:,6) = Jupiter(:,6);

Phobos = table2array(readtable('Phobos42033.txt'));
state_p(:,1) = Phobos(:,1);
state_p(:,2) = Phobos(:,2);
state_p(:,3) = Phobos(:,3);
state_p(:,4) = Phobos(:,4);
state_p(:,5) = Phobos(:,5);
state_p(:,6) = Phobos(:,6);

Moon = table2array(readtable('Moon42033.txt'));
state_l(:,1) = Moon(:,1);
state_l(:,2) = Moon(:,2);
state_l(:,3) = Moon(:,3);
state_l(:,4) = Moon(:,4);
state_l(:,5) = Moon(:,5);
state_l(:,6) = Moon(:,6);

% Deimos = table2array(readtable('Deimos42033.txt'));
% state_d(:,1) = Deimos(:,1);
% state_d(:,2) = Deimos(:,2);
% state_d(:,3) = Deimos(:,3);
% state_d(:,4) = Deimos(:,4);
% state_d(:,5) = Deimos(:,5);
% state_d(:,6) = Deimos(:,6);

%%

% Define Constants
G      = 6.67408*10^-20;    % Gravitational Constant [km^3/kgs^2]
g_o    = .00981;            % Earth's gravitational paramater [km/s^2]

m_e    = 5.9723*10^24;      % Mass of Earth [kg]
m_m    = 641.9*10^21;       % Mass of Mars [kg]
m_p    = 1.0659*10^16;      % Mass of Phobos [kg]
m_j    = 1.898*10^27;       % Mass of Jupiter [kg]
m_s    = 1.989*10^30;       % Mass of Sun [kg]

r_e    = 6378;              % Radius of Earth [km]
ri     = r_e + 500;         % Radius of initial Earth parking orbit(500 km) 
r_m    = 3396;              % Radius of Mars [km]
rf     = r_m*10;            % Radius of final Mars parking orbit [km]
r_j    = 71492;             % Radius of Jupiter [km]
r_p    = 11;                % Radius of Phobos [km] (mean radius because of its irregular shape)

u_s    = G*m_s;             % mu Sun [km^3/s^2]
u_e    = G*m_e;             % mu Earth [km^3/s^2]
u_m    = G*m_m;             % mu Mars [km^3/s^2]
u_j    = G*m_j;             % mu Jupiter [km^3/s^2]
u_p    = G*m_p;             % mu Phobos [km^3/s^2]

d2r    = pi/180;            % degrees to radians
r2d    = 180/pi;            % radians to degrees

r_p_ES = 147166462;         % Periapsis of Earth to Sun [km]
r_a_ES = 152171522;         % Apoapsis of Earth to Sun [km]
a_ES   = (r_p_ES+r_a_ES)/2; % Semimajor axis Earth to Sun [km]
e_ES   = 0.0167;            % Eccentricity of earth around Sun
b_ES   = a_ES*(1-e_ES)^.5;  % Semiminor axis Earth to Sun [km]

a_JS   = 778000000;         % Semimajor axis Jupiter to Sun [km]
e_JS   = 0.048;             % Eccentriity of Jupiter around Sun
b_JS   = a_JS*(1-e_JS)^.5;  % Semiminor axis Jupiter to Sun [km]
r_p_JS = a_JS*(1-e_JS);     % Periapsis of Jupiter to Sun [km]

a_L    = 384400;            % Semimajor axis of Moon to Earth [km]
e_L    = 0.0549;            % Eccentricity of Moon orbit 
r_p_L  = a_L*(1-e_L);       % Perapsis of Moon to Earth [km]
b_L    = a_L*(1-e_L)^.5;    % Semiminor axis Moon to Earth [km]

a_MS   = 227800000;         % Semimajor axis of Mars to sun [km]
e_MS   = 0.0935;            % Eccentricity of Mars around sun
r_p_MS = a_MS*(1-e_MS);     % Periapsis of Mars to Sun [km]
b_MS   = a_MS*(1-e_MS)^.5;  % Semiminor axis Mars to Sun [km]

AU     = 1.496e+08;         % 1 Astronomical Unit [km]

%%

% Preliminary Trajectory
% Possible issue with length/timestep of ode45 result for spacecraft trajectory and
% timestep of the ephemeride data.
% Can sidestep this by utilizing the same acceleration/velocity/distance
% method as messed with in four.m (no need for ode45 funky timesteps that
% way).
% Unfortunately I'm unsure if these equations for distance and velocity end
% up being correct, but they should suffice for now to at least get the
% program up and running. After that we can screw around with ODE45's
% inability to set the time step.

                     %%%%%%%%%% Non-ODE45 method %%%%%%%%%%
% Jesse Owens initial conditions
% x0 = 
% y0 = 
% z0 = 
% vx0 = 
% vy0 = 
% vz0 = 
% 
% endtime = %total time covered by imported ephemeride data
% dt = %timestep used for ephemeride data
% i = endtime*dt;
% 
% Xjo = zeros(1,i);
% Yjo = zeros(1,i);
% Zjo = zeros(1,i);
% VXjo = zeros(1,i);
% VYjo = zeros(1,i);
% VZjo = zeros(1,i);
% 
% Xjo(1) = x0;
% Yjo(1) = y0;
% Zjo(1) = z0;
% VXjo(1) = vx0;
% VYjo(1) = vy0;
% VZjo(1) = vz0;
% 
% for j = 1:(i-1)
%    dxe = state_e(j,1) - Xjo(j); %X distance from the spacecraft to Earth
%    dxm = state_m(j,1) - Xjo(j); %X distance from the spacecraft to Mars
%    dxj = state_j(j,1) - Xjo(j); %X distance from the spacecraft to Jupiter
%    dxp = state_p(j,1) - Xjo(j); %X distance from the spacecraft to Phobos
%    dxs = 0 - Xjo(j);            %X distance from the spacecraft to the Sun
% 
%    
%    dye = state_e(j,2) - Yjo(j); %Y distance from the spacecraft to Earth
%    dym = state_m(j,2) - Yjo(j); %Y distance from the spacecraft to Mars
%    dyj = state_j(j,2) - Yjo(j); %Y distance from the spacecraft to Jupiter
%    dyp = state_p(j,2) - Yjo(j); %Y distance from the spacecraft to Phobos
%    dys = 0 - Yjo(j);            %Y distance from the spacecraft to the Sun
%    
%    dze = state_e(j,3) - Zjo(j); %Z distance from the spacecraft to Earth
%    dzm = state_m(j,3) - Zjo(j); %Z distance from the spacecraft to Mars
%    dzj = state_j(j,3) - Zjo(j); %Z distance from the spacecraft to Jupiter
%    dzp = state_p(j,3) - Zjo(j); %Z distance from the spacecraft to Phobos
%    dzs = 0 - Zjo(j);            %Z distance from the spacecraft to the Sun
%    
%    % Position vectors relative to the spacecraft
%    dre = [dxe dye dze];
%    drm = [dxm dym dzm];
%    drj = [dxj dyj dzj];
%    drp = [dxp dyp dzp];
%    drs = [dxs dys dzs];
%    
%    % Magnitudes
%    mre = norm(dre);
%    mrm = norm(drm);
%    mrj = norm(drj);
%    mrp = norm(drp);
%    mrs = norm(drs);
%    
%    % Acceleration
%    % going to need to mess with this equation as we go; currently the signs
%    % might not be right
%    ajo = u_s*drs/(mrs^3) + u_e*dre/(mre^3) + u_m*drm/(mrm^3) + u_j*drj/(mrj^3) + u_p*drp/(mrp^3);
%    
%    % Calculation and assignment of new spacecraft states
%    Xjo(j+1) = Xjo(j) + VXjo(j)*dt + 0.5*ajo(1)*dt^2;
%    Yjo(j+1) = Yjo(j) + VYjo(j)*dt + 0.5*ajo(2)*dt^2;
%    Zjo(j+1) = Zjo(j) + VZjo(j)*dt + 0.5*ajo(3)*dt^2;
%    VXjo(j+1) = VXjo(j) + ajo(1)*dt;
%    VYjo(j+1) = VYjo(j) + ajo(2)*dt;
%    VZjo(j+1) = VZjo(j) + ajo(3)*dt;
% end   
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R_e    = 149597870.7;       % Radius of Earth orbit around sun [km]                     
                     
V_e = sqrt(u_s/R_e);            % Velocity around sun [km/s]
V_excessANT = sqrt(23);         % Excess velocity (C3~23) available from Antares [km/s]
V_launch = V_e+V_excessANT;     % Initial launch velocity with Antares rocket [km/s]  
r_e_SOI = 925000;               % Radius of Earth SOI
r_launch = R_e+r_e_SOI;         % Launch Radius at burnout of Antares

r_J_SOI = 48200000;             % Radius of Jupiter SOI
r_p_cap = 10*r_j;               % Perigee radius of capture orbit
e_cap = 0.866;                  % Eccentricity of capture orbit
V_cap = sqrt(u_j*(1+e_cap)/r_p_cap);     % Capture velocity for capture orbit at SOI

R_J    = R_e*5.2;           % Radius of Jupiter orbit around sun [km]
r_f_des = R_J - r_J_SOI;        % Radius magnitude at capture

                     %%%%%%%%%% RUNS FOR 2 YEARS %%%%%%%%%%
% currently arbitrary`
t = [0:10:63093600];                % Duplicate time index

thrust = linspace(0,0,length(t));   % Populate Thrust vector

MaxMag = 111;                       % Maximum Thrust Magnitude [kN]

% Thrust model used by previous group
thrust_1min = MaxMag*[.4 .8 .98 .98 .8 .4];               % Thrust vector modeled as cosine wave
thrust_1_5min = MaxMag*[.4 .8 .98 1 1 1 .98 .8 .4];               % Thrust vector modeled as cosine wave
thrust_100s = MaxMag*[.4 .8 .98 1 1 1 1 .98 .8 .4];               % Thrust vector modeled as cosine wave
thrust_2min = MaxMag*[.2 .6 .85 .9 .98 .98 .98 .98 .9 .86 .6 .2];               % Thrust vector modeled as cosine wave
thrust_25min = MaxMag*[.2 .6 .85 .9 .98 .98 .99 1 .99 .98 .98 .9 .86 .6 .2];               % Thrust vector modeled as cosine wave
thrust_3min = MaxMag*[.1 .5 .9 .98 .98 .99 .99 1 1 1 1 .99 .99 .98 .98 .9 .5 .15];               % Thrust vector modeled as cosine wave
thrust_35min = MaxMag*[.1 .5 .9 .98 .98 .99 .99 1 1 1 1 1 1 1 .99 .99 .98 .98 .9 .5 .15];               % Thrust vector modeled as cosine wave
thrust_4min = MaxMag*[.1 .5 .7 .75 .85 .9 .98 .98 .99 1 1 1 1 1 1 .99 .98 .98 .9 .85 .75 .7 .5 .1];               % Thrust vector modeled as cosine wave
thrust_5min = MaxMag*[.1 .5 .7 .75 .85 .9 .98 .98 .99 1 1 1 1 1 1 1 1 1 1 1 1 .99 .98 .98 .9 .85 .75 .7 .5 .1];               % Thrust vector modeled as cosine wave
thrust_55min = MaxMag*[.1 .5 .7 .75 .85 .9 .98 .98 .99 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 .99 .98 .98 .9 .85 .75 .7 .5 .1];               % Thrust vector modeled as cosine wave
thrust_6min = MaxMag*[.1 .5 .7 .75 .85 .9 .98 .98 .99 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 .99 .98 .98 .9 .85 .75 .7 .5 .1];               % Thrust vector modeled as cosine wave
thrust_8min = MaxMag*[.1 .5 .7 .75 .85 .9 .98 .98 .99 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 .99 .98 .98 .9 .85 .75 .7 .5 .1];               % Thrust vector modeled as cosine wave
thrust_85min = MaxMag*[.1 .5 .7 .75 .85 .9 .98 .98 .99 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 .99 .98 .98 .9 .85 .75 .7 .5 .1];               % Thrust vector modeled as cosine wave
thrust_9min = MaxMag*[.1 .5 .7 .75 .85 .9 .98 .98 .99 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 .99 .98 .98 .9 .85 .75 .7 .5 .1];               % Thrust vector modeled as cosine wave
thrust_10min = MaxMag*[.1 .5 .7 .75 .85 .9 .98 .98 .99 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 .99 .98 .98 .9 .85 .75 .7 .5 .1];               % Thrust vector modeled as cosine wave
thrust_12min = MaxMag*[.1 .5 .7 .75 .85 .9 .98 .98 .99 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 .99 .98 .98 .9 .85 .75 .7 .5 .1];               % Thrust vector modeled as cosine wave
thrust_15min = MaxMag*[.1 .5 .7 .75 .85 .9 .98 .98 .99 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 .99 .98 .98 .9 .85 .75 .7 .5 .1];               % Thrust vector modeled as cosine wave
thrust_18min = MaxMag*[.1 .5 .7 .75 .85 .9 .98 .98 .99 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 .99 .98 .98 .9 .85 .75 .7 .5 .1];               % Thrust vector modeled as cosine wave
thrust_20min = MaxMag*[.1 .5 .7 .75 .85 .9 .98 .98 .99 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 .99 .98 .98 .9 .85 .75 .7 .5 .1];               % Thrust vector modeled as cosine wave

%%

j = 1;
q = 1;
s = 1;
d = 1;

for i=1:length(t)
    if i > 0 && i< 52
    THRUST(i) = thrust_85min(j);      % s burn @ departure
        j = j+1;
%     elseif i > 8640 && i< 8647    
%         THRUST(i) = thrust_1min(q);      % burn @ 1 day
%         q = q+1;    
    elseif i > 2937600 && i< 2937611   
        THRUST(i) = -thrust_100s(s);     % s burn @ 340 days
        s = s+1;
%     elseif i > 2937612 && i< 2937631        
%         THRUST(i) = -thrust_3min(d);         % s burn @ 331 days
%         d = d+1;
    else 
        THRUST(i) = 0;       
    end
end

%%

% T = thrust [N]
% E = error/noise

% constants
Isp = 850;
g0 = 9.81*10^-3;
%mu = 1.3275*10^11;

% reading the OED variable from Y
y1 = Y(1);  % x
y2 = Y(2);  % y
y3 = Y(3);  % z
y4 = Y(4);  % xdot
y5 = Y(5);  % ydot
y6 = Y(6);  % zdot
m  = Y(7);  % m

% Computing velocity, position 
V = sqrt(y4^2+y5^2+y6^2); % velocity
r = sqrt(y1^2+y2^2+y3^2); % position

r_earth(1) = y1 - state_e(k,1);
r_earth(2) = y2 - state_e(k,2);
r_earth(3) = y3 - state_e(k,3);
r_jupiter(1) = y1 - state_j(k,1);
r_jupiter(2) = y2 - state_j(k,1);
r_jupiter(3) = y3 - state_j(k,1);
r_mars(1) = y1 - state_m(k,1);
r_mars(2) = y2 - state_m(k,1);
r_mars(3) = y3 - state_m(k,1);

%
% Computing Ydot (ODE implementation)
% DOUBLE CHECK THAT THIS IS CORRECT
Ydot(k+1,:) = [y4;
        y5;
        y6;
       -u_s*y1/r^3 -u_e*y1/r^3 -u_j*y1/r^3 -u_m*y1/r^3 + T/m*y4/V + E;   
       -u_s*y2/r^3 -u_e*y2/r^3 -u_j*y2/r^3 -u_m*y2/r^3 + T/m*y5/V + E;
       -u_s*y3/r^3 -u_e*y3/r^3 -u_j*y3/r^3 -u_m*y3/r^3 + T/m*y6/V + E;
       -abs(T)/(Isp*g0)
       ];
   
   %%

simin = [t',THRUST'];

SimOut = sim('PrevThrustTestingNew');

total_days = Time(end)/86400;                % Total time ran for [days]
fprintf('        %.f DAYS\n',total_days)

v_mag = sqrt(vx.^2+vy.^2+vz.^2);    % Velocity magnitude string [km/s]
r_mag = sqrt(x.^2+y.^2+z.^2);    % Radius magnitude string [km]

indx = find(r_mag>=r_f_des);

first = indx(1);

v_final = sqrt(vx(first)^2+vy(first)^2+vz(first)^2);    % Final velocity magnitude [km/s]
r_final = sqrt(x(first)^2+y(first)^2+z(first)^2);   % Final radius magnitude [km]
t_final = Time(first)/86400;                % Time on transit [days]
fprintf('Desired radius at rendezvous: %.4f [AU]\n',r_f_des/AU)
fprintf('Desired velocity at rendezvous: %.2f [km/s]\n',V_cap)


mp_o = m(1)-700;              % Initial mass of fuel [kg]
mp_spent = 9700-m(first);     % Mass of fuel spent [kg]
mp_left = 9000-mp_spent;    % Mass of fuel left [kg]
fprintf('Initial fuel mass: %.2f [kg]\n\n',mp_o)

fprintf('Velocity at rendezvous: %.2f [km/s]\n',v_final)
fprintf('Radius at rendezvous: %.4f [AU]\n',r_final)
fprintf('Propellant mass left: %.2f [kg]\n',mp_left)
fprintf('Time on transit: %.2f [days]\n',t_final)

coe = coe_from_sv([x(end) y(end) z(end)],[vx(end) vy(end) vz(end)],u_s);
h = coe(1); e = coe(2); RA = coe(3); incl = coe(4); w = coe(5); TA = coe(6); a = coe(7);
A = [e a h];

dV_p_cap = sqrt(V_cap^2+(2*u_j/r_p_cap))-sqrt(u_j*(1+e_cap)/r_p_cap);

e_cap_actual = ((v_final^2)*r_p_cap/u_j)-1;      % Actual eccentricity of capture orbit
dV_p_cap_actual = sqrt(v_final^2+(2*u_j/r_p_cap))-sqrt(u_j*(1+e_cap_actual)/r_p_cap);   % Actual dV required to ellipticize [km/s]

fprintf('New eccentricity of Jovian elliptical orbit based off of hyperbolic excess velocity: %.4f\n',e_cap_actual)
fprintf('Actual dV required at perigee: %.3f [km/s]\n',dV_p_cap_actual)

MR = exp(-dV_p_cap_actual/(850*g_o));
m_req = MR*mp_left;
m_left_cap = mp_left-m_req; 

fprintf('Approximate fuel remaining after capture: %.2f [kg]\n',m_left_cap)

%%

% 3D Simulation
% Outsource visualization magic to Kevin?
% The variables here are currently from four.m
% You can find four.m on GitHub  to see how this plot function works.
% Unsure whether this or the previously way we were plotting is more
% efficient. We'll have to experiment (or replace it with Kevin magic).

% posit=plot3(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,x1,y1,z1,'ob',x2,y2,z2,'or',x3,y3,z3,'oy',x4,y4,z4,'ok');
posit = plot3(state_e(:,1),state_e(:,2),state_e(:,3),state_m(:,1),state_m(:,2),state_m(:,3),state_j(:,1),state_j(:,2),state_j(:,3),state_p(:,1),state_p(:,2),state_p(:,3),state_l(:,1),state_l(:,2),state_l(:,3),x,y,z);%state_d(:,1),state_d(:,2),state_d(:,3),x,y,z);
title('Positions as a function of t');
xlabel('x');
ylabel('y');
zlabel('z');
daspect([1 1 1]);
set(gca,'PlotBoxAspectRatio',[1 1 1]);
grid on;

view(-37.5,30);

% for rd=1:i/300:i
%     
%     rd=floor(rd);
%     
%     % reset all values for position and velocity, then redraw
%     
%     set(posit(1),'XData',X1(1:rd));
%     set(posit(1),'YData',Y1(1:rd));
%     set(posit(1),'ZData',Z1(1:rd));
%     set(posit(2),'XData',X2(1:rd));
%     set(posit(2),'YData',Y2(1:rd));
%     set(posit(2),'ZData',Z2(1:rd));
%     set(posit(3),'XData',X3(1:rd));
%     set(posit(3),'YData',Y3(1:rd));
%     set(posit(3),'ZData',Z3(1:rd));
%     set(posit(4),'XData',X4(1:rd));
%     set(posit(4),'YData',Y4(1:rd));
%     set(posit(4),'ZData',Z4(1:rd));
%     
%     set(posit(5),'XData',X1(rd));
%     set(posit(5),'YData',Y1(rd));
%     set(posit(5),'ZData',Z1(rd));
%     set(posit(6),'XData',X2(rd));
%     set(posit(6),'YData',Y2(rd));
%     set(posit(6),'ZData',Z2(rd));
%     set(posit(7),'XData',X3(rd));
%     set(posit(7),'YData',Y3(rd));
%     set(posit(7),'ZData',Z3(rd));
%     set(posit(8),'XData',X4(rd));
%     set(posit(8),'YData',Y4(rd));
%     set(posit(8),'ZData',Z4(rd));
%     
%     drawnow
% end