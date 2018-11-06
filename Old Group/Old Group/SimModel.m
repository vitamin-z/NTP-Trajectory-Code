%% Jesse Owens Trajctories
% This script will explore some basic trajectories to Jupiter including a
% Hohmann Transfer and other simple burn paths.
%
% Date Updated: 16 October 2017
% Written by:
%
% * Zach Allen
% * Skylar Dierker
% * Craig Heim
% * Ryan Howard
% * Taylor Huneycutt
% * Stephen Ioas
 
%% Initial Conditions and Constants
clear all; close all; clc
 
G      = 6.67408*10^-20;    % Gravitational Constant [km^3/kgs^2]
g_o    = .00981;            % Earth's gravitational paramater [km/s^2
m_e    = 5.9723*10^24;      % Mass of Earth [kg]
m_J    = 1.898*10^27;       % Mass of Jupiter [kg]
m_s    = 1.989*10^30;       % Mass of Sun [kg]
u_s    = G*m_s;             % mu Sun [km^3/s^2]
R_e    = 149597870.7;       % Radius of Earth orbit around sun [km]
V_e    = sqrt(u_s/R_e);     % Velocity of Earth [km/s]
r_e    = 6378;              % Radius of Earth [km]
ri     = r_e + 500;         % Radius of Earth parking orbit(500 km) 
r_J    = 71492;             % Radius of Jupiter [km]
rf     = r_J*10;            % Radius of Jupiter parking orbit [km]
r_o    = 42146;
R_J    = R_e*5.2;           % Radius of Jupiter orbit around sun [km]
u_e    = G*m_e;             % mu Erth [km^3/s^2]
u_J    = G*m_J;             % mu Jupiter [km^3/s^2]
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
r_focus_earth = a_ES-r_p_ES;      % Location of focus in Earth orbit [km]
r_focus_jupiter = a_JS-r_p_JS;    % Location of focus in Jupiter orbit [km]

%%
V_e = sqrt(u_s/R_e);            % Velocity around sun [km/s]
V_excessANT = sqrt(23);         % Excess velocity (C3~23) available from Antares [km/s]
V_launch = V_e+V_excessANT;     % Initial launch velocity with Antares rocket [km/s]  
r_e_SOI = 925000;               % Radius of Earth SOI
r_launch = R_e+r_e_SOI;         % Launch Radius at burnout of Antares

r_J_SOI = 48200000;             % Radius of Jupiter SOI
r_p_cap = 10*r_J;               % Perigee radius of capture orbit
e_cap = 0.866;                  % Eccentricity of capture orbit
V_cap = sqrt(u_J*(1+e_cap)/r_p_cap);     % Capture velocity for capture orbit at SOI

r_f_des = R_J - r_J_SOI;        % Radius magnitude at capture

                     %%%%%%%%%% RUNS FOR 341 DAYS %%%%%%%%%%
t = [0:10:29462400];                % Duplicate time index

THRUST = linspace(0,0,length(t));   % Populate Thrust vector

MaxMag = 111;                       % Maximum Thrust Magnitude [kN]

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

simin = [t',THRUST'];

SimOut = sim('Prev_ThrustTesting');

total_days = Time(end)/86400;                % Total time ran for [days]
fprintf('        %.f DAYS\n',total_days)

v_mag = sqrt(vx.^2+vy.^2);    % Velocity magnitude string [km/s]
r_mag = sqrt(x.^2+y.^2);    % Radius magnitude string [km]

indx = find(r_mag>=r_f_des);

first = indx(1);

v_final = sqrt(vx(first)^2+vy(first)^2);    % Final velocity magnitude [km/s]
r_final = sqrt(x(first)^2+y(first)^2)/AU;   % Final radius magnitude [AU]
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

figure
hold on
earth = ellipse(a_ES/AU,b_ES/AU,0,r_focus_earth/AU,0,'b'); % Earth Orbit
jupiter = ellipse(a_JS/AU,b_JS/AU,0,r_focus_jupiter/AU,0,'g'); % Jupiter Orbit

plot(x/AU,y/AU,'r')
    
sun = circle(0,0,30000000/AU);  % Sun Location **NOT TO SCALE**

legend('Earth','Jupiter','Analytical')
axis square
title('Jesse Owens Transit to Jupiter')
xlabel('AU')
ylabel('AU')

coe = coe_from_sv([x(end) y(end) z(end)],[vx(end) vy(end) vz(end)],u_s);
h = coe(1); e = coe(2); RA = coe(3); incl = coe(4); w = coe(5); TA = coe(6); a = coe(7);
A = [e a h];

dV_p_cap = sqrt(V_cap^2+(2*u_J/r_p_cap))-sqrt(u_J*(1+e_cap)/r_p_cap);

e_cap_actual = ((v_final^2)*r_p_cap/u_J)-1;      % Actual eccentricity of capture orbit
dV_p_cap_actual = sqrt(v_final^2+(2*u_J/r_p_cap))-sqrt(u_J*(1+e_cap_actual)/r_p_cap);   % Actual dV required to ellipticize [km/s]

fprintf('New eccentricity of Jovian elliptical orbit based off of hyperbolic excess velocity: %.4f\n',e_cap_actual)
fprintf('Actual dV required at perigee: %.3f [km/s]\n',dV_p_cap_actual)

MR = exp(-dV_p_cap_actual/(850*g_o));
m_req = MR*mp_left;
m_left_cap = mp_left-m_req; 

fprintf('Approximate fuel remaining after capture: %.2f [kg]\n',m_left_cap)
