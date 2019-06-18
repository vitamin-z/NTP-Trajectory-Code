%% Jesse Owens Nuclear Thermal Propulsion Spacecraft Trajectories
% Primary Author: Zach Zoloty
% Preliminary n-body Sim and Trajectory code
% Possible Departure Dates assuming a 2033 window:
% 4/6/33
% 4/28/33
% 4/20/33
% 1/26/33

%initial conditions currently for 4/20/33 launch date

% Simulink previously ran for 34560000 seconds?
% Accelerator or Rapid Accelerator for Sim?

% Add a switch to run both the high energy trajectory [1] and Hohmann
% transfer [2]

% Change it from 1 second to something else
% Change the thrust vectors to match the time interval?

% Have multiple time intervals? that way we can move through space faster
% when we aren't thrusting

%use sv to coe to check final JO orbit

clc;
clear;
%%

% Ephemerides
% Get them from 4/20/33 until about 3 years later. That should be a
% sufficient amount of time to have regardless of the wait time at Mars.
% As for the time step between data, we want it as small as possible
% obviously, but too small makes the files way too large and too large
% will make the spacecraft jump around through space in a very disjointed
% manner.

% Columns 2-4 are X,Y,Z and columns 5-7 are U,V,W
% We need X, Y, and Z to plot so those are the state arrays

addpath('D:\School\Capstone\NTP Code 2018-2019\SPICE\mice\src\mice');
addpath('D:\School\Capstone\NTP Code 2018-2019\SPICE\mice\lib');

cspice_furnsh('de432s.bsp');
cspice_furnsh('mar097.bsp');
cspice_furnsh('naif0010.tls');

% Summary for: de432s.bsp
%  
% Bodies: MERCURY BARYCENTER (1) w.r.t. SOLAR SYSTEM BARYCENTER (0)
%         VENUS BARYCENTER (2) w.r.t. SOLAR SYSTEM BARYCENTER (0)
%         EARTH BARYCENTER (3) w.r.t. SOLAR SYSTEM BARYCENTER (0)
%         MARS BARYCENTER (4) w.r.t. SOLAR SYSTEM BARYCENTER (0)
%         JUPITER BARYCENTER (5) w.r.t. SOLAR SYSTEM BARYCENTER (0)
%         SATURN BARYCENTER (6) w.r.t. SOLAR SYSTEM BARYCENTER (0)
%         URANUS BARYCENTER (7) w.r.t. SOLAR SYSTEM BARYCENTER (0)
%         NEPTUNE BARYCENTER (8) w.r.t. SOLAR SYSTEM BARYCENTER (0)
%         PLUTO BARYCENTER (9) w.r.t. SOLAR SYSTEM BARYCENTER (0)
%         SUN (10) w.r.t. SOLAR SYSTEM BARYCENTER (0)
%         MERCURY (199) w.r.t. MERCURY BARYCENTER (1)
%         VENUS (299) w.r.t. VENUS BARYCENTER (2)
%         MOON (301) w.r.t. EARTH BARYCENTER (3)
%         EARTH (399) w.r.t. EARTH BARYCENTER (3)
%         Start of Interval (ET)              End of Interval (ET)
%         -----------------------------       -----------------------------
%         1949 DEC 14 00:00:00.000            2050 JAN 02 00:00:00.000
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Summary for: mar097.bsp
%  
% Bodies: EARTH BARYCENTER (3) w.r.t. SOLAR SYSTEM BARYCENTER (0)
%         MARS BARYCENTER (4) w.r.t. SOLAR SYSTEM BARYCENTER (0)
%         SUN (10) w.r.t. SOLAR SYSTEM BARYCENTER (0)
%         EARTH (399) w.r.t. EARTH BARYCENTER (3)
%         PHOBOS (401) w.r.t. MARS BARYCENTER (4)
%         DEIMOS (402) w.r.t. MARS BARYCENTER (4)
%         MARS (499) w.r.t. MARS BARYCENTER (4)
%         Start of Interval (ET)              End of Interval (ET)
%         -----------------------------       -----------------------------
%         1900 JAN 04 00:00:41.184            2100 JAN 01 00:01:07.183
%

tprime = 10540800; %four months %use 2 months for GTO initial conditions
tinterval = 5*60;%25 second interval%60*60*6; %time interval in seconds
tinstance = tprime/tinterval;
start = cspice_str2et('April 20 2033');
et = (0:tinstance)*tinterval + start;

[state_e,ltime_e] = cspice_spkezr('EARTH BARYCENTER',et,'J2000','LT+S','SOLAR SYSTEM BARYCENTER');
[state_m,ltime_m] = cspice_spkezr('MARS BARYCENTER',et,'J2000','LT+S','SOLAR SYSTEM BARYCENTER');
[state_j,ltime_j] = cspice_spkezr('JUPITER BARYCENTER',et,'J2000','LT+S','SOLAR SYSTEM BARYCENTER');
[state_pa,ltime_pa] = cspice_spkezr('PHOBOS',et,'J2000','LT+S','MARS BARYCENTER');
[state_da,ltime_da] = cspice_spkezr('DEIMOS',et,'J2000','LT+S','MARS BARYCENTER');
[state_la,ltime_la] = cspice_spkezr('MOON',et,'J2000','LT+S','EARTH BARYCENTER');
state_p = state_pa + state_m;
state_d = state_da + state_m;
state_l = state_la + state_e;

state_e = state_e';
state_m = state_m';
state_j = state_j';
state_p = state_p';
state_d = state_d';
state_l = state_l';

% assuming approximate 4 month trip, use 8 month arbitrary "Tom Petty Wait
% Time"
% tprime2 = 2*10540800; %8months
% tinterval2 = 5*60;
% tinstance2 = tprime2/tinterval2;
% start2 = cspice_str2et('August 20 2033');
% et2 = (0:tinstance2)*tinterval2 + start2;
% 
% [between_e,lbetween_e] = cspice_spkezr('EARTH BARYCENTER',et2,'J2000','LT+S','SOLAR SYSTEM BARYCENTER');
% [between_m,lbetween_m] = cspice_spkezr('MARS BARYCENTER',et2,'J2000','LT+S','SOLAR SYSTEM BARYCENTER');
% [between_j,lbetween_j] = cspice_spkezr('JUPITER BARYCENTER',et2,'J2000','LT+S','SOLAR SYSTEM BARYCENTER');
% between_e = between_e';
% between_m = between_m';
% between_j = between_j';

    tprime2 = 3*10540800; %12months %for mission one
    tinterval2 = 5*60;
    tinstance2 = tprime2/tinterval2;
    start2 = cspice_str2et('June 20 2033');
    et2 = (0:tinstance2)*tinterval2 + start2;

    [between_e,lbetween_e] = cspice_spkezr('EARTH BARYCENTER',et2,'J2000','LT+S','SOLAR SYSTEM BARYCENTER');
    [between_m,lbetween_m] = cspice_spkezr('MARS BARYCENTER',et2,'J2000','LT+S','SOLAR SYSTEM BARYCENTER');
    [between_j,lbetween_j] = cspice_spkezr('JUPITER BARYCENTER',et2,'J2000','LT+S','SOLAR SYSTEM BARYCENTER');
    between_e = between_e';
    between_m = between_m';
    between_j = between_j';

tprime3 = 10540800;
tinterval3 = 5*60;
tinstance3 = tprime3/tinterval3;
start3 = cspice_str2et('April 20 2034');
et3 = (0:tinstance3)*tinterval3 + start3;

[return_e,lreturn_e] = cspice_spkezr('EARTH BARYCENTER',et3,'J2000','LT+S','SOLAR SYSTEM BARYCENTER');
[return_m,lreturn_m] = cspice_spkezr('MARS BARYCENTER',et3,'J2000','LT+S','SOLAR SYSTEM BARYCENTER');
[return_j,lreturn_j] = cspice_spkezr('JUPITER BARYCENTER',et3,'J2000','LT+S','SOLAR SYSTEM BARYCENTER');
return_e = return_e';
return_m = return_m';
return_j = return_j';

cspice_kclear();



ma = length(state_m(:,1));

%%

% Define Constants
G      = 6.67408*10^-20;    % Gravitational Constant [km^3/kgs^2]
g_o    = .00981;            % Earth's gravitational parameter [km/s^2]

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

R_e    = 149597870.7;       % Radius of Earth orbit around sun [km]                     
                     
V_e = sqrt(u_s/R_e);            % Velocity around sun [km/s]
V_excessANT = sqrt(23);         % Excess velocity (C3~23) available from Antares [km/s]
V_launch = V_e+V_excessANT;     % Initial launch velocity with Antares rocket [km/s]  
r_e_SOI = 925000;               % Radius of Earth SOI
r_launch = R_e+r_e_SOI;         % Launch Radius at burnout of Antares

% r_J_SOI = 48200000;             % Radius of Jupiter SOI
r_p_cap = 2*r_m;               % Perigee radius of capture orbit
e_cap = 1;%0.866;                  % Eccentricity of capture orbit %relatively arbitrary for now
V_cap = sqrt(u_m*(1+e_cap)/r_p_cap);     % Capture velocity for capture orbit at SOI
% 
% R_J    = R_e*5.2;           % Radius of Jupiter orbit around sun [km]

r_M_SOI = 577000; %km
R_M = 227.9e6; %km

r_f_des = sqrt(state_m(ma,1)^2+state_m(ma,2)^2+state_m(ma,3)^2) - r_M_SOI;        % Radius magnitude at capture

                     %%%%%%%%%% RUNS FOR 2 YEARS %%%%%%%%%%
% currently arbitrary
%t = [0:10:63093600];                % Duplicate time index
t = [0:tinterval:tprime];                % Duplicate time index

thrust = linspace(0,0,length(t));   % Populate Thrust vector

MaxMag = 111;                       % Maximum Thrust Magnitude [kN]

%from last year's code, approx 6183.30635kg fuel spent in 5200s of burns
%this gives an mdot of approximately 1.1891 kg/s
mdotavg = 1.1891; %kg/s

% Thrust model
thrust_1min = MaxMag*[.4 .8 .98 .98 .8 .4];               % Thrust vector modeled as cosine wave
thrust_1_5min = MaxMag*[.4 .8 .98 1 1 1 .98 .8 .4];               % Thrust vector modeled as cosine wave
thrust_100s = MaxMag*[.4 .8 .98 1 1 1 1 .98 .8 .4];               % Thrust vector modeled as cosine wave
thrust_2min = MaxMag*[.2 .6 .85 .9 .98 .98 .98 .98 .9 .86 .6 .2];               % Thrust vector modeled as cosine wave
thrust_25min = MaxMag*[.2 .6 .85 .9 .98 .98 .99 1 .99 .98 .98 .9 .86 .6 .2];               % Thrust vector modeled as cosine wave
thrust_3min = MaxMag*[.1 .5 .9 .98 .98 .99 .99 1 1 1 1 .99 .99 .98 .98 .9 .5 .15];               % Thrust vector modeled as cosine wave
thrust_4min = MaxMag*[.1 .5 .7 .75 .85 .9 .98 .98 .99 1 1 1 1 1 1 .99 .98 .98 .9 .85 .75 .7 .5 .1];               % Thrust vector modeled as cosine wave
thrust_5min = MaxMag*[.1 .5 .7 .75 .85 .9 .98 .98 .99 1 1 1 1 1 1 1 1 1 1 1 1 .99 .98 .98 .9 .85 .75 .7 .5 .1];               % Thrust vector modeled as cosine wave
thrust_6min = MaxMag*[.1 .5 .7 .75 .85 .9 .98 .98 .99 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 .99 .98 .98 .9 .85 .75 .7 .5 .1];               % Thrust vector modeled as cosine wave
thrust_8min = MaxMag*[.1 .5 .7 .75 .85 .9 .98 .98 .99 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 .99 .98 .98 .9 .85 .75 .7 .5 .1];               % Thrust vector modeled as cosine wave
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

%change this; they have it thrusting for the first X number of indices
%rather than seconds
% for i=1:length(t)
%     if i > 0 && i< 52
%     THRUST(i) = thrust_25min(j);      % s burn @ departure
%         j = j+1;
% %     elseif i > 8640 && i< 8647    
% %         THRUST(i) = thrust_1min(q);      % burn @ 1 day
% %         q = q+1;    
%     elseif i > 2937600 && i< 2937611   
%         THRUST(i) = -thrust_100s(s);     % s burn @ 340 days
%         s = s+1;
% %     elseif i > 2937612 && i< 2937631        
% %         THRUST(i) = -thrust_3min(d);         % s burn @ 331 days
% %         d = d+1;
%     else 
%         THRUST(i) = 0;       
%     end
% end
%if tinterval == 1

burntime = 0;
del_v_out = 0;

fun = @(x) 111./(9700-mdotavg*(x));

   for i=1:length(t)
      if i>=0 && i<= 5%1500 %25 minute burn
          THRUST(i) = MaxMag; %assumed instantaneous reactor
          del_v_out = del_v_out + integral(fun,0,300);
          burntime = burntime + 5*60; %5 minute time intervals for burns
          
      elseif i>=17000 && i<=17002 %this one is for mission3
          THRUST(i) = MaxMag;
          del_v_out = del_v_out + integral(fun,0,300);
          burntime = burntime + 5*60; %5 minute time intervals for burns
          
%       elseif i>=17568 && i<=17573
%           THRUST(i) = MaxMag;
%       elseif i>=172800 && i<=174300
%           THRUST(i) = MaxMag;
      else
          THRUST(i) = 0; %reactor cooldown period to remove neutron poisons
      end
   end
%end
%293/298 for a day
%22233
%initcondits = [14.3013 -34.5845 -21]

% run multiple simulink things?
% one for small time steps, then bigger time steps, etc
% would need to have initial conditions change as we go
% THRUST = zeros(length(t));
% THRUST(1) = MaxMag;
% THRUST(58) = MaxMag;
   %%
%init_condit = [state_e(1,1)+ri,state_e(1,2),state_e(1,3),0,34.5845,0,9700];
   
simin = [t',THRUST'];
sim_e1 = [t',state_e(:,1)];
sim_e2 = [t',state_e(:,2)];
sim_e3 = [t',state_e(:,3)];
sim_j1 = [t',state_j(:,1)];
sim_j2 = [t',state_j(:,2)];
sim_j3 = [t',state_j(:,3)];
sim_m1 = [t',state_m(:,1)];
sim_m2 = [t',state_m(:,2)];
sim_m3 = [t',state_m(:,3)];

%Put it initial conditions
% posssibly assume circular orbit at 500km around earth to get velocities
% and positions and initial mass (initial mass can be found below)
%sim_x0 = 

% y10 = -131169500;
% y20 = -68577000;
% y30 = -29701000;
% y40 = 14.3013;
% y50 = -31.4484;
% y60 = -10.3323;

SimOut = sim('PrevThrustTestingUpdate');
xout = x;
yout = y;
zout = z;
vxout = vx;
vyout = vy;
vzout = vz;
%SimOut = sim('Prev_ThrustTesting');

%return trip
t2 = [0:tinterval3:tprime3];                % Duplicate time index

   for i=1:length(t2)
      if i>=0 && i<= 5%1500 %25 minute burn
          THRUST2(i) = MaxMag; %assumed instantaneous reactor
      elseif i>=293 && i<=298 %this one is for mission3
          THRUST2(i) = MaxMag;
%       elseif i>=17568 && i<=17573
%           THRUST(i) = MaxMag;
%       elseif i>=172800 && i<=174300
%           THRUST(i) = MaxMag;
      else
          THRUST2(i) = 0; %reactor cooldown period to remove neutron poisons
      end
   end
   
simin = [t2',THRUST2'];
sim_e1 = [t2',return_e(:,1)];
sim_e2 = [t2',return_e(:,2)];
sim_e3 = [t2',return_e(:,3)];
sim_j1 = [t2',return_j(:,1)];
sim_j2 = [t2',return_j(:,2)];
sim_j3 = [t2',return_j(:,3)];
sim_m1 = [t2',return_m(:,1)];
sim_m2 = [t2',return_m(:,2)];
sim_m3 = [t2',return_m(:,3)];

SimOut = sim('PrevThrustTestingReturn');
xreturn = x;
yreturn = y;
zreturn = z;
vxreturn = vx;
vyreturn = vy;
vzreturn = vz;

outbound_days = tprime/86400;                % Total time ran for [days]
fprintf('Outbound Trajectory Length:        %.f Days\n',outbound_days)
% 
v_mag = sqrt(vxout.^2+vyout.^2+vzout.^2);    % Velocity magnitude string [km/s]
r_mag = sqrt(xout.^2+yout.^2+zout.^2);    % Radius magnitude string [km]

% indx = find(r_mag>=r_f_des);
% 
% first = indx(1);
% 
% v_final = sqrt(vxout(first)^2+vyout(first)^2+vzout(first)^2);    % Final velocity magnitude [km/s]
% r_final = sqrt(xout(first)^2+yout(first)^2+zout(first)^2);   % Final radius magnitude [km]
%t_final = t(first)/86400;                % Time on transit [days]
fprintf('Desired radius at rendezvous: %.4f [AU]\n',r_f_des/AU)
fprintf('Desired velocity at rendezvous: %.2f [km/s]\n',V_cap)

%equation with mass flow rate and Isp and whatnot
%use the total time burning during outbound trip
mp_o = m(1)-700;              % Initial mass of fuel [kg]
mp_spent = mdotavg*burntime; %9700-m(length(t));     % Mass of fuel spent [kg]
mp_left = 9000-mp_spent;    % Mass of fuel left [kg]
fprintf('Initial fuel mass: %.2f [kg]\n\n',mp_o)

% fprintf('Velocity at rendezvous: %.2f [km/s]\n',v_final)
% fprintf('Radius at rendezvous: %.4f [AU]\n',r_final/AU)
fprintf('Propellant mass left: %.2f [kg]\n',mp_left)
fprintf('Time on transit: %.2f [days]\n',outbound_days)
% 
% coe = coe_from_sv([xout(end) yout(end) zout(end)],[vxout(end) vyout(end) vzout(end)],u_s);
% h = coe(1); e = coe(2); RA = coe(3); incl = coe(4); w = coe(5); TA = coe(6); a = coe(7);
% A = [e a h];
% 
dV_p_cap = sqrt(V_cap^2+(2*u_m/r_p_cap))-sqrt(u_m*(1+e_cap)/r_p_cap);

% e_cap_actual = ((v_final^2)*r_p_cap/u_m)-1;      % Actual eccentricity of capture orbit
% dV_p_cap_actual = sqrt(v_final^2+(2*u_m/r_p_cap))-sqrt(u_m*(1+e_cap_actual)/r_p_cap);   % Actual dV required to ellipticize [km/s]

% fprintf('New eccentricity of Martian elliptical orbit based off of hyperbolic excess velocity: %.4f\n',e_cap_actual)
% fprintf('Actual dV required at perigee: %.3f [km/s]\n',dV_p_cap_actual)
% 
% MR = exp(-dV_p_cap_actual/(850*g_o));
% m_req = MR*mp_left;
% m_left_cap = mp_left-m_req; 
% 
% fprintf('Approximate fuel remaining after capture: %.2f [kg]\n',m_left_cap)

fprintf('Interplanetary Transit Delta-V: %.3f [km/s]\n',del_v_out)

%%

% 3D Simulation
% Outsource visualization magic to Kevin?
% The variables here are currently from four.m
% You can find four.m on GitHub  to see how this plot function works.
% Unsure whether this or the previously way we were plotting is more
% efficient. We'll have to experiment (or replace it with Kevin magic).

% posit=plot3(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,x1,y1,z1,'ob',x2,y2,z2,'or',x3,y3,z3,'oy',x4,y4,z4,'ok');
positfull = plot3(state_e(:,1),state_e(:,2),state_e(:,3),'-ob',state_m(:,1),state_m(:,2),state_m(:,3),'-or',state_p(:,1),state_p(:,2),state_p(:,3),'-xk',state_l(:,1),state_l(:,2),state_l(:,3),'-xc',state_d(:,1),state_d(:,2),state_d(:,3),'-xy',xout,yout,zout,'og',between_e(:,1),between_e(:,2),between_e(:,3),'ob',between_m(:,1),between_m(:,2),between_m(:,3),'or',return_e(:,1),return_e(:,2),return_e(:,3),'ob',return_m(:,1),return_m(:,2),return_m(:,3),'or',xreturn,yreturn,zreturn,'og',state_j(:,1),state_j(:,2),state_j(:,3),'-om',between_j(:,1),between_j(:,2),between_j(:,3),'om',return_j(:,1),return_j(:,2),return_j(:,3),'om',return_j(:,1),return_j(:,2),return_j(:,3),'om');
hold on
[s1,s2,s3] = sphere;
yellow = [1,1,0];
surf(s1*695510*20,s2*695510*20,s3*695510*20) %20 times larger than actual
colormap(yellow)
% plot3(0,0,0,'oy','MarkerSize',10,'MarkerFaceColor','y')
hold off
title('Positions as a function of t (Heliocentric)');
xlabel('x');
ylabel('y');
zlabel('z');
daspect([1 1 1]);
set(gca,'PlotBoxAspectRatio',[1 1 1]);
grid on;
view(-37.5,30);

figure
%issues with sizes of matrices when trying to center plot at Mars

%change this to plot mars as a circle of appropriate radius
positmars = plot3(state_p(:,1)-state_m(:,1),state_p(:,2)-state_m(:,2),state_p(:,3)-state_m(:,3),'-xk',state_d(:,1)-state_m(:,1),state_d(:,2)-state_m(:,2),state_d(:,3)-state_m(:,3),'-xy',xout(ma)-state_m(ma,1),yout(ma)-state_m(ma,2),zout(ma)-state_m(ma,3),'og');
hold on
[s1,s2,s3] = sphere;
red = [1,0,0];
surf(s1*r_m,s2*r_m,s3*r_m)
colormap(red)
hold off
title('Positions as a function of t (Mars-Centered)');
xlabel('x');
ylabel('y');
zlabel('z');
daspect([1 1 1]);
set(gca,'PlotBoxAspectRatio',[1 1 1]);
grid on;
view(-37.5,30);

% state_m(ma,1)
% xout(ma)
% state_m(ma,1)-xout(ma)
% state_m(ma,2)
% yout(ma)
% state_m(ma,2)-yout(ma)
% state_m(ma,3)
% zout(ma)
% state_m(ma,3)-zout(ma)
% sqrt((state_m(ma,1)-xout(ma))^2+(state_m(ma,2)-yout(ma))^2+(state_m(ma,3)-zout(ma))^2)

% figure
% posit = plot3(state_e(:,1),state_e(:,2),state_e(:,3),'-ob',state_m(:,1),state_m(:,2),state_m(:,3),'-or',state_j(:,1),state_j(:,2),state_j(:,3),'-om',state_p(:,1),state_p(:,2),state_p(:,3),'-xk',state_l(:,1),state_l(:,2),state_l(:,3),'-xc',state_d(:,1),state_d(:,2),state_d(:,3),'-xy',xout,yout,zout,'og',0,0,0,'oy');
% title('Positions as a function of t (Heliocentric)');
% xlabel('x');
% ylabel('y');
% zlabel('z');
% daspect([1 1 1]);
% set(gca,'PlotBoxAspectRatio',[1 1 1]);
% grid on;
% view(-37.5,30);
% 
% rd = length(t);
% for rd=1:i/300:i
%     
%     rd=floor(rd);
%     
%     % reset all values for position and velocity, then redraw
%     
%     set(posit(1),'XData',state_e(1:rd,1));
%     set(posit(1),'YData',state_e(1:rd,2));
%     set(posit(1),'ZData',state_e(1:rd,3));
%     set(posit(2),'XData',state_m(1:rd,1));
%     set(posit(2),'YData',state_m(1:rd,2));
%     set(posit(2),'ZData',state_m(1:rd,3));
%     set(posit(3),'XData',state_j(1:rd,1));
%     set(posit(3),'YData',state_j(1:rd,2));
%     set(posit(3),'ZData',state_j(1:rd,3));
%     set(posit(4),'XData',state_p(1:rd,1));
%     set(posit(4),'YData',state_p(1:rd,2));
%     set(posit(4),'ZData',state_p(1:rd,3));
%     set(posit(5),'XData',state_l(rd,1));
%     set(posit(5),'YData',state_l(rd,2));
%     set(posit(5),'ZData',state_l(rd,3));
%     set(posit(6),'XData',state_d(rd,1));
%     set(posit(6),'YData',state_d(rd,2));
%     set(posit(6),'ZData',state_d(rd,3));
%     set(posit(7),'XData',xout(1:rd));
%     set(posit(7),'YData',yout(1:rd));
%     set(posit(7),'ZData',zout(1:rd));
%     set(posit(8),'XData',0);
%     set(posit(8),'YData',0);
%     set(posit(8),'ZData',0);
%     
%     drawnow
% end