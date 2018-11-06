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
 
%% Hohmann Transfer with desired parking orbits
% This calculation finds the total delta V needed to complete a Hohmann
% transfer to Jupiter with a circular orbit of 10 radi starting in an
% circular orbit around earth at an altitude of 390 km
 
% Earth Calculation
V_e = sqrt(u_s/R_e);  % Velocity around sun [km/s]
dV_D = sqrt(u_s/R_e)*(sqrt((2*R_J)/(R_e+R_J))-1);   % Initial burn needed for transfer orbit
 
% Jupiter Calculation
dV_A = sqrt(u_s/R_J)*(1-sqrt(2*R_e/((R_e*(R_e+R_J)))));   %Velocity needed at apogee to circularize orbit
 
%% Propagating Jupiter orbits 
figure()
hold on
 
r_focus_earth = a_ES-r_p_ES;      % Location of focus in Earth orbit [km]
r_focus_jupiter = a_JS-r_p_JS;    % Location of focus in Jupiter orbit [km]
r_focus_mars = a_MS-r_p_MS;       % Location of focus in Mars orbit [km]
r_focus_moon = a_L-r_p_L;         % Location of focus in Moon orbit [km]
r_launch = (r_focus_earth+a_ES+42146); % Radius of departure burn from  [km]
r_500 = (r_focus_earth+a_ES+ri);   % Launching from z = 300 [km]
 
earth = ellipse(a_ES/AU,b_ES/AU,0,r_focus_earth/AU,0,'b'); % Earth Orbit
jupiter = ellipse(a_JS/AU,b_JS/AU,0,r_focus_jupiter/AU,0,'g'); % Jupiter Orbit
mars = ellipse(a_MS/AU,b_MS/AU,0,r_focus_mars/AU,0,'r');
V_excessF9 = sqrt(47.94);  % Excess velocity (C3~48) available from Falcon 9 [km/s]
V_excessANT = sqrt(23);    % Excess velocity (C3~23) available from Antares [km/s]

% Trajectory with only C3 Antares
V_i = [0 (V_e+V_excessANT) 0];        % Velocity of earth + C3
r_launch_vector = [r_launch 0 0];  % Radius of departure burn from 500 [km]
[RadiusC3, VelocityC3, tC3] = prop_sun(r_launch_vector,V_i,[0:360]); % Generate hypothetical orbit
plot(RadiusC3(:,1)/AU,RadiusC3(:,2)/AU,'k--');  % Axes scaled in AU
 
% 1) Departure burn with dV = +5 [km/s]
V1 = [0 (V_e+V_excessANT+5) 0];                % Total Velocity [km/s]
[Radius5, Velocity5, t5] = prop_sun(r_launch_vector,V1,[0:143]); % Generating orbit
plot(Radius5(:,1)/AU,Radius5(:,2)/AU,'c');  % Axes scaled in AU

% 1) Braking Burn at Jupiter with dV = -1 [km/s]
Vbo = Velocity5(141,:); % Arbitrary point to begin braking [km/s]
Rbo = Radius5(141,:);   % Arbitrary point to begin braking [km]
Vb_hat = Vbo/norm(Vbo); % Unit vector in direction of motion [km/s]
Vbrake = Vbo-Vb_hat*1;  % Resulting velocity vector after braking **5 KM/S ** in tangential direction [km/s]
[Radius2, Velocity2] = prop_sun(Rbo,Vbrake,[0:360]); % Generating new orbit
plot(Radius2(:,1)/AU,Radius2(:,2)/AU,'c--'); % Axes scaled in AU
 
% 2) Departure burn with dV = +8 [km/s]
Vburn1 = [0 (V_e+V_excessANT+8) 0];            % Total Velocity [km/s]
[Radius8, Velocity8, t8] = prop_sun(r_launch_vector,Vburn1,[0:122]); % Generating orbit
plot(Radius8(:,1)/AU,Radius8(:,2)/AU,'m');  % Axes scaled in AU
 
% Braking Burn with dV = -3.5 [km/s]
Vbo2 = Velocity8(122,:);     % Arbitrary point to begin braking [km/s]
Rbo2 = Radius8(122,:);       % Arbitrary point to begin braking [km]
Vb_hat2 = Vbo2/norm(Vbo2);   % Unit vector in direction of motion [km/s]
Vbrake2 = Vbo2-Vb_hat2*3.5; % Resulting velocity vector after braking **10.5 KM/S** in tangential direction [km/s]
[Radius4, Velocity4] = prop_sun(Rbo2,Vbrake2,[0:360]); % Generating new orbit 
plot(Radius4(:,1)/AU,Radius4(:,2)/AU,'m--'); % Axes scaled in AU
 
% 3) Departure burn with dV = +20 [km/s]
Vburn1 = [0 (V_e+V_excessANT+20) 0];            % Total Velocity [km/s]
[Radius20, Velocity20, t20] = prop_sun(r_launch_vector,Vburn1,[0:97]); % Generating orbit
plot(Radius20(:,1)/AU,Radius20(:,2)/AU,'y');  % Axes scaled in AU
 
% Hohmann Transfer 
r1_hoh = [R_e 0 0];       % Launching from 
v1_hoh = [0 V_e+dV_D 0];    % Velocity vector to begin Hohmann transfer [km/s]
[r2_hoh, v2_hoh, t_hoh] = prop_sun(r1_hoh,v1_hoh,[0:180]); % Generate transfer orbit
plot(r2_hoh(:,1)/AU,r2_hoh(:,2)/AU,'k') % Axes scaled in AU
 
r3_hoh = r2_hoh(180,:);          % Begin to circularize orbit at apogee of transfer orbit [km]
dV_c = dV_A+v2_hoh(180,2);    % Magnitude of burn needed to circularize orbit [km/s]
v3_hoh = [0 -dV_A 0];         % Initial circularized orbit velocity at apogee [km/s]
[r4_hoh, v4_hoh] = prop_sun(r3_hoh, v3_hoh,[-180:180]); % Generate final Hohmann orbit
plot(r4_hoh(:,1)/AU, r4_hoh(:,2)/AU,'k') % Axes scaled in AU

%%%%%%%%% DISCREPENCIES DUE TO ELLEPTICAL NATURE OF JUPITER ORBIT BEING
%%%%%%%%% PLOTTED OFF CENTER WITH THE SUN AT FOCUS
 
sun = circle(0,0,30000000/AU);  % Sun Location **NOT TO SCALE**
axis square
suptitle('Jesse Owens Trajectories')
xlabel('AU')
ylabel('AU')
legend('Earth','Jupiter','Mars','C3 Antares','+5 dV Burn','-1 dV Burn','+8 dV Burn'...
    ,'-3.5 dV Burn','+20 dV Burn','Hohmann','location','northeast')
grid on
hold off

%% Propagating Mars orbits 
figure()
hold on

r_focus_earth = a_ES-r_p_ES;      % Location of focus in Earth orbit [km]
r_focus_mars = a_MS-r_p_MS;       % Location of focus in Mars orbit [km]
r_launch = (r_focus_earth+a_ES+r_o); % Radius of departure burn from GEO [km]
r_300 = (r_focus_earth+a_ES+ri+200);   % Launching from z = 500 [km]
r_launch_vector = [r_500 0 0];  % Radius of departure burn from 500 [km] 

earth = ellipse(a_ES/AU,b_ES/AU,0,r_focus_earth/AU,0,'b'); % Earth Orbit
mars = ellipse(a_MS/AU,b_MS/AU,0,r_focus_mars/AU,0,'r');
V_excessANT = sqrt(23);    % Excess velocity (C3~23) available from Antares [km/s]

% Departure burn with dV = +2 [km/s]
V2 = [0 (V_e+V_excessANT+2) 0];                % Total Velocity [km/s]
[Radius2, Velocity2, T2] = prop_sun(r_launch_vector,V2,(0:76)); % Generating orbit 76
plot(Radius2(:,1)/AU,Radius2(:,2)/AU,'c','LineWidth',1);  % Axes scaled in AU

% Braking Burn at Mars with dV = -10 [km/s]
V2_bo = Velocity2(76,:); % Arbitrary point to begin braking [km/s]
R2_b = Radius2(76,:);   % Arbitrary point to begin braking [km]
V2_bhat = V2_bo/norm(V2_bo); % Unit vector in direction of motion [km/s]
V2_b = V2_bo-V2_bhat*10;  % Resulting velocity vector after braking **10 KM/S** in tangential direction [km/s]
[R2_new, V2_new, T2_trash] = prop_sun(R2_b,V2_b,[0:9]); % Generating new orbit
[R2_trash, V2_trash, T2_new] = prop_sun(R2_b,V2_b,[0:360]); % Generating new orbit
plot(R2_new(:,1)/AU,R2_new(:,2)/AU,'c--'); % Axes scaled in AU

% Departure burn with dV = +5 [km/s]
V5 = [0 (V_e+V_excessANT+5) 0];            % Total Velocity [km/s]
[Radius5, Velocity5, T5] = prop_sun(r_launch_vector,V5,(0:65)); % Generating orbit 65
plot(Radius5(:,1)/AU,Radius5(:,2)/AU,'m','LineWidth',1);  % Axes scaled in AU

% Braking Burn at Mars with dV = -15 [km/s]
V5_bo = Velocity5(65,:); % Arbitrary point to begin braking [km/s]
R5_b = Radius5(65,:);   % Arbitrary point to begin braking [km]
V5_bhat = V5_bo/norm(V5_bo); % Unit vector in direction of motion [km/s]
V5_b = V5_bo-V5_bhat*15;  % Resulting velocity vector after braking **15 KM/S** in tangential direction [km/s]
[R5_new, V5_new, T5_trash] = prop_sun(R5_b,V5_b,[0:11]); % Generating new orbit
[R5_trash, V5_trash, T5_new] = prop_sun(R5_b,V5_b,[0:360]); % Generating new orbit
plot(R5_new(:,1)/AU,R5_new(:,2)/AU,'m--'); % Axes scaled in AU
 
% Departure burn with dV = +10 [km/s]
Vburn10 = [0 (V_e+V_excessANT+10) 0];            % Total Velocity [km/s]
[Radius10, Velocity10, T10] = prop_sun(r_launch_vector,Vburn10,[0:50]); % Generating orbit
plot(Radius10(:,1)/AU,Radius10(:,2)/AU,'g','LineWidth',1);  % Axes scaled in AU

% Braking Burn at Mars with dV = -18.5 [km/s]
V10_bo = Velocity10(50,:); % Arbitrary point to begin braking [km/s]
R10_b = Radius10(50,:);   % Arbitrary point to begin braking [km]
V10_bhat = V10_bo/norm(V10_bo); % Unit vector in direction of motion [km/s]
V10_b = V10_bo-V10_bhat*18.5;  % Resulting velocity vector after braking **18.5 KM/S** in tangential direction [km/s]
[R10_new, V10_new, T10_trash] = prop_sun(R10_b,V10_b,[0:23]); % Generating new orbit
[R10_trash, V10_trash, T10_new] = prop_sun(R10_b,V10_b,[0:360]); % Generating new orbit
plot(R10_new(:,1)/AU,R10_new(:,2)/AU,'g--'); % Axes scaled in AU

% Time Calculations
T_tot2 = T2(end)+(T2_new(85+76)-T2_new(76+76));
T_tot5 = T5(end)+T5_new(76+55)-T5_new(65+55);
T_tot10 = T10(end)+T10_new(50+23)-T5_new(27+23);

sun = circle(0,0,10000000/AU);  % Sun Location **NOT TO SCALE**
axis equal
suptitle('Jesse Owens Trajectories to Mars - Antares')
xlabel('AU') 
ylabel('AU')
legend('Earth','Mars',['+2 km/s   \DeltaV \rightarrow '...
    num2str(T_tot2) ' Days'], '-10 km/s',['+5 km/s   \DeltaV \rightarrow '...
    num2str(T_tot5) ' Days'], '-15 km/s',['+10 km/s \DeltaV \rightarrow '...
    num2str(T_tot10) ' Days']', '-18.5 km/s', 'location','northeast')
grid on
hold off
print -dpdf Jesse_Owens_Mars_Trajectories

figure()
hold on

earth = ellipse(a_ES/AU,b_ES/AU,0,r_focus_earth/AU,0,'b'); % Earth Orbit
mars = ellipse(a_MS/AU,b_MS/AU,0,r_focus_mars/AU,0,'r');
% Trajectory with only C3 Antares
% V_i = [0 (V_e+V_excessANT) 0];        % Velocity of earth + C3
% r_launch_vector = [r_launch 0 0];  % Radius of departure burn from GEO [km]
% [RadiusC3, VelocityC3, tC3] = prop_sun(r_launch_vector,V_i,[0:360]); % Generate hypothetical orbit
% plot(RadiusC3(:,1)/AU,RadiusC3(:,2)/AU,'r--');  % Axes scaled in AU

% Departure burn with dV = +0.5 [km/s]
Vburn0_5 = [0 (V_e+V_excessANT+0.5) 0];            % Total Velocity [km/s]
[Radius0_5, Velocity0_5, t0_5] = prop_sun(r_launch_vector,Vburn0_5,(0:90)); % Generating orbit
plot(Radius0_5(:,1)/AU,Radius0_5(:,2)/AU,'y','LineWidth',1);  % Axes scaled in AU

% Departure burn with dV = +1 [km/s]
Vburn1 = [0 (V_e+V_excessANT+1) 0];            % Total Velocity [km/s]
[Radius1, Velocity1, t1] = prop_sun(r_launch_vector,Vburn1,(0:88)); % Generating orbit
plot(Radius1(:,1)/AU,Radius1(:,2)/AU,'b','LineWidth',1);  % Axes scaled in AU

% Departure burn with dV = +1.5 [km/s]
Vburn1_5 = [0 (V_e+V_excessANT+1.5) 0];            % Total Velocity [km/s]
[Radius1_5, Velocity1_5, t1_5] = prop_sun(r_launch_vector,Vburn1_5,(0:85)); % Generating orbit
plot(Radius1_5(:,1)/AU,Radius1_5(:,2)/AU,'g','LineWidth',1);  % Axes scaled in AU

% Departure burn with dV = +2 [km/s]
Vburn2 = [0 (V_e+V_excessANT+2) 0];            % Total Velocity [km/s]
[Radius2, Velocity2, t2] = prop_sun(r_launch_vector,Vburn2,(0:83)); % Generating orbit
plot(Radius2(:,1)/AU,Radius2(:,2)/AU,'k','LineWidth',1);  % Axes scaled in AU

% Departure burn with dV = +2.5 [km/s]
Vburn2_5 = [0 (V_e+V_excessANT+2.5) 0];            % Total Velocity [km/s]
[Radius2_5, Velocity2_5, t2_5] = prop_sun(r_launch_vector,Vburn2_5,(0:81)); % Generating orbit
plot(Radius2_5(:,1)/AU,Radius2_5(:,2)/AU,'b','LineWidth',1);  % Axes scaled in AU

% Departure burn with dV = +3 [km/s]
Vburn3 = [0 (V_e+V_excessANT+3) 0];            % Total Velocity [km/s]
[Radius3, Velocity3, t3] = prop_sun(r_launch_vector,Vburn3,(0:79)); % Generating orbit
plot(Radius3(:,1)/AU,Radius3(:,2)/AU,'g','LineWidth',1);  % Axes scaled in AU

% Departure burn with dV = +3.5 [km/s]
Vburn3_5 = [0 (V_e+V_excessANT+3.5) 0];            % Total Velocity [km/s]
[Radius3_5, Velocity3_5, t3_5] = prop_sun(r_launch_vector,Vburn3_5,(0:78)); % Generating orbit
plot(Radius3_5(:,1)/AU,Radius3_5(:,2)/AU,'y','LineWidth',1);  % Axes scaled in AU

% Departure burn with dV = +4 [km/s]
Vburn4 = [0 (V_e+V_excessANT+4) 0];            % Total Velocity [km/s]
[Radius4, Velocity4, t4] = prop_sun(r_launch_vector,Vburn4,(0:76)); % Generating orbit
plot(Radius4(:,1)/AU,Radius4(:,2)/AU,'b','LineWidth',1);  % Axes scaled in AU

% Departure burn with dV = +4.5 [km/s]
Vburn4_5 = [0 (V_e+V_excessANT+4.5) 0];            % Total Velocity [km/s]
[Radius4_5, Velocity4_5, t4_5] = prop_sun(r_launch_vector,Vburn4_5,(0:75)); % Generating orbit
plot(Radius4_5(:,1)/AU,Radius4_5(:,2)/AU,'g','LineWidth',1);  % Axes scaled in AU

% Departure burn with dV = +5 [km/s]
V5 = [0 (V_e+V_excessANT+5) 0];                % Total Velocity [km/s]
[Radius5, Velocity5, t5] = prop_sun(r_launch_vector,V5,(0:74)); % Generating orbit
plot(Radius5(:,1)/AU,Radius5(:,2)/AU,'c','LineWidth',1);  % Axes scaled in AU
 
% Departure burn with dV = +5.5 [km/s]
Vburn5_5 = [0 (V_e+V_excessANT+5.5) 0];            % Total Velocity [km/s]
[Radius5_5, Velocity5_5, t5_5] = prop_sun(r_launch_vector,Vburn5_5,(0:73)); % Generating orbit
plot(Radius5_5(:,1)/AU,Radius5_5(:,2)/AU,'k','LineWidth',1);  % Axes scaled in AU

% Departure burn with dV = +6 [km/s]
Vburn6 = [0 (V_e+V_excessANT+6) 0];            % Total Velocity [km/s]
[Radius6, Velocity6, t6] = prop_sun(r_launch_vector,Vburn1,(0:74)); % Generating orbit
plot(Radius6(:,1)/AU,Radius6(:,2)/AU,'b','LineWidth',1);  % Axes scaled in AU

% Departure burn with dV = +6.5 [km/s]
Vburn6_5 = [0 (V_e+V_excessANT+6.5) 0];            % Total Velocity [km/s]
[Radius6_5, Velocity6_5, t6_5] = prop_sun(r_launch_vector,Vburn6_5,(0:71)); % Generating orbit
plot(Radius6_5(:,1)/AU,Radius6_5(:,2)/AU,'g','LineWidth',1);  % Axes scaled in AU

% Departure burn with dV = +7 [km/s]
Vburn7 = [0 (V_e+V_excessANT+7) 0];            % Total Velocity [km/s]
[Radius7, Velocity7, t7] = prop_sun(r_launch_vector,Vburn7,(0:69)); % Generating orbit
plot(Radius7(:,1)/AU,Radius7(:,2)/AU,'y','LineWidth',1);  % Axes scaled in AU

% Departure burn with dV = +7.5 [km/s]
Vburn7_5 = [0 (V_e+V_excessANT+7.5) 0];            % Total Velocity [km/s]
[Radius7_5, Velocity7_5, t7_5] = prop_sun(r_launch_vector,Vburn7_5,(0:69)); % Generating orbit
plot(Radius7_5(:,1)/AU,Radius7_5(:,2)/AU,'b','LineWidth',1);  % Axes scaled in AU

% Departure burn with dV = +8 [km/s]
Vburn8 = [0 (V_e+V_excessANT+8) 0];            % Total Velocity [km/s]
[Radius8, Velocity8, t8] = prop_sun(r_launch_vector,Vburn8,(0:68)); % Generating orbit
plot(Radius8(:,1)/AU,Radius8(:,2)/AU,'m','LineWidth',1);  % Axes scaled in AU
 
% Departure burn with dV = +8.5 [km/s]
Vburn8_5 = [0 (V_e+V_excessANT+8.5) 0];            % Total Velocity [km/s]
[Radius8_5, Velocity8_5, t8_5] = prop_sun(r_launch_vector,Vburn8_5,(0:67)); % Generating orbit
plot(Radius8_5(:,1)/AU,Radius8_5(:,2)/AU,'g','LineWidth',1);  % Axes scaled in AU

% Departure burn with dV = +9[km/s]
Vburn9 = [0 (V_e+V_excessANT+9) 0];            % Total Velocity [km/s]
[Radius9, Velocity9, t9] = prop_sun(r_launch_vector,Vburn9,(0:67)); % Generating orbit
plot(Radius9(:,1)/AU,Radius9(:,2)/AU,'k','LineWidth',1);  % Axes scaled in AU

% Departure burn with dV = +9.5 [km/s]
Vburn9_5 = [0 (V_e+V_excessANT+9.5) 0];            % Total Velocity [km/s]
[Radius9_5, Velocity9_5, t9_5] = prop_sun(r_launch_vector,Vburn9_5,(0:66)); % Generating orbit
plot(Radius9_5(:,1)/AU,Radius9_5(:,2)/AU,'b','LineWidth',1);  % Axes scaled in AU

% Departure burn with dV = +10 [km/s]
Vburn10 = [0 (V_e+V_excessANT+10) 0];            % Total Velocity [km/s]
[Radius10, Velocity10, t10] = prop_sun(r_launch_vector,Vburn10,(0:66)); % Generating orbit
plot(Radius10(:,1)/AU,Radius10(:,2)/AU,'g','LineWidth',1);  % Axes scaled in AU

%%%%%%%%% DISCREPENCIES DUE TO ELLEPTICAL NATURE OF JUPITER ORBIT BEING
%%%%%%%%% PLOTTED OFF CENTER WITH THE SUN AT FOCUS
 
sun = circle(0,0,30000000/AU);  % Sun Location **NOT TO SCALE**
axis equal
suptitle('Jesse Owens Trajectories to Mars - Antares')
xlabel('AU') 
ylabel('AU')
% legend('Earth','Mars','\Delta V = 0.5 km/s','\Delta V = 0.5 km/s',...
%     '\Delta V = 1 km/s','\Delta V = 1.5 km/s','\Delta V = 2 km/s',...
%     '\Delta V = 2.5 km/s','\Delta V = 3 km/s','\Delta V = 3.5 km/s',...
%     '\Delta V = 4 km/s','\Delta V = 4.5 km/s','\Delta V = 5 km/s',...
%     '\Delta V = 5.5 km/s','\Delta V = 6 km/s','\Delta V = 6.5 km/s',...
%     '\Delta V = 7 km/s','\Delta V = 7.5 km/s','\Delta V = 8 km/s',...
%     '\Delta V = 8.5 km/s','\Delta V = 9 km/s','\Delta V = 9.5 km/s',...
%     '\Delta V = 10 km/s','Location','Northwest')
grid on
hold off

% figure()
% hold on 
% Earth = circle(0,0,30000);  % Sun Location **NOT TO SCALE**
% moon = ellipse(a_L,b_L,0,r_focus_moon,0,'k');

%% Mars Fuel 
%{
    NOTE: All times depended on hard coding the angle range of each 
    propagation. The plot reflects this by showing a non linear line. 
    Possible solution to this is to either fine LOBF or solve for exact 
    angle range to reach Mars. 
%}



Isp = 850;             % Specific Impulse [s]
Md = 680.3886;         % Dry mass of spacecraft [kg]
DeltaVs = (0.5:.5:10); % Delta Vs       
g = 0.00981;           % Gravity [km/s^2]

Mw = Md.*exp(DeltaVs./(g*Isp)); % Wet Mass [kg]
Mf = Mw - Md;                   % Mass of Fuel [kg]


t0_5 = t0_5(end); t1 = t1(end); t1_5 = t1_5(end); t2 = t2(end); % Extracting
t2_5 = t2_5(end); t3 = t3(end); t3_5 = t3_5(end); t4 = t4(end); % times from
t4_5 = t4_5(end); t5 = t5(end); t5_5 = t5_5(end); t6 = t6(end); % each departure
t6_5 = t6_5(end); t7 = t7(end); t7_5 = abs(t7_5(end)); t8 = abs(t8(end)); 
t8_5 = abs(t8_5(end)); t9 = abs(t9(end)); t9_5 = abs(t9_5(end)); 
t10 = abs(t10(end)); 

tt=[t0_5 t1 t1_5 t2 t2_5 t3 t3_5 t4 t4_5 t5 t5_5 t6 t6_5 t7 t7_5 t8 t8_5...
    t9 t9_5 t10]; % time vectors 

p = polyfit(DeltaVs,tt,1);
tt1 = polyval(p,DeltaVs);
 
figure()
yyaxis left
plot(DeltaVs,Mf)
title('JO Requirments for Mars')
ylabel('Mass of Fuel per Burn (kg)')
xlabel('Departure Burn \Delta V (km/s)')
yyaxis right
plot(DeltaVs,tt1)
ylabel('Time (days)')
grid on
print -dpdf Jesse_Owens_Mars_Requirments


%% Jupiter Fuel
% Ve = (10:1:70);       % Varying Exit Velocity [km/s]
% V  = 19;              % Delta V [km/s]
% MR = exp(V./Ve);      % Mass ratio [kg]
% Md = 150;             % Dry mass of spacecraft [kg]
% M_fuel = Md.*MR - Md; % Mass of fuel [kg]
%  
% % figure()
% % plot(Ve,M_fuel)
% % title('JO Fuel Requirments for Jupiter')
% % ylabel('Mass of Fuel (kg)')
% % xlabel('Exhaust Velocity (km/s)')
% % grid on







