% Zach Zoloty
% JPL Ephemerides Test
% Does not take multiple bodies into account

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Earth Vector Ephemerides w/r/t the Sun:
% 2458409.500000000 = A.D. 2018-Oct-18 00:00:00.0000 TDB 
%  X = 9.079361611992910E-01 Y = 4.174196284373897E-01 Z =-9.320364656857048E-05
%  VX=-7.374461784272602E-03 VY= 1.561335766240027E-02 VZ= 1.196212616732829E-07
%
% Mars Vector Ephemerides w/r/t the Sun:
% 2458409.500000000 = A.D. 2018-Oct-18 00:00:00.0000 TDB 
%  X = 1.384840563405629E+00 Y =-9.437889789572897E-02 Z =-3.595859583226332E-02
%  VX= 1.486059872073195E-03 VY= 1.515696619906733E-02 VZ= 2.811372207049188E-04
%
% (used JPL Horizons web-app: https://ssd.jpl.nasa.gov/horizons.cgi)
%
% Ephemeride position vector given in AU.
% Ephemeride velocity vector given in AU/day.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;

%Define Constants
au_m = 1.496e11; %constant to convert the au of ephemerides to meters
day_s = 86400; %constant to convert the days of ephermides to seconds

mu = 1.327124e20; %m^3/s^2

% Earth Ephemerides
Xe = 9.079361611992910e-1;
Ye = 4.174196284373897e-1;
Ze = -9.320364656857048e-5;
VXe = -7.374461784272602e-3;
VYe = 1.561335766240027e-2;
VZe = 1.196212616732829e-7;

r0e = [Xe*au_m; Ye*au_m; Ze*au_m]; %m
v0e = [VXe*au_m/day_s; VYe*au_m/day_s; VZe*au_m/day_s]; %m/s
x0e = [r0e; v0e]; %state variables
Ee = norm(v0e)^2/2 - mu/norm(r0e); %km^2/s^2
ae = -mu/(2*Ee); %km
Te = 2*pi*sqrt(ae^3/mu); %s

% Mars Ephemerides
Xm = 1.384840563405629;
Ym = -9.437889789572897e-2;
Zm = -3.595859583226332e-2;
VXm = 1.486059872073195e-3;
VYm = 1.515696619906733e-2;
VZm = 2.811372207049188e-4;

r0m = [Xm*au_m; Ym*au_m; Zm*au_m]; %m
v0m = [VXm*au_m/day_s; VYm*au_m/day_s; VZm*au_m/day_s]; %m/s
x0m = [r0m; v0m]; %state variables
Em = norm(v0m)^2/2 - mu/norm(r0m); %km^2/s^2
am = -mu/(2*Em); %km
Tm = 2*pi*sqrt(am^3/mu); %s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define Equations of Motion and Solve
EoM = @(t,x) [zeros(3,3), eye(3); -(mu/norm(x(1:3,1))^3)*eye(3), zeros(3,3)]*x;
options = odeset('AbsTol', 1e-9, 'RelTol', 1e-6);

     %Earth
[time_e, state_e] = ode45(EoM,[0,2*Te],x0e,options);

     %Mars
[time_m, state_m] = ode45(EoM,[0,2*Tm],x0m,options);
     
%Plot resulting orbits
figure(1)
hold on;
for j = 1:length(time_m)
    plot3(state_e(1:j,1), state_e(1:j,2), state_e(1:j,3), 'b', 'linewidth', 1);
    hold on;
    plot3(state_m(1:j,1), state_m(1:j,2), state_m(1:j,3), 'r', 'linewidth', 1);
    
    hold on;
    plot3(state_e(j,1), state_e(j,2), state_e(j,3), 'go', 'Markersize', 8, 'markerfacecolor', 'b');
    plot3(state_m(j,1), state_m(j,2), state_m(j,3), 'yo', 'Markersize', 6, 'markerfacecolor', 'r');
    
    plot3(0, 0, 0, 'ro', 'Markersize', 16, 'markerfacecolor', 'y');
    
    axis([-2.5e11 2.5e11 -2.25e11 2.5e11 -0.2e11 0.2e11]);
    set(gca, 'fontsize', 14, 'fontweight', 'bold');
    view(100,24);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    grid on;
    hold off;
    getframe(figure(1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%