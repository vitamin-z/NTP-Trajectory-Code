% Zach Zoloty
% AAE 5626
% Project 2
% 11/27/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc

%Define Constants
%Obtained from Table A.1 and A.2
Rearth = 149.6e6; %km
Rmars = 227.9e6; %km
mu_earth = 398600; %km^3/s^2
mu_mars = 42828; %km^3/s^2
mu_sun = 132.712e9; %km^3/s^2
Tearth = 365.256; %days
Tmars = 1.881*365; %days
parking1 = 500; %km
mearth = 5.974e24; %kg
mmars = 641.9e21; %kg
msun = 1.989e30; %kg
radearth = 6378; %km
radmars = 3396; %km
rpmars = 2*radmars; %km

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Part 1
%Synodic Period between Earth and Mars
Tsyn = (Tearth*Tmars)/abs(Tearth-Tmars); %Equation from Example 8.1
    %Equation 8.10

%Phase Angle at departure from Earth
t12 = (pi/sqrt(mu_sun))*((Rearth+Rmars)/2)^(3/2); %Equation from Example 8.2
    %Equation 8.11
t12day = t12/(60*60*24);
n_earth = 2*pi/Tearth; %rad/day
n_Mars = 2*pi/Tmars; %rad/day

phi0 = pi - n_Mars*t12day; %rad
    %Equation 8.12

%Phase Angle at arrival at Mars
phif = pi - n_earth*t12day; %rad
    %Equation 8.13

%Time of flight to execute outbound Hohmann Transfer
a_t = 0.5*(Rearth+Rmars); %km
    %Eq9.68 in SMAD

transfertime = pi*(a_t^(3/2))*(mu_sun^(-1/2))/(60*60*24); %days
    %This is the same equation as 8.11


%Wait time required after arrival before a return Hohmann Transfer
%trajectory can be initiated
    %Equation 8.16 because n1>n2
for i=1:101
    N = i-1;
    twait = (-2*phif-2*pi*N)/(n_Mars-n_earth);
    
    if twait > 0
        break
    end
end

%Delta-V required between the departure planet's orbital velocity and
%solar-system Hohmann transfer ellipse
delta_vL = sqrt(mu_sun)*(sqrt((2/Rearth)-(1/a_t))-sqrt(1/Rearth)); %Eq9.69 in SMAD
    %Same as Equation 8.3

%Orbital Characteristics: semi-major axis, eccentricity, angular momentum
%of the Hohmann Transfer orbit
%semi-major axis already solved for, given by a_t

e_transfer = (Rmars-Rearth)/(Rmars+Rearth); %substituted for perihelion and aphelion
h_transfer = sqrt(Rearth*mu_sun*(1+e_transfer)); %km^2/s

%Hyperbolic Excess Velocity needed for the departure hyperbola in the
%planet's sphere of influence
%First, calculate the sphere of influence of the Earth
rsoi_earth = Rearth*(mearth/msun)^(2/5);

vdepart_excess = sqrt(mu_sun/Rearth)*(sqrt(2*Rmars/(Rearth+Rmars))-1); %km/s
    %Equation 8.3

%Beta angle required for the departure burn to hop from the parking orbit
%to the hyperbolic departure orbit
beta_dephyp = acosd(1/(1+(radearth+parking1)*(vdepart_excess^2)/mu_earth)); %degrees
    %Equation 8.43

%Orbital Characteristics: semi-major axis, eccentricity, angular momentum
%of the hyperbolic departure orbit

e_hyp = 1 + (parking1+radearth)*vdepart_excess^2/mu_earth;%1/cosd(beta_dephyp);

%h_hyp = (parking1+radearth)*sqrt(vdepart_excess^2 + 2*mu_earth)/(parking1+radearth);%sqrt((radearth+parking1)*mu_earth*(1+e_hyp)); %km^2/s
h_hyp = sqrt((radearth+parking1)*mu_earth*(1+e_hyp)); %km^2/s

a_hyp = (radearth+parking1)/(e_hyp-1); %km

%%%%%%%%%%%%%%%%%Plot the hyperbolic departure orbit and parking orbit
t = linspace(0,2*pi);

Xpark = (parking1+radearth)*cos(t);
Ypark = (parking1+radearth)*sin(t);

%Hyperbolic trajectory was giving me trouble, so I'll draw that in
%afterwards

figure(1);
plot(0,0,'ob',Xpark,Ypark,'g');
%grid on;
legend('Earth','Parking Orbit');
xlabel('X-Position [km]');
ylabel('Y-Position [km]');
title('Geocentric Parking Orbit and Hyperbolic Departure Orbit');

%%%%%%%%%%%%%%%%%%Plot the solar-system Hohmann transfer ellipse
btrans = a_t*sqrt(1-e_transfer^2);
Xtrans = a_t*cos(t);
Ytrans = btrans*sin(t);

figure(2);
plot(0,0,'ok',Xtrans,Ytrans);
%grid on
%axis([-2.5*10^9 2.5*10^9 -2.5*10^9 2.5*10^9]);
legend('Sun','Hohmann Transfer Orbit');
xlabel('X-Position [km]');
ylabel('Y-Position [km]');
title('Solar-System Hohmann Transfer Ellipse');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Part 2
%Enter into an orbit which has a periapsis radius of 5 times the radius of
%the arrival planet using either a leading or trailing approach
%Excess Inbound Hyperbolic Velocity upon arrival into planet's SOI
%rsoi_Mars = RMars*(mMars/msun)^(2/5);

vinbound_excess = sqrt(mu_sun/Rmars)*(1-sqrt(2*Rearth/(Rearth+Rmars))); %km/s
    %Equation 8.4; same as 8.50

%Eccentricity of the capture orbit that minimizes the delta-v required to
%gain this periapsis radius
e_capture = ((2*mu_mars/(rpmars*vinbound_excess^2))-1)/((2*mu_mars/(rpmars*vinbound_excess^2))+1);
    %rearranged Equation 8.67 (for optimal periapse radius as far as fuel
    %expenditure is concerned)

%Periapsis Velocities prior to and after capture burn, and the Resulting
%Delta-V
vp_hyp = sqrt(vinbound_excess^2 + 2*mu_mars/rpmars);
    %Equation 8.58
vp_cap = sqrt(mu_mars*(1+e_capture)/rpmars);
    %Equation 8.59
delta_v_cap = vp_hyp - vp_cap;
    %Equation 8.60
    %can also be found using Equation 8.70
    %(vinbound_excess*sqrt((1-e_capture)/2)

%Aim Distance (Delta)
delta_aim = rpmars*sqrt(2/(1-e_capture)); %km
    %Equation 8.71

%Orbital Elements: Period, Angular Momentum, Semi-major axis of the
%resulting capture orbit
a_cap = rpmars/(1-e_capture); %km

T_cap = (2*pi*a_cap^(3/2))/sqrt(mu_mars)/(60*60*24); %days

h_cap = sqrt(rpmars*mu_mars*(1+e_capture)); %km^2/s

%Periapse Angle (Beta) for this capture orbit
beta = acosd(1/(1+rpmars*vinbound_excess^2/mu_mars)); %degrees
    %Equation 8.43

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Oops! Your engine didn't fire!
%%%%%%%%%%%%%%%Excess Outbound Hyperbolic Velocity leaving the planet's SOI
%For this problem, the spacecraft will be flying a leading approach.
vperp1 = (mu_sun/h_transfer)*(1+e_transfer*cosd(180)); %km/s
    %Equation 8.78a
vrad1 = (mu_sun/h_transfer)*e_transfer*sind(180); %km/s
    %Equation 8.78b

v1v = norm([vperp1 -1*vrad1]);
    %Equation 8.75

vMars = sqrt(mu_sun/Rmars);%vdepart_excess - v1v;
    %Equation 8.1
    
alpha1 = acosd(vperp1/v1v);
    %rearranged Equation 8.76a

vexcess_v = v1v*cosd(alpha1) - vMars;
    %Equation 8.81a
vexcess_s = v1v*sind(alpha1);
    %Equation 8.81b

vexcess = norm([vexcess_v vexcess_s]);
    %Equation 8.80

%%%%%%%%%%%The Turn angle Delta through which the spacecraft will fly
h_flyby = rpmars*sqrt(vexcess^2 + 2*mu_mars/rpmars);
    %Equation 8.39

e_flyby = 1 + (rpmars*vexcess^2)/mu_mars;
    %Equation 8.38

delta_turn = 2*asind(1/e_flyby); %degrees; from Equation in Example 8.6

%Resulting Heliocentric Velocity in radial/tangential solar coordinates
%immediately after flyby
phi1 = atand(vexcess_s/vexcess_v); %from Equation in Example 8.6
phi2 = phi1 + delta_turn; %from Equation in Example 8.6
vperp2 = vMars + vexcess*cosd(phi2); %from Equation in Example 8.6
vrad2 = -1*vexcess*sind(phi2); %from Equation in Example 8.6

%Characterize the resulting solar orbit: eccentricity, semi-major axis,
%perihelion, period, and angular momentum
h_final = Rmars*vperp2;
    %Equation 8.90

e_final = 1 - (h_final^2)/(Rmars*mu_sun);
    %rearranged Equation 8.91

a_final = Rmars/(1+e_final);

rp_final = a_final*(1-e_final);
    %can also use Equation 2.50

T_final = (2*pi*a_final^(3/2))/sqrt(mu_sun)/(60*60*24);

%How does this compare to the original Hohmann Transfer orbit?
%Specifically, has the solar velocity increased or decreased?
%vis-viva
vpt = sqrt(mu_sun*((2/(a_t*(1-e_transfer)))-1/a_t));
vat = sqrt(mu_sun*((2/(a_t*(1+e_transfer)))-1/a_t));

vph = sqrt(mu_sun*((2/rp_final-1/a_final)));
vah = sqrt(mu_sun*((2/(a_final*(1+e_final)))-1/a_final));

%Plot the Earth's Orbit, the arrival planet's orbit, the original Hohman
%transfer orbit, and the final heliocentric orbit of the spacecraft
Xearth = Rearth*cos(t);
Yearth = Rearth*sin(t);

XMars = Rmars*cos(t);
YMars = Rmars*sin(t);

% Xfinal = 
% Yfinal = 

figure(3);
plot(0,0,'ok',Xearth,Yearth,'b',XMars,YMars,'b',Xtrans,Ytrans,'r');
xlabel('X-Position [km]');
ylabel('Y-Position [km]');
title('Orbits of Earth and Mars with the Spacecraft''s Final and Transfer Orbits');