% Approximation of Three Body Motion

%Initial ephemerides are taken at 12AM on Oct18, 2018

%Need to take out current pos/vel/accel equations and replace with
%equations of motion so that the spacecraft doesn't just travel in a
%straight line

clc;

% Universal Constants

G=6.674e-11; %N/kg^2*m^2

% initial conditions

cc=input('Case?');

switch cc
    case 1 %earth/mars/sun/phobos
        %Earth
        m1=5.97219e24; %kg
        x1=1.358519213625654E+11; %m
        y1=6.135912462718222E+10; %m
        z1=-3.183658747591078E+6; %m
        vx1=-1.275541057982471E+4; %m/s
        vy1=2.702943044653523E+4;
        vz1=-1.216464121185368E-1;
        
        %Mars
        m2=6.4171e23; %kg
        x2=2.071691995444704E+11;
        y2=-1.411888216421376E+10;
        z2=-5.379329369868487E+9;
        vx2=2.573048525403522E+3;
        vy2=2.624363275060584E+4;
        vz2=4.867769628700440E+2;
        
        %Sun
        m3=1.989e30; %kg
        x3=0;
        y3=0;
        z3=0;
        vx3=0;
        vy3=0;
        vz3=0;
        
        %Phobos
        m4=10e20; %kg
        x4=2.071356927243748E+11;
        y4=-1.303814318469472E+10;
        z4=-5.386790737800404E+09;
        vx4=3.584369434387757E+03;
        vy4=2.445376574911518E+04;
        vz4=-1.546067254295771E+02;
        
        endtime=4*365.25*86400;
        i=1e5;
        
        
    case 2 %earth/mars/sun/spacecraft
        %use 11/14/26 as launch date ephemerides
        %Earth
        m1=5.97219e24; %kg
        x1=9.265568400303160E+10; %m
        y1=1.154504238774003E+11; %m
        z1=-7.413929998397827E+6; %m
        vx1=-2.372704237275350E+4; %m/s
        vy1=1.853055826394773E+4;
        vz1=-1.806665656604700E+0;
        
        %Mars
        m2=6.4171e23; %kg
        x2=-6.854726215799426E+10;
        y2=2.304624765567116E+11;
        z2=6.510417970610857E+9;
        vx2=-2.230674236863950E+4;
        vy2=-4.849166862197805E+3;
        vz2=4.453309903383578E+2;
        
        %Sun
        m3=1.989e30; %kg
        x3=0;
        y3=0;
        z3=0;
        vx3=0;
        vy3=0;
        vz3=0;
        
        %Spacecraft (couldn't get it to work so currently set to Earth
        %ephemerides)
        m4=5000; %kg
        x4=(9.265568400303160E+10) + 2000000; %m
        y4=1.154504238774003E+11; %m
        z4=-7.413929998397827E+6; %m
        vx4=(-2.372704237275350E+4)-2e3; %m/s
        vy4=(1.853055826394773E+4)+2e3;%(2.5+0.7+0.6+0.9+0.2+0.3+0.9+4.1);
        vz4=(-1.806665656604700E+0)+2e3;
        
        endtime=4*365.25636*86400;
        i=1e5;
        
        
    case 3 %earth/mars/sun/jupiter
        %Earth
        m1=5.97219e24; %kg
        x1=1.358519213625654E+11; %m
        y1=6.135912462718222E+10; %m
        z1=-3.183658747591078E+6; %m
        vx1=-1.275541057982471E+4; %m/s
        vy1=2.702943044653523E+4;
        vz1=-1.216464121185368E-1;
        
        %Mars
        m2=6.4171e23; %kg
        x2=2.071691995444704E+11;
        y2=-1.411888216421376E+10;
        z2=-5.379329369868487E+9;
        vx2=2.573048525403522E+3;
        vy2=2.624363275060584E+4;
        vz2=4.867769628700440E+2;
        
        %Sun
        m3=1.989e30; %kg
        x3=0;
        y3=0;
        z3=0;
        vx3=0;
        vy3=0;
        vz3=0;
        
        %Jupiter
        m4=1898.13e24; %kg
        x4=-3.941040234284360E+11;
        y4=-7.000842112606339E+11;
        z4=1.172597970064026E+10;
        vx4=1.123915193788393E+4;
        vy4=-5.799311513986108E+3;
        vz4= -2.273837956723392E+2;
        
        endtime=15*365.25*86400;
        i=1e5;
end

% time interval

dt=endtime/i;

X1=zeros(1,i);
Y1=zeros(1,i);
Z1=zeros(1,i);
X2=zeros(1,i);
Y2=zeros(1,i);
Z2=zeros(1,i);
X3=zeros(1,i);
Y3=zeros(1,i);
Z3=zeros(1,i);
X4=zeros(1,i);
Y4=zeros(1,i);
Z4=zeros(1,i);

VX1=zeros(1,i);
VY1=zeros(1,i);
VZ1=zeros(1,i);
VX2=zeros(1,i);
VY2=zeros(1,i);
VZ2=zeros(1,i);
VX3=zeros(1,i);
VY3=zeros(1,i);
VZ3=zeros(1,i);
VX4=zeros(1,i);
VY4=zeros(1,i);
VZ4=zeros(1,i);

% assignment of initial condition to the first elements of parametrics

X1(1)=x1;
Y1(1)=y1;
Z1(1)=z1;
X2(1)=x2;
Y2(1)=y2;
Z2(1)=z2;
X3(1)=x3;
Y3(1)=y3;
Z3(1)=z3;
X4(1)=x4;
Y4(1)=y4;
Z4(1)=z4;

VX1(1)=vx1;
VY1(1)=vy1;
VZ1(1)=vz1;
VX2(1)=vx2;
VY2(1)=vy2;
VZ2(1)=vz2;
VX3(1)=vx3;
VY3(1)=vy3;
VZ3(1)=vz3;
VX4(1)=vx4;
VY4(1)=vy4;
VZ4(1)=vz4;

% computation loop for positions and velocities of the three bodies as functions of t

c=1;

while c<i
    % calculation of acceleration for interval
    
    % distances between particles, dist is + when x1<x2<x3
    %Case 3: x3<x1<x2<x4
    
    dx12=X2(c)-X1(c);
    dx13=X3(c)-X1(c);
    dx23=X3(c)-X2(c);
    dx14 = X4(c) - X1(c);
    dx24 = X4(c) - X2(c);
    dx34 = X3(c) - X4(c);%%%%%
    
    dy12=Y2(c)-Y1(c);
    dy13=Y3(c)-Y1(c);
    dy23=Y3(c)-Y2(c);
    dy14 = Y4(c) - Y1(c);
    dy24 = Y4(c) - Y2(c);
    dy34 = Y3(c) - Y4(c);%%%%%%
    
    dz12=Z2(c)-Z1(c);
    dz13=Z3(c)-Z1(c);
    dz23=Z3(c)-Z2(c);
    dz14 = Z4(c) - Z1(c);
    dz24 = Z4(c) - Z2(c);
    dz34 = Z3(c) - Z4(c);%%%%%%
    
    r12=[dx12 dy12 dz12];
    r13=[dx13 dy13 dz13];
    r23=[dx23 dy23 dz23];
    r14 = [dx14 dy14 dz14];
    r24 = [dx24 dy24 dz24];
    r34 = [dx34 dy34 dz34];
    
    rm12= norm(r12);%(dx12^2+dy12^2+dz12^2)^(1/2);
    rm13= norm(r13);%(dx13^2+dy13^2+dz13^2)^(1/2);
    rm23= norm(r23);%(dx23^2+dy23^2+dz23^2)^(1/2);
    rm14 = norm(r14);
    rm24 = norm(r24);
    rm34 = norm(r34);
    
    rh12=r12/rm12;
    rh13=r13/rm13;
    rh23=r23/rm23;
    rh14 = r14/rm14;
    rh24 = r24/rm24;
    rh34 = r34/rm34;
    
    % acceleration vectors;
    
   if cc==3
       a1=G*((m2/rm12^2)*rh12+(m3/rm13^2)*rh13+(m4/rm14^2)*rh14);
       a2=G*((-m1/rm12^2)*rh12+(m3/rm23^2)*rh23+(m4/rm24^2)*rh24);
       a3=[0 0 0];%G*((-m1/rm13^2)*rh13-(m2/rm23^2)*rh23);
       a4 = G*((m3/rm34^2)*rh34);%G*((-m1/rm14^2)*r14+(m2/rm24^2)*rh24+(-m3/rm34^2)*rh34);
   elseif cc==2
       a1=G*((m2/rm12^2)*rh12+(m3/rm13^2)*rh13+(m4/rm14^2)*rh14);
       a2=G*((-m1/rm12^2)*rh12+(m3/rm23^2)*rh23+(m4/rm24^2)*rh24);
       a3=[0 0 0];%G*((-m1/rm13^2)*rh13-(m2/rm23^2)*rh23);
       a4 = G*(-(m1/rm14^2)*rh14);%G*((m3/rm34^2)*rh34+(m2/rm24^2)*rh24);%G*((-m1/rm14^2)*r14+(m2/rm24^2)*rh24+(-m3/rm34^2)*rh34);
   %issue above with a4; adding in Earth makes it go haywire
   elseif cc==1
              a1=G*((m2/rm12^2)*rh12+(m3/rm13^2)*rh13+(m4/rm14^2)*rh14);
       a2=G*((-m1/rm12^2)*rh12+(m3/rm23^2)*rh23+(m4/rm24^2)*rh24);
       a3=[0 0 0];%G*((-m1/rm13^2)*rh13-(m2/rm23^2)*rh23);
       a4 = G*((m2/rm24^2)*rh24) + rh23*(norm(VX2(c)+VY2(c)+VZ2(c))^2)/rm23;
   else
       a1=G*((m2/rm12^2)*rh12+(m3/rm13^2)*rh13);
       a2=G*((-m1/rm12^2)*rh12+(m3/rm23^2)*rh23);
       a3=[0 0 0];%G*((-m1/rm13^2)*rh13-(m2/rm23^2)*rh23);
       a4 = G*((m3/rm34^2)*rh34);%G*((-m1/rm14^2)*rh14+(m2/rm24^2)*rh24);%G*((-m1/rm14^2)*r14+(m2/rm24^2)*rh24+(-m3/rm34^2)*rh34);
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculation and assignment of new positions
    
    X1(c+1)=X1(c)+VX1(c)*dt+(1/2)*a1(1)*dt^2;
    X2(c+1)=X2(c)+VX2(c)*dt+(1/2)*a2(1)*dt^2;
    X3(c+1)=X3(c)+VX3(c)*dt+(1/2)*a3(1)*dt^2;
    X4(c+1) = X4(c)+VX4(c)*dt+(1/2)*a4(1)*dt^2;
    
    Y1(c+1)=Y1(c)+VY1(c)*dt+(1/2)*a1(2)*dt^2;
    Y2(c+1)=Y2(c)+VY2(c)*dt+(1/2)*a2(2)*dt^2;
    Y3(c+1)=Y3(c)+VY3(c)*dt+(1/2)*a3(2)*dt^2;
    Y4(c+1) = Y4(c)+VY4(c)*dt+(1/2)*a4(2)*dt^2;
    
    Z1(c+1)=Z1(c)+VZ1(c)*dt+(1/2)*a1(3)*dt^2;
    Z2(c+1)=Z2(c)+VZ2(c)*dt+(1/2)*a2(3)*dt^2;
    Z3(c+1)=Z3(c)+VZ3(c)*dt+(1/2)*a3(3)*dt^2;
    Z4(c+1) = Z4(c)+VZ4(c)*dt+(1/2)*a4(3)*dt^2;
    
    % calculation and assignment of velocities to be used for next interval
    
    VX1(c+1)=VX1(c)+a1(1)*dt;
    VX2(c+1)=VX2(c)+a2(1)*dt;
    VX3(c+1)=VX3(c)+a3(1)*dt;
    VX4(c+1) = VX4(c)+a4(1)*dt;
    
    VY1(c+1)=VY1(c)+a1(2)*dt;
    VY2(c+1)=VY2(c)+a2(2)*dt;
    VY3(c+1)=VY3(c)+a3(2)*dt;
    VY4(c+1) = VY4(c)+a4(2)*dt;
    
    VZ1(c+1)=VZ1(c)+a1(3)*dt;
    VZ2(c+1)=VZ2(c)+a2(3)*dt;
    VZ3(c+1)=VZ3(c)+a3(3)*dt;
    VZ4(c+1) = VZ4(c)+a4(3)*dt;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % increase counter by one before beginning the next loop
    
    c=c+1;
end

% subplot(1,2,1);
posit=plot3(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,x1,y1,z1,'ob',x2,y2,z2,'or',x3,y3,z3,'oy',x4,y4,z4,'ok');
title('Positions as a function of t');
xlabel('x');
ylabel('y');
zlabel('z');
daspect([1 1 1]);
set(gca,'PlotBoxAspectRatio',[1 1 1]);
grid on;

if cc == 3
    axis([-10e11 10e11 -10e11 10e11 -2e11 2e11]);
else
    axis([-2.5e11 2.5e11 -2.25e11 2.5e11 -1e11 1e11]);
end
view(-37.5,30);

% subplot(1,2,2);
% velo=plot3(VX1,VY1,VZ1,VX2,VY2,VZ2,VX3,VY3,VZ3,vx1,vy1,vz1,'o',vx2,vy2,vz2,'o',vx3,vy3,vz3,'o');
% legend('m1','m2','m3','m1 current','m2 current','m3 current','Location','best');
% title('Velocities as a function of t');
% xlabel('vx');
% ylabel('vy');
% zlabel('vz');
% daspect([1 1 1]);
% set(gca,'PlotBoxAspectRatio',[1 1 1]);
% grid on;
% set(gcf,'Position',[0 356 1440 900]);

for rd=1:i/300:i
    
    rd=floor(rd);
    
    % reset all values for position and velocity, then redraw
    
    set(posit(1),'XData',X1(1:rd));
    set(posit(1),'YData',Y1(1:rd));
    set(posit(1),'ZData',Z1(1:rd));
    set(posit(2),'XData',X2(1:rd));
    set(posit(2),'YData',Y2(1:rd));
    set(posit(2),'ZData',Z2(1:rd));
    set(posit(3),'XData',X3(1:rd));
    set(posit(3),'YData',Y3(1:rd));
    set(posit(3),'ZData',Z3(1:rd));
    set(posit(4),'XData',X4(1:rd));
    set(posit(4),'YData',Y4(1:rd));
    set(posit(4),'ZData',Z4(1:rd));
    
    set(posit(5),'XData',X1(rd));
    set(posit(5),'YData',Y1(rd));
    set(posit(5),'ZData',Z1(rd));
    set(posit(6),'XData',X2(rd));
    set(posit(6),'YData',Y2(rd));
    set(posit(6),'ZData',Z2(rd));
    set(posit(7),'XData',X3(rd));
    set(posit(7),'YData',Y3(rd));
    set(posit(7),'ZData',Z3(rd));
    set(posit(8),'XData',X4(rd));
    set(posit(8),'YData',Y4(rd));
    set(posit(8),'ZData',Z4(rd));
    
%     set(velo(1),'XData',VX1(1:rd));
%     set(velo(1),'YData',VY1(1:rd));
%     set(velo(1),'ZData',VZ1(1:rd));
%     set(velo(2),'XData',VX2(1:rd));
%     set(velo(2),'YData',VY2(1:rd));
%     set(velo(2),'ZData',VZ2(1:rd));
%     set(velo(3),'XData',VX3(1:rd));
%     set(velo(3),'YData',VY3(1:rd));
%     set(velo(3),'ZData',VZ3(1:rd));
%     set(velo(4),'XData',VX1(rd));
%     set(velo(4),'YData',VY1(rd));
%     set(velo(4),'ZData',VZ1(rd));
%     set(velo(5),'XData',VX2(rd));
%     set(velo(5),'YData',VY2(rd));
%     set(velo(5),'ZData',VZ2(rd));
%     set(velo(6),'XData',VX3(rd));
%     set(velo(6),'YData',VY3(rd));
%     set(velo(6),'ZData',VZ3(rd));
    
    drawnow
end