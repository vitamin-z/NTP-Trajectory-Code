% Zach Zoloty and Leland Klein
% JPL Ephemerides Non-Equation Test
% Uses only ephemeride data to plot the trajectories of Earth and Mars
% Does not take the gravitational effect of multiple bodies into account

clc, clear;

%Import Ephimeride data and convert to array form
%Columns 2-4 are X,Y,Z and columns 5-7 are U,V,W
%We need X, Y, and Z to plot so those are the state arrays
Mars = table2array(readtable('Mars.txt'));
state_m(:,1) = Mars(:,2);
state_m(:,2) = Mars(:,3);
state_m(:,3) = Mars(:,4);
Earth = table2array(readtable('Earth.txt'));
state_e(:,1) = Earth(:,2);
state_e(:,2) = Earth(:,3);
state_e(:,3) = Earth(:,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Attempt at modeling a spacecraft trajectory between Earth and Mars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot resulting orbits
figure(1)
plot3(state_e(:,1),state_e(:,2),state_e(:,3),'b', 'linewidth', 1);
hold on;
plot3(state_m(:,1),state_m(:,2),state_m(:,3),'r', 'linewidth', 1);
hold on;
plot3(0, 0, 0, 'ro', 'Markersize', 16, 'markerfacecolor', 'y');

axis([-2.5e8 2.5e8 -2.25e8 2.5e8 -0.2e8 0.2e8]);
set(gca, 'fontsize', 14, 'fontweight', 'bold');
xlabel('x');
ylabel('y');
zlabel('z');
grid on;
%hold off;

% Moving "Plot"
figure(2)
hold on;
for j = 1:1:length(state_m)
    plot3(state_e(1:j,1), state_e(1:j,2), state_e(1:j,3), 'b', 'linewidth', 1);
    hold on;
    plot3(state_m(1:j,1), state_m(1:j,2), state_m(1:j,3), 'r', 'linewidth', 1);
    
    hold on;
    plot3(state_e(j,1), state_e(j,2), state_e(j,3), 'go', 'Markersize', 8, 'markerfacecolor', 'b');
    plot3(state_m(j,1), state_m(j,2), state_m(j,3), 'yo', 'Markersize', 6, 'markerfacecolor', 'r');
    
    plot3(0, 0, 0, 'ro', 'Markersize', 16, 'markerfacecolor', 'y');
    
    axis([-2.5e8 2.5e8 -2.25e8 2.5e8 -0.2e8 0.2e8]);
    set(gca, 'fontsize', 14, 'fontweight', 'bold');
    view(100,24);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    grid on;
    hold off;
    getframe(figure(2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%