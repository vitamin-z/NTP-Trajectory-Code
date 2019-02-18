m = mobiledev;
t = 1;
while t<51
    time(t) = t;
    ornt(t,1:3) = m.Orientation + 180; %+180 gives all angles >0
    t = t+1;
    pause(.1)
end
figure(1)
plot(time/10,ornt)
xlabel('Time (s)')
ylabel('Angle (deg)')
legend('Azimuth','Pitch','Roll','location','northwest')