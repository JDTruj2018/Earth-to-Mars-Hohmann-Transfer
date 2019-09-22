%% Author: Jered Dominguez-Trujillo
% Date: March 6, 2017
clear vars; close all;
set(gcf,'position',get(0,'screensize'))

%% Declare Variables
AU = 149597870700;
G = 6.674275935628303e-11;
MSolarSystem = 1.999*10^30;
RMars = 3.3899*10^6;
EarthYear = 365.24237; EarthDay = 86400;
dt = 1; inx = 1;
StartYear = 2010; years = 30;
length = years * EarthYear * EarthDay;
dates = datetime(StartYear,01,01):caldays(1):datetime(StartYear+years+1,01,01);
t = zeros(1, ceil(length/dt));
t(inx) = 0;

axEarth = zeros(1, ceil(length/dt)); ayEarth = zeros(1, ceil(length/dt));
azEarth = zeros(1, ceil(length/dt));
vxEarth = zeros(1, ceil(length/dt)); vyEarth = zeros(1, ceil(length/dt));
vzEarth = zeros(1, ceil(length/dt));
xEarth = zeros(1, ceil(length/dt)); yEarth = zeros(1, ceil(length/dt));
zEarth = zeros(1, ceil(length/dt));
rEarth = zeros(1, ceil(length/dt));
alphaEarth = zeros(1, ceil(length/dt)); betaEarth = zeros(1, ceil(length/dt));
gammaEarth = zeros(1, ceil(length/dt));

axMars = zeros(1, ceil(length/dt)); ayMars = zeros(1, ceil(length/dt));
azMars = zeros(1, ceil(length/dt));
vxMars = zeros(1, ceil(length/dt)); vyMars = zeros(1, ceil(length/dt));
vzMars = zeros(1, ceil(length/dt));
xMars = zeros(1, ceil(length/dt)); yMars = zeros(1, ceil(length/dt));
zMars = zeros(1, ceil(length/dt));
rMars = zeros(1, ceil(length/dt));
alphaMars = zeros(1, ceil(length/dt)); betaMars = zeros(1, ceil(length/dt));
gammaMars = zeros(1, ceil(length/dt));

axRocket = zeros(1, ceil(length/dt)); ayRocket = zeros(1, ceil(length/dt));
azRocket = zeros(1, ceil(length/dt));
vxRocket = zeros(1, ceil(length/dt)); vyRocket = zeros(1, ceil(length/dt));
vzRocket = zeros(1, ceil(length/dt));
xRocket = zeros(1, ceil(length/dt)); yRocket = zeros(1, ceil(length/dt));
zRocket = zeros(1, ceil(length/dt));
rRocket = zeros(1, ceil(length/dt));
alphaRocket = zeros(1, ceil(length/dt)); betaRocket = zeros(1, ceil(length/dt));
gammaRocket = zeros(1, ceil(length/dt));

%% Data Gathered for Earth and Mars for January 1, 2010 from JPL
axEarth(inx) = 0; ayEarth(inx) = 0; azEarth(inx) = 0;
vxEarth(inx) = -2.978405751621624e+04; vyEarth(inx) = -5.451137243323289e+03;
vzEarth(inx) = 1.551580077431058;
xEarth(inx) = -2.689245210784379e+10; yEarth(inx) = 1.451618583042868e+11;
zEarth(inx) = -2.608728054337204e+06;
rEarth(inx) = sqrt(xEarth(inx)^2+yEarth(inx)^2+zEarth(inx)^2);
alphaEarth(inx) = acos(xEarth(inx)/rEarth(inx)); betaEarth(inx) = acos(yEarth(inx)/rEarth(inx));
gammaEarth(inx) = acos(zEarth(inx)/rEarth(inx));                                         axMars(inx) = 0; ayMars(inx) = 0; azMars(inx) = 0;
vxMars(inx) = -2.074332689615267e+04; vyMars(inx) = -8.817501559049284e+03;
vzMars(inx) = 3.248028553828797e+02;
xMars(inx) = -1.097178350768191e+11; yMars(inx) = 2.179958733988708e+11;
zMars(inx) =  7.239687840140940e+09;
rMars(inx) = sqrt(xMars(inx)^2+yMars(inx)^2+zMars(inx)^2);
alphaMars(inx) = acos(xMars(inx)/rMars(inx)); betaMars(inx) = acos(yMars(inx)/rMars(inx));
gammaMars(inx) = acos(zMars(inx)/rMars(inx));
distanceEarthtoMars = zeros(1, ceil(length/dt));
distanceEarthtoMars(inx) = sqrt((xMars(inx)-xEarth(inx))^2+(yMars(inx)-yEarth(inx))^2+(zMars(inx)-zEarth(inx))^2);                                         

%% Calculate Earth and Mars Trajectories
while t(inx) < length
    inx = inx + 1;
    alphaEarth(inx)=acos(xEarth(inx-1)/rEarth(inx-1));
    betaEarth(inx)=acos(yEarth(inx-1)/rEarth(inx-1));
    gammaEarth(inx) = acos(zEarth(inx-1)/rEarth(inx-1));
    rEarth(inx) = sqrt(xEarth(inx-1)^2+yEarth(inx-1)^2+zEarth(inx-1)^2);
    axEarth(inx) = -(G*MSolarSystem)/(rEarth(inx)^2)*cos(alphaEarth(inx));
    ayEarth(inx) = -(G*MSolarSystem)/(rEarth(inx)^2)*cos(betaEarth(inx));
    azEarth(inx) = -(G*MSolarSystem)/(rEarth(inx)^2)*cos(gammaEarth(inx));
    vxEarth(inx)=vxEarth(inx-1)+axEarth(inx)*dt; vyEarth(inx)=vyEarth(inx-1)+ayEarth(inx)*dt;
    vzEarth(inx)=vzEarth(inx-1)+azEarth(inx)*dt;
    xEarth(inx)=xEarth(inx-1)+vxEarth(inx)*dt; yEarth(inx)=yEarth(inx-1)+vyEarth(inx)*dt;
    zEarth(inx)=zEarth(inx-1)+vzEarth(inx)*dt;
    alphaMars(inx)=acos(xMars(inx-1)/rMars(inx-1));betaMars(inx)=acos(yMars(inx-1)/rMars(inx-1));
    gammaMars(inx) = acos(zMars(inx-1)/rMars(inx-1));
    rMars(inx) = sqrt(xMars(inx-1)^2+yMars(inx-1)^2+zMars(inx-1)^2);
    axMars(inx) = -(G*MSolarSystem)/(rMars(inx)^2)*cos(alphaMars(inx));
    ayMars(inx) = -(G*MSolarSystem)/(rMars(inx)^2)*cos(betaMars(inx));
    azMars(inx) = -(G*MSolarSystem)/(rMars(inx)^2)*cos(gammaMars(inx));
    vxMars(inx)=vxMars(inx-1)+axMars(inx)*dt; vyMars(inx)=vyMars(inx-1)+ayMars(inx)*dt;
    vzMars(inx)=vzMars(inx-1)+azMars(inx)*dt;
    xMars(inx)=xMars(inx-1)+vxMars(inx)*dt; yMars(inx) = yMars(inx-1) + vyMars(inx) * dt;
    zMars(inx) = zMars(inx-1) + vzMars(inx) * dt;
    distanceEarthtoMars(inx) = sqrt((xMars(inx)-xEarth(inx))^2+(yMars(inx)-yEarth(inx))^2+(zMars(inx)-zEarth(inx))^2);

    t(inx) = t(inx-1) + dt;
end

%% Find Time and Distance of Closest Approach
closestApproach = find(distanceEarthtoMars==min(distanceEarthtoMars))*dt/EarthDay/EarthYear;
ApproachYear = StartYear + floor(closestApproach);
ApproachDay = floor(rem(closestApproach, 1) * EarthYear);
datesApproach = datetime(ApproachYear,01,01):caldays(1):datetime(ApproachYear,12,31);
fprintf('Mars Closest Approach to Earth between %s and %s occurs on %s at a distance of %.4f AU\n',...
    num2str(StartYear), num2str(StartYear + years), datestr(datesApproach(ApproachDay)), ...
    min(distanceEarthtoMars)/AU);

%% Calculating Transit Time and Target Distance if Launched at Each Time Step
alphaEarth = round(atan2(real(yEarth), real(xEarth)), 3);
alphaMars = round(atan2(real(yMars), real(xMars)), 3);
TransitTime = zeros(1, ceil(inx)); TargetDistance = zeros(1, ceil(inx));
for ii = 1:0.5*inx
    if alphaEarth(ii) > 0
        value = alphaEarth(ii) - round(pi, 3);
    else
        value = alphaEarth(ii) + round(pi, 3);
    end
    MarsOppositeInx = find(abs(alphaMars - value)<0.0001);
    if isempty(MarsOppositeInx) == 0
        d = sqrt((xEarth(ii) - xMars(MarsOppositeInx(1)))^2+(yEarth(ii) - yMars(MarsOppositeInx(1)))^2+(zEarth(ii) - zMars(MarsOppositeInx(1)))^2);
        TransitTime(ii) = ((sqrt((((d/AU))/2)^3))/2)*EarthYear;
        TimeinSeconds=TransitTime(ii)*EarthDay; TimeSteps=real(ceil(TimeinSeconds/dt));
        TargetDistance(ii) = sqrt(((xMars(MarsOppositeInx(1)) - xMars(ii+TimeSteps))^2+...
            (yMars(MarsOppositeInx(1)) - yMars(ii+TimeSteps))^2+...
            (zMars(MarsOppositeInx(1)) - zMars(ii+TimeSteps))^2));
    else
        TransitTime(ii) = TransitTime(ii - 1); TargetDistance(ii) = TargetDistance(ii - 1);
    end
end

%% Print Appropriate Launch Dates
LaunchDates = find(TargetDistance(1:floor(0.5*inx))/AU < 5*RMars/AU);
datecount = size(LaunchDates);
for ii = 1:datecount(2)
    Val = LaunchDates(ii)*dt/EarthDay/EarthYear;
    LaunchHour = round(rem(LaunchDates(ii)/24, 1) * 24, 0);
    LaunchYear = StartYear + floor(Val); LaunchDay = floor(rem(Val, 1) * EarthYear);
    datesLaunch = datetime(LaunchYear,01,01):caldays(1):datetime(LaunchYear,12,31);
    fprintf('Launch to Mars: %s. Closest Approach to Mars: %.4f Mars Radii.  Travel Time: %.4f Days\n',...
        datestr(datesLaunch(LaunchDay)), TargetDistance(LaunchDates(ii))/RMars, ...
        TransitTime(LaunchDates(ii)));
end

%% Calculate Rocket Start and End Positions
LaunchofInterest = LaunchDates(end);  % 2018 Launch
if alphaEarth(LaunchofInterest) > 0
    val = alphaEarth(LaunchofInterest) - round(pi, 3);
else
    val = alphaEarth(LaunchofInterest) + round(pi, 3);
end
MarsOppositeInx = find(abs(alphaMars - val)<0.0001);
dEtoM = sqrt(((xEarth(LaunchofInterest)-xMars(MarsOppositeInx(1)))^2)+...
    ((yEarth(LaunchofInterest)-yMars(MarsOppositeInx(1)))^2)+...
    ((zEarth(LaunchofInterest)-zMars(MarsOppositeInx(1)))^2));
Vp = sqrt(G*MSolarSystem*((2/rEarth(LaunchofInterest))-(1/(dEtoM/2))));

%% Initial Rocket Parameters Calculated to Complete Transfer Orbit
inx = 1; t = zeros(1, ceil(length/dt)); t(inx) = 0; rad = 0 ; vtest = 100000;              axRocket(inx) = 0; ayRocket(inx) = 0; azRocket(inx) = 0;
while vtest > Vp
    vxRocket(inx) = ((-1*vxEarth(LaunchofInterest))/abs(vxEarth(LaunchofInterest)))*...
        (-abs(vxEarth(LaunchofInterest)) + (Vp-abs(vxEarth(LaunchofInterest)))*(1-cos(rad)));
    vyRocket(inx) = ((-1*vyEarth(LaunchofInterest))/abs(vyEarth(LaunchofInterest)))*...
        (abs(vyEarth(LaunchofInterest))+(Vp+abs(vyEarth(LaunchofInterest)))*cos(pi-rad));
    vzRocket(inx) = ((vzEarth(LaunchofInterest))/abs(vzEarth(LaunchofInterest)))*...
        (abs(vzEarth(LaunchofInterest))-(Vp-abs(vzEarth(LaunchofInterest)))*cos(rad));
    vtest = sqrt(vxRocket(inx)^2 + vyRocket(inx)^2 + vzRocket(inx)^2); rad = rad + 0.00001;
end
xRocket(inx) = xEarth(LaunchofInterest); yRocket(inx) = yEarth(LaunchofInterest);
zRocket(inx) =  zEarth(LaunchofInterest);
rRocket(inx) = sqrt(xRocket(inx)^2+yRocket(inx)^2+zRocket(inx)^2);
alphaRocket(inx)=acos(xRocket(inx)/rRocket(inx));betaRocket(inx)=acos(yRocket(inx)/rRocket(inx));
gammaRocket(inx) = acos(zRocket(inx)/rRocket(inx));

%% Calculate Rocket Trajectory
while t(inx) < length
    inx = inx + 1;
    alphaRocket(inx) = acos(xRocket(inx-1)/rRocket(inx-1));
    betaRocket(inx) = acos(yRocket(inx-1)/rRocket(inx-1));
    gammaRocket(inx) = acos(zRocket(inx-1)/rRocket(inx-1));
    rRocket(inx) = sqrt(xRocket(inx-1)^2+yRocket(inx-1)^2+zRocket(inx-1)^2);
    axRocket(inx) = -(G*MSolarSystem)/(rRocket(inx)^2)*cos(alphaRocket(inx));
    ayRocket(inx) = -(G*MSolarSystem)/(rRocket(inx)^2)*cos(betaRocket(inx));
    azRocket(inx) = -(G*MSolarSystem)/(rRocket(inx)^2)*cos(gammaRocket(inx));
    vxRocket(inx) = vxRocket(inx-1) + axRocket(inx) * dt;
    vyRocket(inx) = vyRocket(inx-1) + ayRocket(inx) * dt;
    vzRocket(inx) = vzRocket(inx-1) + azRocket(inx) * dt;
    xRocket(inx)=xRocket(inx-1)+vxRocket(inx)*dt; yRocket(inx)=yRocket(inx-1)+vyRocket(inx)*dt;
    zRocket(inx) = zRocket(inx-1) + vzRocket(inx) * dt;
    t(inx) = t(inx - 1) + dt;
end

%% Plot Animation from 2010 Until Arrival at Mars
flag =false;
for ii = 1:LaunchofInterest + (real(max(TransitTime))*EarthDay/dt)
     if rem(t(ii), EarthDay) == 0 && t(ii) > 0
         clf; hold on; view(3)
         plot3(real(xEarth(ii))/AU, real(yEarth(ii))/AU, real(zEarth(ii))/AU, 'bo', ...
             'MarkerFaceColor','b')
         plot3(real(xMars(ii))/AU, real(yMars(ii))/AU, real(zMars(ii))/AU, 'ro', ...
             'MarkerFaceColor','r')
         plot3(0, 0, 0, 'go', 'MarkerFaceColor','g')
         if ii >= LaunchofInterest
             plot3(real(xRocket(ii-LaunchofInterest))/AU, real(yRocket(ii-LaunchofInterest))/AU,...
                 real(zRocket(ii-LaunchofInterest))/AU, 'ko', 'MarkerFaceColor', 'k')
             plot3(real(xRocket)/AU, real(yRocket)/AU, real(zRocket)/AU, 'k', 'linewidth', 0.25)
             flag = true;
         end
         plot3(real(xEarth/AU), real(yEarth/AU), real(zEarth/AU))
         plot3(real(xMars/AU), real(yMars/AU), real(zMars/AU))
         xlim([-2, 2]); ylim([-2, 2]); zlim([-1.5, 1.5])
         xlabel('AU'); ylabel('AU'); zlabel('AU')
         title('Orbital Simulation')
         text(0, 0.2, 0.2, datestr(dates(ceil(t(ii)/EarthDay))));
         if flag == true
            legend({'Earth', 'Mars', 'Solar System Barycenter', 'Rocket'})
         else
             legend({'Earth', 'Mars', 'Solar System Barycenter'})
         end
         grid on
         print(['Simulation\Time' num2str(ii)], '-dpng')
         pause(0.01)
     end
end

%% Plot and Save Figures
clf; close all;
set(gcf,'position',get(0,'screensize'))
figure(1); hold on; view(3); 
plot3(xEarth/AU, yEarth/AU, zEarth/AU, 'b');
plot3(xMars/AU, yMars/AU, zMars/AU, 'r'); plot3(0, 0, 0, 'go', 'MarkerFaceColor','g')
plot3(xRocket/AU, yRocket/AU, zRocket/AU, 'k');
plot3(xEarth(LaunchofInterest)/AU, yEarth(LaunchofInterest)/AU, zEarth(LaunchofInterest)/AU, 'bo', 'MarkerFaceColor', 'b');
plot3(xMars(MarsOppositeInx(1))/AU, yMars(MarsOppositeInx(1))/AU, zMars(MarsOppositeInx(1))/AU, 'ro', 'MarkerFaceColor', 'r');
legend({'Earth', 'Mars', 'Solar System Barycenter', 'Rocket'})
title('3D Transfer Orbit');
xlabel('X Distance [AU]'); ylabel('Y Distance [AU]'); zlabel('Z Distance [AU'); grid on;
print('3D Transfer Orbit', '-dpng');

figure(2); set(gcf,'position',get(0,'screensize')); hold on; plot(xEarth/AU, yEarth/AU, 'b'); plot(xMars/AU, yMars/AU, 'r'); 
plot(0, 0, 'go', 'MarkerFaceColor', 'g'); plot(xRocket/AU, yRocket/AU, 'k');
plot(xEarth(LaunchofInterest)/AU, yEarth(LaunchofInterest)/AU, 'bo', 'MarkerFaceColor', 'b')
plot(xMars(MarsOppositeInx(1))/AU, yMars(MarsOppositeInx(1))/AU, 'ro', 'MarkerFaceColor', 'r')
legend({'Earth', 'Mars', 'Solar System Barycenter', 'Rocket'})
title('XY Transfer Orbit');
xlabel('X Distance [AU]'); ylabel('Y Distance [AU]'); grid on;
print('XY Transfer Orbit', '-dpng');

figure(3); set(gcf,'position',get(0,'screensize')); hold on; plot(yEarth/AU, zEarth/AU, 'b'); plot(yMars/AU, zMars/AU, 'r'); 
plot(0, 0, 'go', 'MarkerFaceColor', 'g'); plot(yRocket/AU, zRocket/AU, 'k');
plot(yEarth(LaunchofInterest)/AU, zEarth(LaunchofInterest)/AU, 'bo', 'MarkerFaceColor', 'b')
plot(yMars(MarsOppositeInx(1))/AU, zMars(MarsOppositeInx(1))/AU, 'ro', 'MarkerFaceColor', 'r')
legend({'Earth', 'Mars', 'Solar System Barycenter', 'Rocket'})
title('YZ Transfer Orbit');
xlabel('Y Distance [AU]'); ylabel('Z Distance [AU]'); grid on;
print('YZ Transfer Orbit', '-dpng');

figure(4); set(gcf,'position',get(0,'screensize')); hold on; plot(xEarth/AU, zEarth/AU, 'b'); plot(xMars/AU, zMars/AU, 'r'); 
plot(0, 0, 'go', 'MarkerFaceColor', 'g'); plot(xRocket/AU, zRocket/AU, 'k');
plot(xEarth(LaunchofInterest)/AU, zEarth(LaunchofInterest)/AU, 'bo', 'MarkerFaceColor', 'b')
plot(xMars(MarsOppositeInx(1))/AU, zMars(MarsOppositeInx(1))/AU, 'ro', 'MarkerFaceColor', 'r')
legend({'Earth', 'Mars', 'Solar System Barycenter', 'Rocket'})
title('XZ Transfer Orbit');
xlabel('X Distance [AU]'); ylabel('Z Distance [AU]'); grid on;
print('XZ Transfer Orbit', '-dpng');

figure(5);
plot(t(1:inx)/EarthDay/EarthYear, real(distanceEarthtoMars(1:inx)) / AU);
ylim([0, 3]);
xlabel('Days since January 1, 2010 [Years]'); ylabel('Distance [AU]');
title('Distance between Earth and Mars');
xlim([0 years])
grid on;
print('Distance Earth and Mars', '-dpng');

figure(6);
plot(t(1:floor(0.5*inx))/EarthDay/EarthYear, real(TargetDistance(1:floor(0.5*inx)))/AU);
xlabel('Time since January 1, 2010 [Years]'); ylabel('Distance [AU]');
title('Launch Window');
grid on;
print('Launch Window', '-dpng');

figure(7);
plot(t(1:floor(0.5*inx))/EarthDay/EarthYear, real(TransitTime(1:floor(0.5*inx))));
xlabel('Time since January 1, 2010 [Years]'); ylabel('Transit Time [Days]');
title('Time of Travel to Mars');
grid on;
print('Time of Travel', '-dpng');
