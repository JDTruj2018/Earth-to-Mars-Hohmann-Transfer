%% Author: Jered Dominguez-Trujillo
% Date: March 6, 2017

% MEarth = 5.97219*10^24; 
% MMars= 6.4185*10^23;
% MMoon = 7.34767309*10^22;

% RSun2Earth = 149597870700; 
% RSun2Mars = 227.9*10^9;
% REearth2Moon = 384.4 * 10^6;

% RSun = 695.7*10^6; 
% REarth = 6.371*10^6; 
% RMoon = 1.737*10^6; 
% RMars = 3.3899*10^6;

% EarthEccentricity = 0.0174;
% EarthPerihelion = RSun2Earth*(1-EarthEccentricity);
% EarthAphelion = RSun2Earth*(1+EarthEccentricity);

% MarsEccentricity = 0.0934;
% MarsPerihelion = RSun2Mars*(1-MarsEccentricity);
% MarsAphelion = RSun2Mars*(1+MarsEccentricity);

% VEarthPerihelion = sqrt(G*MSun*((2/EarthPerihelion)-(1/RSun2Earth)));
% VEarthAphelion = sqrt(G*MSun*((2/EarthAphelion)-(1/RSun2Earth)));

% VMarsPerihelion = sqrt(G*MSun*((2/MarsPerihelion)-(1/RSun2Mars)));
% VMarsAphelion = sqrt(G*MSun*((2/MarsAphelion)-(1/RSun2Mars)));

clear all; close all;
set(gcf,'position',get(0,'screensize'))

AU = 149597870700;
G = 6.674275935628303e-11;

MSolarSystem = 1.999*10^30; 

EarthYear = 365.24237; EarthDay = 86400;
dt = 3600; inx = 1;
StartYear = 2010; years = 2;
length = years * EarthYear * EarthDay;
dates = datetime(StartYear,01,01):caldays(1):datetime(StartYear+years+1,01,01);

t = zeros(1, ceil(length/dt));
t(inx) = 0;

xEarth = zeros(1, ceil(length/dt));
yEarth = zeros(1, ceil(length/dt));
zEarth = zeros(1, ceil(length/dt));
rEarth = zeros(1, ceil(length/dt));

axEarth = 0; ayEarth = 0; azEarth = 0;
vxEarth = -2.978405751621624e+04; vyEarth = -5.451137243323289e+03; vzEarth = 1.551580077431058;
xEarth(inx) = -2.689245210784379e+10; yEarth(inx) = 1.451618583042868e+11; zEarth(inx) = -2.608728054337204e+06;
rEarth(inx) = sqrt(xEarth(inx)^2+yEarth(inx)^2+zEarth(inx)^2); 
alphaEarth = acos(xEarth(inx)/rEarth(inx));
betaEarth = acos(yEarth(inx)/rEarth(inx));
gammaEarth = acos(zEarth(inx)/rEarth(inx));


xMars = zeros(1, ceil(length/dt));
yMars = zeros(1, ceil(length/dt));
zMars = zeros(1, ceil(length/dt));
rMars = zeros(1, ceil(length/dt));

axMars = 0; ayMars = 0; azMars = 0;
vxMars = -2.074332689615267e+04; vyMars = -8.817501559049284e+03; vzMars = 3.248028553828797e+02;
xMars(inx) = -1.097178350768191e+11; yMars(inx) = 2.179958733988708e+11; zMars(inx) =  7.239687840140940e+09;
rMars(inx) = sqrt(xMars(inx)^2+yMars(inx)^2+zMars(inx)^2); 
alphaMars = acos(xMars(inx)/rMars(inx));
betaMars = acos(yMars(inx)/rMars(inx));
gammaMars = acos(zMars(inx)/rMars(inx));

distanceEarthtoMars = zeros(1, ceil(length/dt));
distanceEarthtoMars(inx) = sqrt((xMars(inx)-xEarth(inx))^2+(yMars(inx)-yEarth(inx))^2+(zMars(inx)-zEarth(inx))^2);

while t(inx) < length
    inx = inx + 1;
    alphaEarth = acos(xEarth(inx-1)/rEarth(inx-1));
    betaEarth = acos(yEarth(inx-1)/rEarth(inx-1));
    gammaEarth = acos(zEarth(inx-1)/rEarth(inx-1));
    rEarth(inx) = sqrt(xEarth(inx-1)^2+yEarth(inx-1)^2+zEarth(inx-1)^2);
    axEarth = -(G*MSolarSystem)/(rEarth(inx)^2)*cos(alphaEarth);
    ayEarth = -(G*MSolarSystem)/(rEarth(inx)^2)*cos(betaEarth);
    azEarth = -(G*MSolarSystem)/(rEarth(inx)^2)*cos(gammaEarth);
    vxEarth = vxEarth + axEarth * dt;
    vyEarth = vyEarth + ayEarth * dt;
    vzEarth = vzEarth + azEarth * dt;
    xEarth(inx) = xEarth(inx-1) + vxEarth * dt;
    yEarth(inx) = yEarth(inx-1) + vyEarth * dt;
    zEarth(inx) = zEarth(inx-1) + vzEarth * dt;
    
    alphaMars = acos(xMars(inx-1)/rMars(inx-1));
    betaMars = acos(yMars(inx-1)/rMars(inx-1));
    gammaMars = acos(zMars(inx-1)/rMars(inx-1));
    rMars(inx) = sqrt(xMars(inx-1)^2+yMars(inx-1)^2+zMars(inx-1)^2);
    axMars = -(G*MSolarSystem)/(rMars(inx)^2)*cos(alphaMars);
    ayMars = -(G*MSolarSystem)/(rMars(inx)^2)*cos(betaMars);
    azMars = -(G*MSolarSystem)/(rMars(inx)^2)*cos(gammaMars);
    vxMars = vxMars + axMars * dt;
    vyMars = vyMars + ayMars * dt;
    vzMars = vzMars + azMars * dt;
    xMars(inx) = xMars(inx-1) + vxMars * dt;
    yMars(inx) = yMars(inx-1) + vyMars * dt;
    zMars(inx) = zMars(inx-1) + vzMars * dt;
    
    distanceEarthtoMars(inx) = sqrt((xMars(inx)-xEarth(inx))^2+(yMars(inx)-yEarth(inx))^2+(zMars(inx)-zEarth(inx))^2);
    
    t(inx) = t(inx-1) + dt;
end

find(distanceEarthtoMars==min(distanceEarthtoMars))*dt/86400/EarthYear

for ii = 1:inx
     if rem(t(ii), EarthDay) == 0 && t(ii) > 0
         clf; hold on; view(3)
         plot3(real(xEarth(ii))/AU, real(yEarth(ii))/AU, real(zEarth(ii))/AU, 'bo', ...
             'MarkerFaceColor','b')
         plot3(real(xMars(ii))/AU, real(yMars(ii))/AU, real(zMars(ii))/AU, 'ro', ...
             'MarkerFaceColor','r')
         plot3(0, 0, 0, 'ko', 'MarkerFaceColor','k')
         plot3(real(xEarth/AU), real(yEarth/AU), real(zEarth/AU))
         plot3(real(xMars/AU), real(yMars/AU), real(zMars/AU))
         xlim([-2, 2]); ylim([-2, 2]); zlim([-0.1, 0.1])
         xlabel('AU'); ylabel('AU'); zlabel('AU')
         title('Orbital Simulation')
         text(0, 0.2, 0.2, datestr(dates(ceil(t(ii)/EarthDay))));
         legend({'Earth', 'Mars', 'Solar System Barycenter'})
         grid on
         pause(0.01)
     end
end

clf;
plot(t(1:inx) / EarthDay, real(distanceEarthtoMars(1:inx)) / AU)
ylim([0, 3])
xlabel('Days'); ylabel('AU')
title('Distance between Earth and Mars');
grid on