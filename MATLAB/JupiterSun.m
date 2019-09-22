clear all; close all;
AU = 149597870700; G = 6.674275935628303e-11;
day = 86400;
dt = 0.1;

inx = 1;

MSun = 1.988544*10^30;
MJupiter = 1898.13 * 10^24;
MTotal = MSun + MJupiter;

t(inx) = 0;

xSun(inx) = -3.747147775505800E-03 * AU;
ySun(inx) = 2.926587035303590E-03 * AU;
zSun(inx) = 4.446807296978120E-06 * AU;

vxSun(inx) = -2.990258994506073E-06 * AU / day;
vySun(inx) = -5.530878190433255E-06 * AU / day;
vzSun(inx) = 6.981801439481428E-08 * AU / day;

axSun(inx) = 0; aySun(inx) = 0; azSun(inx) = 0;

xJupiter(inx) = 4.505316299005550E+00 * AU;
yJupiter(inx) = -2.163555523654128E+00 * AU;
zJupiter(inx) = -9.190511451296739E-02 * AU;

vxJupiter(inx) = 3.174273032056012E-03 * AU / day;
vyJupiter(inx) = 7.161028024427482E-03 * AU / day;
vzJupiter(inx) = -1.007935797387438E-04 * AU / day;

axJupiter(inx) = 0; ayJupiter(inx) = 0; azJupiter(inx) = 0;

rSun(inx) = sqrt(xSun(inx)^2+ySun(inx)^2+zSun(inx)^2);
rJupiter(inx) = sqrt(xJupiter(inx)^2+yJupiter(inx)^2+zJupiter(inx)^2);

alphaSun(inx) = acos(xSun(inx)/rSun(inx));
betaSun(inx) = acos(ySun(inx)/rSun(inx));
gammaSun(inx) = acos(zSun(inx)/rSun(inx));

alphaJupiter(inx) = acos(xJupiter(inx)/rJupiter(inx));
betaJupiter(inx) = acos(yJupiter(inx)/rJupiter(inx));
gammaJupiter(inx) = acos(zJupiter(inx)/rJupiter(inx));

while t(inx) < 36*86400
    inx = inx + 1;
    
    alphaSun(inx) = acos(xSun(inx-1) / rSun(inx-1));
    betaSun(inx) = acos(ySun(inx-1) / rSun(inx-1));
    gammaSun(inx) = acos(zSun(inx-1) / rSun(inx-1));
    
    axSun(inx) = -(G*MTotal)/(rSun(inx-1)^2)*cos(alphaSun(inx));
    aySun(inx) = -(G*MTotal)/(rSun(inx-1)^2)*cos(betaSun(inx));
    azSun(inx) = -(G*MTotal)/(rSun(inx-1)^2)*cos(gammaSun(inx));
    vxSun(inx) = vxSun(inx-1) + axSun(inx) * dt;
    vySun(inx) = vySun(inx-1) + aySun(inx) * dt;
    vzSun(inx) = vzSun(inx-1) + azSun(inx) * dt;
    xSun(inx) = xSun(inx-1) + vxSun(inx) * dt;
    ySun(inx) = ySun(inx-1) + vySun(inx) * dt;
    zSun(inx) = zSun(inx-1) + vzSun(inx) * dt;
    
    alphaJupiter(inx) = acos(xJupiter(inx-1) / rJupiter(inx-1));
    betaJupiter(inx) = acos(yJupiter(inx-1) / rJupiter(inx-1));
    gammaJupiter(inx) = acos(zJupiter(inx-1) / rJupiter(inx-1));
    
    axJupiter(inx) = -(G*MTotal)/(rJupiter(inx-1)^2)*cos(alphaSun(inx));
    ayJupiter(inx) = -(G*MTotal)/(rJupiter(inx-1)^2)*cos(betaSun(inx));
    azJupiter(inx) = -(G*MTotal)/(rJupiter(inx-1)^2)*cos(gammaSun(inx));
    vxJupiter(inx) = vxJupiter(inx-1) + axJupiter(inx) * dt;
    vyJupiter(inx) = vyJupiter(inx-1) + ayJupiter(inx) * dt;
    vzJupiter(inx) = vzJupiter(inx-1) + azJupiter(inx) * dt;
    xJupiter(inx) = xJupiter(inx-1) + vxJupiter(inx) * dt;
    yJupiter(inx) = yJupiter(inx-1) + vyJupiter(inx) * dt;
    zJupiter(inx) = zJupiter(inx-1) + vzJupiter(inx) * dt;

    rSun(inx) = sqrt(xSun(inx)^2+ySun(inx)^2+zSun(inx)^2);
    rJupiter(inx) = sqrt(xJupiter(inx)^2+yJupiter(inx)^2+zJupiter(inx)^2);
    
    t(inx) = t(inx-1) + dt;
    
     if rem(t(inx), 1000) == 0
         clf
         hold on
         view(3)
         plot3(xSun(inx), ySun(inx), zSun(inx), 'bo')
         plot3(xJupiter(inx), yJupiter(inx), zJupiter(inx), 'ro')
         xlim([-6*AU, 6*AU])
         ylim([-6*AU, 6*AU])
         zlim([-2*AU, 2*AU])
         grid on
         pause(0.01)
     end
end