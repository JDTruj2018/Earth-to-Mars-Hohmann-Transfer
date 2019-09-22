%% Author: Jered Dominguez-Trujillo
% Description: Aerospace Propulsion Homework 3

%% Assign Variables
clear; clc; clf; close all;
WaterM = 18.01528;
WaterR = 8.3144598/WaterM; AirR = 0.2869;
g = 9.80665; Ttp = 216.65; Ptp = 22573; Pref = 101325; Tref = 288.15;
To = 3500; Po = 300*10^5 + Pref;
NozzleArea = 0.062912356383544;
AeRatio = 200;
% AeRatio = 40;
AExit = AeRatio * NozzleArea; AThroat = AExit / AeRatio;

%% Define Cp, Cv, and Gamma as functions of temperature
Cp = @(T) ((T >= 0.5 && T <= 1.7).*((30.092 + (6.832514.*T) +...
    (6.793435.*(T.^2)) - (2.53448.*(T.^3)) + (0.082139./(T.^2)))...
    ./ WaterM)) + ((T > 1.7 && T <= 6).*((41.96426 + (8.622053.*T) -...
    (1.49978.*(T.^2)) + (0.098119.*(T.^3)) + (-11.1576./(T.^2)))...
    ./ WaterM));

Cv = @(T) Cp(T) - WaterR;

Gamma = @(T) Cp(T)/Cv(T);

%% Establish A/A* as a function of Pe/Po and Me; accomodating variation of Gamma (T)
set(gcf, 'Position', get(0, 'Screensize'));
hold on

TempArray = [0.5, 1:1:6];
leg = cell(length(TempArray), 1);
count = 1;

for i = TempArray
    G = Gamma(i);
    leg{count} = ['Gamma = ' num2str(G)];
    Me = (1:0.001:6);
    A = (1./Me).*((2+(G - 1).*(Me.^2))./(G + 1)).^((G + 1)./(2.*(G-1)));
    count = count + 1;
    plot(A, Me)
end

title('$$\frac{A}{A*}\ as\ a\ Function\ of\ Mach\ Number$$',...
    'interpreter', 'latex', 'fontsize', 20)
ylabel('$$Mach\ Number$$', 'interpreter', 'latex', 'fontsize', 16)
xlabel('$$\frac{A}{A*}$$', 'interpreter', 'latex', 'fontsize', 16)
legend(leg, 'location', 'southeast')
set(findall(gcf,'Type','Line'),'LineWidth', 2)
grid on

print('Area Ratio as Function of Mach Number', '-dpng')

close all; clf;

%% Calculating Design Exit Pressure and Altitude
criticalRatio = 300;
Pa = 200;
while criticalRatio > AeRatio
    Pa = Pa + 10;
    Temp = To / (((Po/Pa)^(Gamma(To/1000)-1)/Gamma(To/1000)));

    G = Gamma(Temp/1000);
    Me = 1; ARatioList = (1:0.025:1000);

    inx = 1;
    guess = (1/Me)*((2+(G-1)*Me^2)/(G+1))^((G+1)/(2*(G-1)));

    for ratio = ARatioList
        while guess < ratio
            Me = Me + 0.0025;
            guess = (1/Me)*((2+(G-1)*Me^2)/(G+1))^((G+1)/(2*(G-1)));
        end

        PRatio(inx) = 1/((1+((G-1)/2)*Me^2)^(G/(G-1)));
        Term1 = 2*G/(G-1);
        Term2 = (2/(G+1))^((G+1)/(G-1));
        Term3 = (1 - (PRatio(inx))^((G-1)/G));
        Term4 = sqrt(Term1 * Term2 * Term3);
        Cf(inx) = Term4 + (PRatio(inx) - Pa/Po) * ratio;
        inx = inx + 1;
    end

    criticalRatio = ARatioList(Cf==max(Cf));
end

height = ((1 - (Pa/Pref)^(1/5.2561))/(0.0065))*Tref;

%% Calculating Thrust
DesignThrust = Po * AThroat * max(Cf);
GroundThrust = DesignThrust - ((Pref - Pa)*AExit);
VacuumThrust = DesignThrust - ((0 - Pa)*AExit);

%% Printing Results
fprintf('Exit Temperature: %.2f K\n', Temp);
fprintf('Gamma: %.4f \n', G);
fprintf('Design Pressure: %.0f Pa\n', Pa);
fprintf('Design Altitude: %.2f m\n', height);
fprintf('Design Thrust: %.2f N\n', DesignThrust);
fprintf('Ground Thrust: %.2f N\n', GroundThrust);
fprintf('Vacuum Thrust: %.2f N\n', VacuumThrust);
fprintf('Pe/Po : %.6f \n', Pa/Po);

%% Plotting Thrust Coefficient of Design as a Function or A/A*
set(gcf, 'Position', get(0, 'Screensize'));

plot(ARatioList, Cf, 'g')

% line([AeRatio AeRatio], [0 max(Cf)], 'Color', 'k');
% text(80, max(Cf)/1.5,'$$ \frac{A}{A*}\ =\ 77.5$$', 'interpreter', 'latex')
% text(60.5, max(Cf)+0.05, ['Max \tau_{c} =' num2str(max(Cf))]);
% text(500, 1.55, '$$Design\ \frac{P_{e}}{P_{o}}$$', 'interpreter', 'latex')

hold on
TGamma = G;

for PaPo = [0.0001, 0.001, 0.002, 0.005, 0.01]
    Me = 1;
    Term = (1/Me)*((2+(TGamma-1)*Me^2)/(TGamma+1))^((TGamma+1)/(2*(TGamma-1)));
    inx = 1;
    for a = ARatioList
        while Term < a
            Me = Me + 0.01;
            Term = (1/Me)*((2+(TGamma-1)*Me^2)/(TGamma+1))^((TGamma+1)/(2*(TGamma-1)));
        end

        Me;

        ratio = 1/((1+((TGamma-1)/2)*Me^2)^(TGamma/(TGamma-1)));

        Term1 = 2*TGamma/(TGamma-1);
        Term2 = (2/(TGamma+1))^((TGamma+1)/(TGamma-1));
        Term3 = (1-(ratio)^((TGamma-1)/TGamma));
        Term4 = sqrt(Term1*Term2*Term3);
        Cf(inx) = Term4+(ratio-(PaPo))*a;
        inx = inx + 1;
    end
    plot(ARatioList, Cf)
    axis tight
end

title('$$Thrust\ Coefficient\ as\ a\ Function\ of\ \frac{A}{A*}$$',...
    'interpreter', 'latex', 'fontsize', 20)
xlabel('$$\frac{A}{A*}$$', 'interpreter', 'latex', 'fontsize', 16)
ylabel('$$Thrust\ Coefficient$$', 'interpreter', 'latex', 'fontsize', 16)
legend({'Design ^{Pe}/_{Po} = .000175', '^{Pe}/_{Po} = 0.0001',...
    '^{Pe}/_{Po} = 0.001', '^{Pe}/_{Po} = 0.002', '^{Pa}/_{Po} = 0.005',...
    '^{Pe}/_{Po} = 0.01'}, 'location', 'southwest')
ylim([0 2])
set(findall(gcf,'Type','Line'),'LineWidth', 2)
grid on

print('Thrust Coefficient 200', '-dpng')