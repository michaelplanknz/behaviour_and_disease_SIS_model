clear 
close all

% Run the SIRS model and plot the results

% Figure folder
figFolder = "figures/";

% Parameter values
par.Gamma = 1;        
par.Alpha = 0.25;            
par.Chi = 20 * par.Alpha;
par.Tau = 5 * par.Alpha;
par.Beta = 2 * par.Gamma; 
par.qs = 0;
par.qc = 0.7;
par.w = 0.1;

% Initial conditions
i0 = 0.001;
IC = [1-i0; 0; i0; 0; 0];
tSpan = [0, 100];

% Parameters for equivalent model with no behaviour
parN = par;
parN.Chi = 0;
parN.Tau = 0;

% Solve ODEs with behaviour
[t1, Y1] = ode45(@(t, y)SIRSmodelODEs(t, y, par), tSpan, IC);
I1 = Y1(:, 3)+Y1(:, 4);
B1 = Y1(:, 5);

% Solve ODEs without behaviour
[t2, Y2] = ode45(@(t, y)SIRSmodelODEs(t, y, parN), tSpan, IC);
I2 = Y2(:, 3)+Y2(:, 4);
B2 = Y2(:, 5);


% Plot results
h = figure(1);
h.Position = [ 2052         688         842         327];
tiledlayout(1, 2, "TileSpacing", "compact");
nexttile(1);
plot(t1, I1, t2, I2)
hold on
xlabel('time')
ylabel('I(t)')
legend('with behaviour', 'without behaviour')
title('(a)')
nexttile(2);
plot(t1, B1, t2, B2)
hold on
xlabel('time')
ylabel('B(t)')
title('(b)')

% Save
fName = figFolder + "SIRS.png";
saveas(h, fName);




