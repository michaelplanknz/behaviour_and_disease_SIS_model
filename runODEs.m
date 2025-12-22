clear 
close all

% Figure folder
figFolder = "figures/";

% Number of parameter combinations 
nCases = 8;

% Initial conditions for B
B01 = 0;
B02 = 0.9;

% Initial condition for S
S0 = 0.95;

% Time range        3000 for slow, 300 for fast
tSpan = 0:300;

% Set ICs
IC1 = [(1-B01)*S0; B01*S0; B01 ];
IC2 = [(1-B02)*S0; B02*S0; B02 ];


h = figure(1);
h.Position = [   564          64        1165         934];
tiledlayout(3, 3, "TileSpacing", "compact");

% Tiles to plot cases 1-98 in and description of dynamics
iPlot = [1, 2, 3, 4, 5, 7, 8, 9];
desc = ["EE low B", "bistable EEs", "EE high B", "bistable EEs", "EE high B", "bistable EE/BDFE", "bistable EE/BDFE", "BDFE" ];

lbls = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)"];

% Run and plot each parameter combination in turn
for iCase = 1:nCases
    par = getPar(iCase);
    
    % Run the two initial conditions
    [t1, Y1] = ode45(@(t, y)myODEs(t, y, par), tSpan, IC1);
    [t2, Y2] = ode45(@(t, y)myODEs(t, y, par), tSpan, IC2);

    % Extract epi variables
    SN1 = Y1(:, 1);
    SB1 = Y1(:, 2);
    B1 = Y1(:, 3);
    S1 = SN1+SB1;

    SN2 = Y2(:, 1);
    SB2 = Y2(:, 2);
    B2 = Y2(:, 3);
    S2 = SN2+SB2;

    % Plot
    nexttile(iPlot(iCase));
    plot(t1, S1, t1, B1, t1, SB1./S1)
    hold on
    ha = gca;
    ha.ColorOrderIndex = 1;
    plot(t2, S2, '--', t2, B2, '--', t2, SB2./S2, '--')
    ylim([0 1])
    grid on
    xlabel('time (days)')
    title(lbls(iCase) + " R_0 = " + sprintf('%.1f', par.Beta/par.Gamma) +  ", \tau = " + sprintf('%.1f', par.Tau) +  ", q = " + sprintf('%.1f\n', par.q) + desc(iCase))

end


l = legend('S', 'B', 'S_B/S', 'Location', '');
l.Layout.Tile = 6;

fName = figFolder + "time_series_diff_susceptibility.png";
saveas(h, fName);

