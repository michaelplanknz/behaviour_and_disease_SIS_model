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

% Time range        1000 for slow, 300 for fast
tSpan = 0:1000;

% Spacing for vector field arrow length plotting
dx = 0.1;

% Power for vector field arrow length plotting
pQuiv = 0.5;      

% Set ICs
IC1 = [(1-B01)*S0; B01*S0; B01 ];
IC2 = [(1-B02)*S0; B02*S0; B02 ];


h = figure(1);
h.Position = [ 89   126   931   864];
tiledlayout(3, 3, "TileSpacing", "compact");

h = figure(2);
h.Position = [     164   409   931   562];
tiledlayout(2, 3, "TileSpacing", "compact");

h = figure(3);
h.Position = [  300   126   931   864];
tiledlayout(3, 3, "TileSpacing", "compact");

% For each parameter combo, specify the figure and tile to put the results
% in
figNum = [ones(1, 3), 2*ones(1, 2), 3*ones(1, 3)];
tileNum = [1, 4, 7, 1, 4, 1, 4, 7];


desc = ["EE low B", "bistable EEs", "EE high B", "bistable EEs", "EE high B", "bistable EE/BDFE", "bistable EE/BDFE", "BDFE" ];
    
lbls = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)" ];

% Run and plot each parameter combination in turn
for iCase = 1:nCases
    par = getPar(iCase);

    h = figure(figNum(iCase));



    % Do the same thing for the susceptibility and transmission versions of the
    % model
    for iType = 1:2
        % Select which ODE model to run
        if iType == 1
            ODE_handle = @(t, y)myODEs_trans(t, y, par);
        else
            ODE_handle = @(t, y)myODEs_suscept(t, y, par);
        end        
        
        % Run the two initial conditions
        [t1, Y1] = ode45(ODE_handle, tSpan, IC1);
        [t2, Y2] = ode45(ODE_handle, tSpan, IC2);
    
        % Extract epi variables
        SN1 = Y1(:, 1);
        SB1 = Y1(:, 2);
        B1 = Y1(:, 3);
        S1 = SN1+SB1;
    
        SN2 = Y2(:, 1);
        SB2 = Y2(:, 2);
        B2 = Y2(:, 3);
        S2 = SN2+SB2;
    

        % Plot phase plane (for trans-modulated model only)
        if iType == 1

            % Evaluate vector field over a grid
            sx = 0:dx:1;
            by = 0:dx:1;
            [SX, BY] = meshgrid(sx, by);
            [U, V] = myVectorField(SX, BY, par);

            % Plot arrows with a power p of their real length so long
            % vectors don't dominate too much
            Z = sqrt(U.^2+V.^2);
            U = U.*Z.^(pQuiv-1);
            V = V.*Z.^(pQuiv-1);

            iTile = tileNum(iCase);
            nexttile(iTile);
            quiver(sx, by, U, V, 'Color', 0.8*[1 1 1]);
            hold on
            ha = gca;
            ha.ColorOrderIndex = 1;
            plot(S1, B1)
            ha.ColorOrderIndex = 1;
            plot(S2, B2)

            xlim([0 1])
            ylim([0 1])
            xlabel('S')
            ylabel('B')
            title("R_0 = " + sprintf('%.2f', par.Beta/par.Gamma) + ", " + desc(iCase) + newline + lbls(iTile))
        end


        % Plot time series
        iTile = tileNum(iCase) + iType; 
        nexttile(iTile);
        plot(t1, 1-S1)
        hold on
        plot(t1, B1)
        if iType == 2
            % Only need to plot S_B/S for the susceptibility model as it is
            % identical to B in the transmission model
            plot(t1, SB1./S1)
        end
        ha = gca;
        ha.ColorOrderIndex = 1;
        plot(t2, 1-S2, '--')
        plot(t2, B2, '--')
        if iType == 2
            % Only need to plot S_B/S for the susceptibility model as it is
            % identical to B in the transmission model
            plot(t2, SB2./S2, '--')
        end
        % Plot equilibrium without behaviour I = 1-1/R0
        yline(1-par.Gamma/par.Beta, 'k:');
        ylim([0 1])
        xlabel('time (days)')
        title(lbls(iTile))
    end



     if ismember(iCase, [1 4 6])
        sgtitle("\tau/\alpha = " + sprintf('%.1f', par.Tau/par.Alpha) + ", q = " +sprintf('%.1f', par.q) );
        nexttile(2);  
        lgd = legend('I', 'B', 'Location', 'northeast');
        nexttile(3);
        lgd = legend('I', 'B', 'S_B/S', 'Location', 'northeast');
    end
end

for iFig = 1:3
    h = figure(iFig);
    fName = figFolder + "fig_cases" + iFig + ".png";
    saveas(h, fName);
end