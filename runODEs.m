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

% Set ICs
IC1 = [(1-B01)*S0; B01*S0; B01 ];
IC2 = [(1-B02)*S0; B02*S0; B02 ];



% Do the same thing for the susceptibility and transmission versions of the
% model
for iType = 1:2
    
    h = figure(iType);
    h.Position = [   564          64        1165         934];
    tiledlayout(3, 3, "TileSpacing", "compact");
    
    % Tiles to plot cases 1-98 in and description of dynamics
    iPlot = [1, 2, 3, 4, 5, 7, 8, 9];
    desc = ["EE low B", "bistable EEs", "EE high B", "bistable EEs", "EE high B", "bistable EE/BDFE", "bistable EE/BDFE", "BDFE" ];
    
    lbls = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)"];
    
    % Run and plot each parameter combination in turn
    for iCase = 1:nCases
        par = getPar(iCase);
    
        % Select which ODE model to run
        if iType == 1
            ODE_handle = @(t, y)myODEs_suscept(t, y, par);
        else
            ODE_handle = @(t, y)myODEs_trans(t, y, par);
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
    
        % Plot
        nexttile(iPlot(iCase));
        plot(t1, S1)
        hold on
        plot(t1, B1)
        if iType == 1
            % Only need to plot S_B/S for the susceptibility model as it is
            % identical to B in the transmission model
            plot(t1, SB1./S1)
        end
        ha = gca;
        ha.ColorOrderIndex = 1;
        plot(t2, S2, '--')
        plot(t2, B2, '--')
        if iType == 1
            % Only need to plot S_B/S for the susceptibility model as it is
            % identical to B in the transmission model
            plot(t2, SB2./S2, '--')
        end
        ylim([0 1])
        grid on
        xlabel('time (days)')
        title(lbls(iCase) + " R_0 = " + sprintf('%.1f', par.Beta/par.Gamma) +  ", \tau = " + sprintf('%.3f', par.Tau) +  ", q = " + sprintf('%.1f\n', par.q) + desc(iCase))
    
    end
    
    
    if iType == 1   
        l = legend('S', 'B', 'S_B/S');
    else
        l = legend('S', 'B');
    end
    l.Layout.Tile = 6;
    
    if iType == 1
        sgtitle('Behaviour-dependent susceptibility')
        fName = figFolder + "time_series_diff_susceptibility.png";
    else
        sgtitle('Behaviour-dependent transmission')
        fName = figFolder + "time_series_diff_transmission.png";
    end
    saveas(h, fName);

end