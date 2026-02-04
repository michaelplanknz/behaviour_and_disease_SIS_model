clear 
close all

% Figure folder
figFolder = "figures/";

% Number of parameter combinations 
nCases = 8;

% Initial conditions for B (all will be plotted in phase plane but only first two will be plotted as time series)
B0 = [0, 0.9, 0.3, 0.6];

% Initial condition for S
S0 = 0.95;

% Time range        1000 for slow, 300 for fast
tSpan = 0:1000;

% Time span and perturbation size for finding stable/unstable manifolds
tManifold = [0, 1e4];
pert = 1e-4;
opts = optimoptions('fsolve', 'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10, 'Display' ,'off', 'Algorithm', 'levenberg-marquardt');

% Spacing for nullcline plotting (dx) and vector field arrow length plotting (dxquiv and dyquiv)
dx = 0.001;
dxquiv = 0.05;
dyquiv = 0.1;


% Power for vector field arrow length plotting (1 for actual length, 0 for
% all the same length)
pQuiv = 0.5;      

% Set ICs
IC = [(1-B0)*S0; B0*S0; B0 ];


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
greyCol =  0.8*[1 1 1];
SSSmarksize = 22;

% Run and plot each parameter combination in turn
for iCase = 1:nCases
    par = getPar(iCase);

    % Calculate R0
    R0 = par.Beta/par.Gamma;

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
        
        % Run the four initial conditions
        nICs = length(B0);
        for iIC = 1:nICs
            traj(iIC) = solveODE(ODE_handle, tSpan, IC(:, iIC));
        end

        % Plot phase plane (for trans-modulated model only)
        if iType == 1

            results = getPhasePlotOutputs(par, dx, dxquiv, dyquiv, pQuiv, tManifold, pert, opts);

            % Plot vector field
            iTile = tileNum(iCase);
            nexttile(iTile);
            quiver(results.sx, results.by, results.U, results.V, 'Color', greyCol, 'HandleVisibility', 'off');
            hold on

            ha = gca;

            % Plot S nullcline
            ha.ColorOrderIndex = 4;
            plot(results.Snull1, results.Bnull1, '--', 'LineWidth', 2, 'HandleVisibility', 'off')
            % Plot the other branch of the S nullcline which is the line S=1
            ha.ColorOrderIndex = 4;
            plot([1 1], [0 1], '--', 'LineWidth', 2, 'HandleVisibility', 'off')

            % Plot B nullcline
            plot(results.Snull2, results.Bnull2, '--', 'LineWidth', 2, 'HandleVisibility', 'off')
            ha.ColorOrderIndex = 5;

            % Plot equilibria
            % Plot NDFE
            plot(1, 0, 'ko', 'HandleVisibility', 'off')
            % Plot BDFEs (if real)
            if isreal(results.Bstar(1))
                if ismember(iCase, [6, 7, 8])
                    % BDFE+ is stable in these cases
                    plot(1, results.Bstar(1), 'k.', 'MarkerSize', SSSmarksize, 'HandleVisibility', 'off')
                else
                    % Otherwise BDFE is unstable in these cases
                    plot(1, results.Bstar(1), 'ko', 'HandleVisibility', 'off')
                end
                % BDFE is always unstable when it exists
                plot(1, results.Bstar(2), 'ko', 'HandleVisibility', 'off')
            end
            % Plot endemic equilibria
            plot(results.EE(1, 1), results.EE(2, 1), 'k.', 'MarkerSize', SSSmarksize, 'HandleVisibility', 'off')            
            plot(results.EE(1, 2), results.EE(2, 2), 'ko', 'HandleVisibility', 'off')            
            plot(results.EE(1, 3), results.EE(2, 3), 'k.', 'MarkerSize', SSSmarksize, 'HandleVisibility', 'off')            

            % Plot stable/unstable manifolds (if the saddle EE exists)
            if ~isnan(results.EE(1, 2))
                ha.ColorOrderIndex = 6;
                plot(results.Ystable1(:, 1), results.Ystable1(:, 2), '-', 'LineWidth', 2, 'HandleVisibility', 'off')
                ha.ColorOrderIndex = 6;
                plot(results.Ystable2(:, 1), results.Ystable2(:, 2), '-', 'LineWidth', 2, 'HandleVisibility', 'off')
                plot(results.Yunstable1(:, 1), results.Yunstable1(:, 2), '-', 'LineWidth', 2, 'HandleVisibility', 'off')
                ha.ColorOrderIndex = 7;
                plot(results.Yunstable2(:, 1), results.Yunstable2(:, 2), '-', 'LineWidth', 2, 'HandleVisibility', 'off')
            end

            % Plot trajectories from the two initial conditions
            %ha.ColorOrderIndex = 1;
            for iIC = 1:nICs
                plot(traj(iIC).S, traj(iIC).B, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off')
            end
            xlim([0.4 1])
            ylim([0 1])
            xlabel('S')
            ylabel('B')
            title("R_0 = " + sprintf('%.2f', R0) + ", " + desc(iCase) + newline + lbls(iTile))
        end


        % Plot time series
        iTile = tileNum(iCase) + iType; 
        nexttile(iTile);
        ls = ["-", "--"];
        for iIC = 1:2
            ha = gca;
            ha.ColorOrderIndex = 1;
            plot(traj(iIC).t, 1-traj(iIC).S, 'LineStyle', ls(iIC), 'LineWidth', 1.5)
            hold on
            plot(traj(iIC).t, traj(iIC).B, 'LineStyle', ls(iIC), 'LineWidth', 1.5)
            if iType == 2
                % Only need to plot S_B/S for the susceptibility model as it is
                % identical to B in the transmission model
                plot(traj(iIC).t, traj(iIC).SB./traj(iIC).S, 'LineStyle', ls(iIC), 'LineWidth', 1.5)
            end
        end
        % Plot equilibrium without behaviour I = 1-1/R0
        yline(1-1/R0, 'k:');
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
    nexttile(1);
    ha = gca;
    quiver(0, 0, 0, 0, 'color', greyCol, 'DisplayName', 'vector field')
    ha.ColorOrderIndex = 4;
    plot(nan, nan, '--', 'LineWidth', 2, 'DisplayName', 'S nullcline')
    plot(nan, nan, '--', 'LineWidth', 2, 'DisplayName', 'B nullcline')
    plot(nan, nan, '-', 'LineWidth', 2, 'DisplayName', 'stable mnfld')
    plot(nan, nan, '-', 'LineWidth', 2, 'DisplayName', 'unstable mnfld')
    plot(nan, nan, 'k-', 'LineWidth', 1.5, 'DisplayName', 'trajectory')
    plot(nan, nan, 'k.', 'MarkerSize', SSSmarksize, 'DisplayName', 'stable eq.')
    plot(nan, nan, 'ko', 'DisplayName', 'unstable eq.')
    lgd = legend('Location', 'westoutside');

    fName = figFolder + "fig_cases" + iFig + ".png";
    saveas(h, fName);
end