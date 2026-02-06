function makePlots(iTile, results, trajTransMod, trajSusMod, par, desc)  

% Make a phase plot and time series plots for this parameter combination in
% the specified tiles of the current figure

% Plot settings and labels
lbls = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)" ];
colOrd = colororder;
greyCol =  0.8*[1 1 1];
ls = ["-", "--"];
SSSmarksize = 22;
% Index of trajcetory solution to place arrow at
iArr = 30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot phase plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile(iTile);

% Plot vector field
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



% Plot stable/unstable manifolds (if the saddle EE exists)
if ~isempty(results.Ystable1)
    ha.ColorOrderIndex = 6;
    plot(results.Ystable1(:, 1), results.Ystable1(:, 2), '-', 'LineWidth', 2, 'HandleVisibility', 'off')
    ar = orbitArrow(flipud(results.Ystable1(:, 1)), flipud(results.Ystable1(:, 2)), "horiz", 0.75, colOrd(6, :));
end
if ~isempty(results.Ystable2)
    ha.ColorOrderIndex = 6;
    plot(results.Ystable2(:, 1), results.Ystable2(:, 2), '-', 'LineWidth', 2, 'HandleVisibility', 'off')
    ar = orbitArrow(flipud(results.Ystable2(:, 1)), flipud(results.Ystable2(:, 2)), "horiz", 0.5, colOrd(6, :));
end
if ~isempty(results.Yunstable1)
    ha.ColorOrderIndex = 7;
    plot(results.Yunstable1(:, 1), results.Yunstable1(:, 2), '-', 'LineWidth', 2, 'HandleVisibility', 'off')
    ar = orbitArrow(results.Yunstable1(:, 1), results.Yunstable1(:, 2), "vert", 0.7, colOrd(7, :));
end
if ~isempty(results.Yunstable2)
    ha.ColorOrderIndex = 7;
    plot(results.Yunstable2(:, 1), results.Yunstable2(:, 2), '-', 'LineWidth', 2, 'HandleVisibility', 'off')
    ar = orbitArrow(results.Yunstable2(:, 1), results.Yunstable2(:, 2), "vert", 0.7, colOrd(7, :));
end



% Plot trajectories from the different initial conditions
nICs = length(trajTransMod);
for iIC = 1:nICs
    plot(trajTransMod(iIC).S, trajTransMod(iIC).B, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
    ar = orbitArrow( trajTransMod(iIC).S, trajTransMod(iIC).B, "ind", iArr, [0 0 0]);
end


% Plot equilibria
% Plot NDFE
plot(1, 0, 'ko', 'HandleVisibility', 'off')
% Plot BDFEs (if real)
R0 = par.Beta/par.Gamma;
R0TC2 = 1/(1-par.qc*results.Bstar(1));
if isreal(results.Bstar(1))
    if R0 < R0TC2
        % BDFE+ is stable when R0 < R0TC2
        plot(1, results.Bstar(1), 'k.', 'MarkerSize', SSSmarksize, 'HandleVisibility', 'off')
    else
        % Otherwise BDFE is unstable 
        plot(1, results.Bstar(1), 'ko', 'HandleVisibility', 'off')
    end
    % BDFE- is always unstable when it exists
    plot(1, results.Bstar(2), 'ko', 'HandleVisibility', 'off')
end
% Plot endemic equilibria
plot(results.EE(1, 1), results.EE(2, 1), 'k.', 'MarkerSize', SSSmarksize, 'HandleVisibility', 'off')            
plot(results.EE(1, 2), results.EE(2, 2), 'ko', 'HandleVisibility', 'off')            
plot(results.EE(1, 3), results.EE(2, 3), 'k.', 'MarkerSize', SSSmarksize, 'HandleVisibility', 'off')  


xlim([0.4 1])
ylim([0 1])
xlabel('S')
ylabel('B')
title("R_0 = " + sprintf('%.2f', R0) + ", " + desc + newline + lbls(iTile))


if iTile == 1
    % Add phase plot legend to the first row in the figure
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot time series for transmission-modulated model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iTile = iTile + 1; 
nexttile(iTile);
for iIC = 1:2
    ha = gca;
    ha.ColorOrderIndex = 1;
    plot(trajTransMod(iIC).t, 1-trajTransMod(iIC).S, 'LineStyle', ls(iIC), 'LineWidth', 1.5)
    hold on
    plot(trajTransMod(iIC).t, trajTransMod(iIC).B, 'LineStyle', ls(iIC), 'LineWidth', 1.5)
end
% Plot equilibrium without behaviour I = 1-1/R0
yline(1-1/R0, 'k:');
ylim([0 1])
xlabel('time')
title(lbls(iTile))
if iTile == 2
    % Add legend and figure title for first row of plots
    sgtitle("\tau/\alpha = " + sprintf('%.1f', par.Tau/par.Alpha) + ", q = " +sprintf('%.1f', par.qc) );
    lgd = legend('I', 'B', 'Location', 'northeast');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot time series for susceptibility-modulated model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iTile = iTile+1;
nexttile(iTile);
for iIC = 1:2
    ha = gca;
    ha.ColorOrderIndex = 1;
    plot(trajSusMod(iIC).t, 1-trajSusMod(iIC).S, 'LineStyle', ls(iIC), 'LineWidth', 1.5)
    hold on
    plot(trajSusMod(iIC).t, trajSusMod(iIC).B, 'LineStyle', ls(iIC), 'LineWidth', 1.5)
    % For the susceptibility model plot S_B/S 
    plot(trajSusMod(iIC).t, trajSusMod(iIC).SB./trajSusMod(iIC).S, 'LineStyle', ls(iIC), 'LineWidth', 1.5)
end
% Plot equilibrium without behaviour I = 1-1/R0
yline(1-1/R0, 'k:');
ylim([0 1])
xlabel('time')
title(lbls(iTile))

% Add legend to time series plots for the first case in each figure
 if iTile == 3
    lgd = legend('I', 'B', 'S_B/S', 'Location', 'northeast');
 end

