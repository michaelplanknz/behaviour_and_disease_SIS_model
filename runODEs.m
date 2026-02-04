clear 
close all

% Figure folder
figFolder = "figures/";

% Number of parameter combinations 
nCases = 8;

% Initial conditions for B (all will be plotted in phase plane but only first two will be plotted as time series)
B01 = 0;
B02 = 0.9;
B03 = 0.3;
B04 = 0.6;

% Initial condition for S
S0 = 0.95;

% Time range        1000 for slow, 300 for fast
tSpan = 0:1000;

% Spacing for nullcline plotting (dx) and vector field arrow length plotting (dxquiv and dyquiv)
dx = 0.01;
dxquiv = 0.05;
dyquiv = 0.1;


% Power for vector field arrow length plotting
pQuiv = 0.5;      

% Set ICs
IC1 = [(1-B01)*S0; B01*S0; B01 ];
IC2 = [(1-B02)*S0; B02*S0; B02 ];
IC3 = [(1-B03)*S0; B03*S0; B03 ];
IC4 = [(1-B04)*S0; B04*S0; B04 ];


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
        [t1, Y1] = ode45(ODE_handle, tSpan, IC1);
        [t2, Y2] = ode45(ODE_handle, tSpan, IC2);
        [t3, Y3] = ode45(ODE_handle, tSpan, IC3);
        [t4, Y4] = ode45(ODE_handle, tSpan, IC4);
    
        % Extract epi variables
        SN1 = Y1(:, 1);
        SB1 = Y1(:, 2);
        B1 = Y1(:, 3);
        S1 = SN1+SB1;
    
        SN2 = Y2(:, 1);
        SB2 = Y2(:, 2);
        B2 = Y2(:, 3);
        S2 = SN2+SB2;

        SN3 = Y3(:, 1);
        SB3 = Y3(:, 2);
        B3 = Y3(:, 3);
        S3 = SN3+SB3;
    
        SN4 = Y4(:, 1);
        SB4 = Y4(:, 2);
        B4 = Y4(:, 3);
        S4 = SN4+SB4;

        % Plot phase plane (for trans-modulated model only)
        if iType == 1

            % Evaluate vector field over a grid
            sx = 0:dxquiv:1;
            by = 0:dyquiv:1;
            [SX, BY] = meshgrid(sx, by);
            [U, V] = myVectorField(SX, BY, par);

            % Scale arrows with a power p of their initial length so long
            % vectors don't dominate too much
            Z = sqrt(U.^2+V.^2);
            U = U.*Z.^(pQuiv-1);
            V = V.*Z.^(pQuiv-1);

            % Get (B,S) coords for the S nullcline
            Snull1 = 0:dx:1;
            Bnull1 = 1/(1-par.q) * (1 - 1./(R0*Snull1));

            % Get (B,S) coords for the B nullcline
            Bnull2 = 0:dx:1;
            Snull2 = 1 + (par.Tau*Bnull2.^2 - par.Alpha*Bnull2./(1-Bnull2))/par.Chi;
            
            % Get coords for BDFE
            Bstar = 0.5 * (1 + [1, -1]*sqrt(1-4*par.Alpha/par.Tau));

            % Find endemic equilibria
            % First get good initial conditions for root-solving by moving along
            % the S nullcline and seeing where it crosses the B hullcline
            nx = length(Snull1);
            xr0 = [];
            leftFlag = 1;
            for ix = 1:nx
                % Only consider values in [0, 1]
                if Bnull1(ix) >= 0 & Bnull1(ix) <= 1
                    leftFlagPrev = leftFlag;
                    % Compute flag variable indicating whether the current
                    % point on the S nullclnie is left of the B nullcline
        	        leftFlag = Snull1(ix) < interp1(Bnull2, Snull2, Bnull1(ix));
                    if leftFlag ~= leftFlagPrev
                        % Crossing detected - append initial condition to
                        % xr0
                        xr0 = [xr0, [Snull1(ix); Bnull1(ix)] ];
                    end
                end
            end
            nICs = min(3, size(xr0, 2));
            EE = nan(2, 3);
            for iIC = 1:nICs
                EE(:, iIC) = fsolve(@(x)mySteadyState(x, par), xr0(:, iIC) );
            end


            % Plot vector field
            iTile = tileNum(iCase);
            nexttile(iTile);
            quiver(sx, by, U, V, 'Color', 0.8*[1 1 1]);
            hold on

            ha = gca;

            % Plot S nullcline
            ha.ColorOrderIndex = 2;
            plot(Snull1, Bnull1, 'LineWidth', 2)
            ha.ColorOrderIndex = 2;
            plot([1 1], [0 1], 'LineWidth', 2)

            % Plot B nullcline
            plot(Snull2, Bnull2, 'LineWidth', 2)
            ha.ColorOrderIndex = 3;

            % Plot equilibria
            % Plot NDFE
            plot(1, 0, 'bo')
            % Plot BDFEs (if real)
            if isreal(Bstar(1))
                if ismember(iCase, [6, 7, 8])
                    plot(1, Bstar(1), 'b.', 'MarkerSize', 25)
                else
                    plot(1, Bstar(1), 'bo')
                end
                plot(1, Bstar(2), 'bo')
            end
            % Plot endemic equilibria
            plot(EE(1, 1), EE(2, 1), 'b.', 'MarkerSize', 25)            
            plot(EE(1, 2), EE(2, 2), 'bo')            
            plot(EE(1, 3), EE(2, 3), 'b.', 'MarkerSize', 25)            


            % Plot trajectories from the two initial conditions
            ha.ColorOrderIndex = 1;
            plot(S1, B1, 'LineWidth', 1.5)
            ha.ColorOrderIndex = 1;
            plot(S2, B2, 'LineWidth', 1.5)
            ha.ColorOrderIndex = 1;
            plot(S3, B3, 'LineWidth', 1.5)
            ha.ColorOrderIndex = 1;
            plot(S4, B4, 'LineWidth', 1.5)

            xlim([0.4 1])
            ylim([0 1])
            xlabel('S')
            ylabel('B')
            title("R_0 = " + sprintf('%.2f', R0) + ", " + desc(iCase) + newline + lbls(iTile))
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
    fName = figFolder + "fig_cases" + iFig + ".png";
    saveas(h, fName);
end