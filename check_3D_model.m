clear 
close all

% Figure folder
figFolder = "figures/";



% Number of parameter combinations 
caseIDs = [2, 7];
nCases = length(caseIDs);

% Initial conditions for B and S (all will be plotted in phase plane but only first two will be plotted as time series)
b0_vec = 0:0.05:1;
s0_vec = 0:0.05:0.95;
nb0 = length(b0_vec);
ns0 = length(s0_vec);

[S0, B0] = meshgrid(s0_vec, b0_vec);

% Number of IC types 
nICs = 3;

% Time range (scaled according to the current value of Alpha) 
tSpan = 0:1:10000;

% Labels for subplots
lbls = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)" ];


h = figure(1); 
h.Position = [  54   532   906   442];
tiledlayout(nCases, nICs);

% Run and plot each parameter combination in turn
for iCase = 1:nCases
    % Get parameter values for this case
    par = getPar(caseIDs(iCase));
    %par.Alpha = par.Alpha/10;

    [Sstar, Bstar] = deal(nan(ns0, nb0, nICs));

    for is0 = 1:ns0
        for ib0 = 1:nb0
            s0 = s0_vec(is0);
            b0 = b0_vec(ib0);
            for iIC = 1:nICs    
    
                if iIC == 1
                    sb0 = s0*b0;
                elseif iIC == 2
                    sb0 = min(s0, b0);
                else
                    sb0 = max(0, s0+b0-1);
                end
                IC = [s0-sb0, sb0, b0];

                % Solve model
                [t, Y] = ode45(@(t, y)fullModelODEs(t, y, par), tSpan, IC);
                Sstar(is0, ib0, iIC) = Y(end, 1)+Y(end, 2);
                Bstar(is0, ib0, iIC) = Y(end, 3);
            end
        end
    end

    Istar = 1-Sstar;

    for iIC = 1:nICs
        nexttile;
        imagesc(s0_vec, b0_vec, Istar(:, :, iIC)');
        colorbar;
        colormap(bone);
        h = gca;
        h.YDir = 'normal';
        xlabel('S(0)')
        ylabel('B(0)')
        xlim([0 1])
        ylim([0 1])
        clim([0 0.5])
        title(lbls(iIC + nICs*(iCase-1)));
    end
    drawnow
end
sgtitle('Equilibrium infection prevalence (I*)')

fName = figFolder + "IC_dependence.png";
saveas(h, fName);
