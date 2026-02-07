clear 
close all

% Figure folder
figFolder = "figures/";

% Number of sets of figures to run, specifying Alpha for each set
nSets = 4;
Alpha = [0.25, 2.5, 0.025, 0.25];

% Number of parameter combinations 
nCases = 8;

% Initial conditions for B and S (all will be plotted in phase plane but only first two will be plotted as time series)
B0 = [0, 0.9, 0.3, 0.6];
S0 = 0.95;

% Time span and perturbation size for finding stable/unstable manifolds
tManifold = [0, 1e3];
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

% For each parameter combo, specify the figure and tile to put the results
% in, and give a short description of parameter regimes being plotted (regimes 2, 3, 4, 6, 7, 8, 9, 10)
figNum = [ones(1, 3), 2*ones(1, 2), 3*ones(1, 3)];
tileNum = [1, 4, 7, 1, 4, 1, 4, 7];
desc = ["EE low B", "bistable EEs", "EE high B", "bistable EEs", "EE high B", "bistable EE/BDFE", "bistable EE/BDFE", "BDFE" ];

for iSet = 1:nSets
    
    % Setup figures
    h = figure(10*(iSet-1)+1);
    h.Position = [ 89   126   931   864];
    tiledlayout(3, 3, "TileSpacing", "compact");
    
    h = figure(10*(iSet-1)+2);
    h.Position = [     164   409   931   562];
    tiledlayout(2, 3, "TileSpacing", "compact");
    
    h = figure(10*(iSet-1)+3);
    h.Position = [  300   126   931   864];
    tiledlayout(3, 3, "TileSpacing", "compact");
    
    
    % Time range (scaled according to the current value of Alpha) 
    tSpan = (0:0.1:100) / (Alpha(iSet)/Alpha(1));

    % Run and plot each parameter combination in turn
    for iCase = 1:nCases
        % Get parameter values for this case
        par = getPar(iCase, Alpha(iSet));
    
        % Set up parameter structure for the susceptibility-modulated model
        parSusMod = par;
        if iSet < 4
            parSusMod.qc = 0;
            parSusMod.qs = par.qc;
        else
        % Split the behaviour effect between transission and susceptibility for
        % a total effect of the same size (par.qc)
            parSusMod.qc = 1 - sqrt(1-par.qc);
            parSusMod.qs = 1 - sqrt(1-par.qc);
        end
            
    
        % Get phase plane results
        results = getPhasePlotOutputs(par, dx, dxquiv, dyquiv, pQuiv, tManifold, pert, opts);
    
        % Run the four initial conditions
        nICs = length(B0);
        for iIC = 1:nICs
            trajTransMod(iIC) = solveODE(@(t, y)fullModelODEs(t, y, par), tSpan, IC(:, iIC));
            trajSusMod(iIC) = solveODE(@(t, y)fullModelODEs(t, y, parSusMod), tSpan, IC(:, iIC));
        end
    
        h = figure(10*(iSet-1)+figNum(iCase));
        iTile = tileNum(iCase);
        makePlots(iTile, results, trajTransMod, trajSusMod, par, desc(iCase))  
    end
end


% Save each figure
for iSet = 1:nSets
    for iFig = 1:3
         % Save figure file
        h = figure(10*(iSet-1)+iFig); 
        % First set is the main figures, label subsequent sets are supplementary
        if iSet == 1
            fName = figFolder + "fig_cases" + iFig + ".png";
        else
            fName = figFolder + "suppfig_set" + string(iSet-1) +"_cases" + string(iFig) + ".png";
        end
        saveas(h, fName);
    end
end