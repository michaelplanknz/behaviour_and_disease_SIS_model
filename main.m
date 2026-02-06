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
    

% Run and plot each parameter combination in turn
for iCase = 1:nCases
    % Get parameter values for this case
    par = getPar(iCase);

    % Set up parameter structure for the susceptibility-modulated model
    parSusMod = par;
    parSusMod.qs = par.qc;
    parSusMod.qc = 0;



    h = figure(figNum(iCase));
   
    
    % Get phase plane results
    results = getPhasePlotOutputs(par, dx, dxquiv, dyquiv, pQuiv, tManifold, pert, opts);

    % Run the four initial conditions
    nICs = length(B0);
    for iIC = 1:nICs
        trajTransMod(iIC) = solveODE(@(t, y)fullModelODEs(t, y, par), tSpan, IC(:, iIC));
        trajSusMod(iIC) = solveODE(@(t, y)fullModelODEs(t, y, parSusMod), tSpan, IC(:, iIC));
    end

    iTile = tileNum(iCase);
    makePlots(iTile, results, trajTransMod, trajSusMod, par, desc(iCase))

    

end



for iFig = 1:3
     % Save figure file
    h = figure(iFig);
    fName = figFolder + "fig_cases" + iFig + ".png";
    saveas(h, fName);
end