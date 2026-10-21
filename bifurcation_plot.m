clear
close all

% Draw two-parameter bifurcation plots for SIS with behaviour model

% Figure folder
figFolder = "figures/";

% Define values of Tau and Chi
% (BDFEs exist when Tau > 4)
Tau_vec = [3.8, 5];              
Chi = 0.2;

% Range of R0 to plot
R0 = 0.5:0.005:2.5;

% Example parameters for nullcline plots and prevalence vs R0 plots
R0_ex = 2;
qc_ex = [0 0.2 0.4 0.6];

% Time span and perturbation size for finding stable/unstable manifolds
tManifold = [0, 1e3];
pert = 1e-2;

% Root-finding options for finding equilibria 
opts = optimoptions('fsolve', 'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10, 'Display' ,'off', 'Algorithm', 'levenberg-marquardt');

% Spacing for nullcline plotting (dx) and vector field arrow length plotting (dxquiv and dyquiv)
dx = 0.001;
dxquiv = 0.05;
dyquiv = 0.1;

% Power for vector field arrow length plotting (1 for actual length, 0 for
% all the same length)
pQuiv = 0.5; 


% Number of values of tau being plotted
nCases = length(Tau_vec);

% Number of parameter combinations to add as points for each case
nPoints = [3, 5];


cols = hsv(5);
lbls = ["(a)", "(b)"];

h = figure(1);
h.Position = [      680         578        1035         420];
tiledlayout(1, nCases, "TileSpacing", "compact");

for iCase = 1:nCases
    Tau = Tau_vec(iCase);

    % Get (R0,q) parameter points to be plotted
    firstPoint = [1, 1+cumsum(nPoints(1:end-1))];
    [R0point, qpoint] = deal(zeros(1, nPoints(iCase)));
    for iPoint = 1:nPoints(iCase)
        parTmp = getPar(firstPoint(iCase) + iPoint-1);
        R0point(iPoint) = parTmp.Beta/parTmp.Gamma;
        qpoint(iPoint) = parTmp.qc;
    end
    
    
    % BDFEs
    Bstar = 0.5 * (1 + [1, -1]*sqrt(1-4/Tau));
    
    % Bifurcation values of q at each of the BDFEs (transcritical
    % bifurcations)
    qplus = (1 - 1./R0)./Bstar(1);
    qminus = (1 - 1./R0)./Bstar(2);
    
    % Find the value of R0 at the SNBs as a function of q numerically
    % Start at q=0 and work up (as it's easier to guess SNB location when
    % q=0 as the S nullcline is vertical)
    qv = 0:0.01:1;
    
    % First do SNB2, which is the SNB with the lower value of B (and higher
    % value of R0), which exists for all b
    RSNB2 = zeros(size(qv));
    BSNB2 = zeros(size(qv));
    SSNB2 = zeros(size(qv));
    % Initial guess for R0 and B at SNB2 when q=1
    x0 = [1.5, 0.1];
    for ii = 1:length(qv)
        % Solve the system of two equtions to find the R0 and B where the nullclines intersect tangentially
        x = fsolve(@(x)myFn(x, qv(ii), Tau, Chi), x0);
        RSNB2(ii) = x(1);
        BSNB2(ii) = x(2);
        % Calculate the corresponding value of S
        SSNB2(ii) = 1/(RSNB2(ii)*(1-qv(ii)*BSNB2(ii)) );
        % Use this solution as the initial guess for the next value of q
        x0 = x;
    end
    
    % Look for SNB1 if Tau < 4
    if Tau < 4
        RSNB1 = zeros(size(qv));
        BSNB1 = zeros(size(qv));
        SSNB1 = zeros(size(qv));
        % Initial guess for R0 and B at SNB1 when q=1
        x0 = [1, 0.5];
        for ii = 1:length(qv)
            % Solve the system of two equtions to find the R0 and B where the nullclines intersect tangentially
            x = fsolve(@(x)myFn(x, qv(ii), Tau, Chi), x0);
            RSNB1(ii) = x(1);
            BSNB1(ii) = x(2);
            % Calculate the corresponding value of S
            SSNB1(ii) = 1/(RSNB1(ii)*(1-qv(ii)*BSNB1(ii)) );
            % Use this solution as the initial guess for the next value of q
            x0 = x;
        end
    
    end
    
    

    
    % Plot two parameter bifurcation diagram (R0, q)

    % Only need to plot curves when R0 >= 1
    ind = (R0 >= 1);
    
    nexttile;
    hold on
    % Plot the transcritical bifurcation of DFE and EE
    xline(1, 'k-', 'LineWidth', 2, 'DisplayName', 'TC0')
    if Tau > 4
        % If Tau > 4 plot the transcritical bifurcations of BDFE and EE (and an empty SNB1 curve for the legend) 
        plot(R0(ind), qminus(ind), 'color', cols(2, :), 'LineWidth', 2, 'DisplayName', 'TC-')
        plot(R0(ind), qplus(ind), 'color', cols(1, :), 'LineWidth', 2, 'DisplayName', 'TC+')
        plot(nan, nan, 'color', cols(3, :), 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'SNB1')
        lgd = legend('Location', 'eastoutside');
    else
        % Otherwise plot SNB1
        plot(RSNB1, qv, 'color', cols(3, :), 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'SNB1')
    end
    % Plot SNB2
    plot(RSNB2, qv, 'color', cols(4, :), 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'SNB2')
    % Plot parameter point
    plot(R0point, qpoint, 'o', 'color', 'k', 'HandleVisibility', 'off')
    xlim([0, max(R0)])
    ylim([0 1])
    xlabel('R_0')
    ylabel('behavioural effect (q_c)')

    % Add region labels
    if iCase == 1
        % Region 1
        text(0.35, 0.5, '1. NDFE')
        % Region 2
        text(1.1, 0.5, '2. EE low B')
        % Region 3
        text(1.55, 0.2, '3. EE+EE')
        % Region 4
        text(2.1, 0.15, '4. EE')
        text(2.15, 0.1, 'high B')
    else
        % Region 5
        text(0.2 ,0.5, '5. NDFE+BDFE')
        % Region 6
        text(1.03, 0.9, '6. EE+')
        text(1.05, 0.85, 'BDFE')
        % Region 7
        text(1.2, 0.57, '7. EE+')
        text(1.25, 0.52, 'BDFE')
        % Region 8
        text(1.9, 0.85, '8. BDFE')
        % Region 9
        text(1.18, 0.15, '9. EE+')
        text(1.23, 0.1, 'EE')
        % Region 10
        text(1.75, 0.15, '10. EE high B')


    end
    title(lbls(iCase) + " \tau = " + sprintf('%.1f', Tau))

end

fName = figFolder + "bifurcation.png";
saveas(h, fName);






% Plot example nullclines for suppmentary
h = figure(2);
h.Position = [      680         578        1035         420];
tiledlayout(1, nCases, "TileSpacing", "compact");
greyCol =  0.8*[1 1 1];
colOrd = colororder;

for iCase = 1:nCases
    par.Gamma = 1;        
    par.Beta = R0_ex*par.Gamma;
    par.Alpha = 0.25;             
    par.Tau = Tau_vec(iCase);
    par.Chi = Chi;
    par.qc = qc_ex(1);
    results = getPhasePlotOutputs(par, dx, dxquiv, dyquiv, pQuiv, tManifold, pert, opts);

     nexttile;
     hold on

    % Plot B nullcline 
    plot(results.Snull2, results.Bnull2, 'LineWidth', 2, 'Color', colOrd(5, :));

    % Plot S nullcline
    plot(results.Snull1, results.Bnull1, 'LineWidth', 2, 'Color', colOrd(4, :));

    % Set up a vector of color dilution factors
    colAlph = linspace(0, 0.8, length(qc_ex));
    % Plot the S nullcline for increasing values of qc
    for iq = 2:length(qc_ex)
        par.qc = qc_ex(iq);
        results = getPhasePlotOutputs(par, dx, dxquiv, dyquiv, pQuiv, tManifold, pert, opts);
        plot(results.Snull1, results.Bnull1, 'LineWidth', 2, 'Color', (1-colAlph(iq))*colOrd(4, :) + colAlph(iq));
    end

    % Plot vertical branch of S nullcline
    plot([1 1], [0 1], 'k-', 'LineWidth', 2)

   
    xlim([0, 1])
    ylim([0, 1])
    xlabel('S')
    ylabel('B')  

    % S nullcline arrow and label
    h = annotation('arrow');
    h.Parent = gca;
    h.X = 1/R0_ex + [0, -0.1]; 
    h.Y = 0.025* [1 1];
    h.Color = [0 0 0];
    text(0.2, 0.07, "increasing R_0");

    if iCase == 2
        legend([ "B nullcline", "S nullcline (q_c = " + string(qc_ex) + ")", "S nullcline (all q_c)"], 'Location', 'eastoutside')
    end
    title(lbls(iCase) + " \tau = " + sprintf('%.1f', par.Tau))

end

fName = figFolder + "nullclines.png";
saveas(h, fName);




% Solve for equilibria as a function of R0 for each value of tau and for an example value of qc 

SMALL = 1e-12;

casesToPlot = [1 4];
nVals = length(R0);

Istar = nan(4, nVals, nCases);
Rplus = nan(nVals, nCases);
R0TC_minus = nan(1, 2);
R0TC_plus = nan(1, 2);
for iCase = 1:nCases
    par = getPar(casesToPlot(iCase));

    if par.Tau >= 4
        % Calculate stability of BDFE+
        Bminus = 0.5 * (1 - sqrt(1-4/par.Tau));
        Bplus = 0.5 * (1 + sqrt(1-4/par.Tau));
        Rplus(:, iCase) = R0.*(1-par.qc*Bplus);
        R0TC_minus(iCase) = 1/(1-par.qc*Bminus);
        R0TC_plus(iCase) = 1/(1-par.qc*Bplus);
    end
    

    bif1Flag = false;
    bif2Flag = false;
    for iVal = 1:nVals
        coeffs = zeros(1, 5);
        % coeffs(j) is the coefficient of the (5-j)^th power or B in the
        % polynomial to be solved
        coeffs(5) = par.Chi*(1-1/R0(iVal));
        coeffs(4) = - ( 1 + par.Chi*(par.qc + 1 - 1/R0(iVal)) );
        coeffs(3) = par.qc*(1+par.Chi) + par.Tau;
        coeffs(2) = -par.Tau*(1+par.qc);
        coeffs(1) = par.Tau*par.qc;
        xr = roots(coeffs);
        xr = xr(abs(imag(xr)) < SMALL);
        Ir = 1 - 1./(R0(iVal).*(1-par.qc*xr));
        nRealRoots = length(xr);
        if nRealRoots == 4
            bif1Flag = true;
        end
        if bif1Flag && nRealRoots == 2
            bif2Flag = true;
        end

        if ~bif2Flag    
            Istar(1:length(xr), iVal, iCase) =  sort(Ir, 'descend')';
        else
            % After second bifurcation, store roots in positions 3 and 4
            % for continuous plotting of branches
            Istar(3:(2+length(xr)), iVal, iCase) = sort(Ir, 'descend')';
        end
    end
end

Istar(Istar < 0) = nan;

% Plot equilibria versus R0

h = figure(3);
h.Position = [      680         578        1035         420];
tiledlayout(1, nCases, "TileSpacing", "compact");
greyCol =  0.8*[1 1 1];
colOrd = colororder;
zz = zeros(size(R0));
for iCase = 1:nCases
    nexttile;
    hold on
    % Plot dummy data for legend
    plot(nan, nan, "Color", colOrd(1, :), 'LineWidth', 2)
    plot(nan, nan, "Color", colOrd(2, :), 'LineWidth', 2)
    plot(nan, nan, "Color", colOrd(3, :), 'LineWidth', 2)
    plot(nan, nan, "Color", colOrd(4, :), 'LineWidth', 2)
    plot(nan, nan, 'Color', [0 0 0])
    plot(nan, nan, "Color", colOrd(3, :), 'Marker', 'o', 'MarkerFaceColor', colOrd(3, :), 'LineStyle', 'none')
    plot(nan, nan, "Color", colOrd(2, :), 'Marker', 'o', 'MarkerFaceColor', colOrd(2, :), 'LineStyle', 'none')
    plot(nan, nan, "Color", colOrd(3, :), 'Marker', 'square', 'MarkerFaceColor', colOrd(3, :), 'LineStyle', 'none')
    plot(nan, nan, "Color", colOrd(4, :), 'Marker', 'square', 'MarkerFaceColor', colOrd(4, :), 'LineStyle', 'none')


    % Plot equilibrium prevalence on each branch
    plot(R0(Rplus(:, iCase) <= 1), zz(Rplus(:, iCase) <= 1), "Color", colOrd(1, :), 'LineWidth', 2 )
    plot(R0(Rplus(:, iCase) >= 1), zz(Rplus(:, iCase) >= 1), "Color", colOrd(1, :), 'LineWidth', 2, 'LineStyle', '--' )

    plot(R0, Istar(1, :, iCase), 'LineWidth', 2, 'Color', colOrd(2, :))
    plot(R0, Istar(2, :, iCase), 'LineWidth', 2, 'Color', colOrd(3, :), 'LineStyle', '--' )
    plot(R0, Istar(3, :, iCase), 'LineWidth', 2, 'Color', colOrd(4, :))



    % Plot no behaviour equilibrium for comparison
    plot(R0, 1-1./R0, 'Color', [0 0 0])


    % Plot bifurcation points
  
    iSNB1 = find(~isnan(Istar(2, :, iCase)), 1, 'first');
    iSNB2 = find(~isnan(Istar(2, :, iCase)), 1, 'last');
    if Tau_vec(iCase) < 4   
        plot(R0(iSNB1), Istar(2, iSNB1, iCase), 'Marker', 'o', 'color', colOrd(3, :), 'MarkerFaceColor', colOrd(3, :), 'LineStyle', 'none')
    end
    plot(R0(iSNB2), Istar(2, iSNB2, iCase), 'Marker', 'o', 'color', colOrd(2, :), 'MarkerFaceColor', colOrd(2, :), 'LineStyle', 'none')

    plot(R0TC_minus(iCase), 0, 'Marker', 'square', 'color', colOrd(3, :), 'MarkerFaceColor', colOrd(3, :), 'LineStyle', 'none')
    plot(R0TC_plus(iCase), 0, 'Marker', 'square', 'color', colOrd(4, :), 'MarkerFaceColor', colOrd(4, :), 'LineStyle', 'none')

    xlim([1 inf])
    ylim([0, 0.5])
    xlabel('R_0')
    ylabel('equilibrium prevalence (I*)')  

    title(lbls(iCase) + " \tau = " + sprintf('%.1f', Tau_vec(iCase)))

end

legend('BDFE_+', 'EE low B', 'EE med B', 'EE high B', 'No behaviour', 'SNB1', 'SNB2', 'TC-', 'TC+', 'Location', 'southeast')

fName = figFolder + "prevalence_vs_R0.png";
saveas(h, fName);


