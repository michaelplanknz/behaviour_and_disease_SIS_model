clear
close all

% Draw two-parameter bifurcation plots for SIS with behaviour model

% Figure folder
figFolder = "figures/";

% Reparameterise as a = tau/alpha  (BDFEs exist when a > 4)
%                   b = chi/alpha
a_vec = [3.8, 5];              
b = 0.2;

% Range of R0 to plot
R0 = 0.5:0.01:2.5;

% Number of values of tau/alpha being plotted
nCases = length(a_vec);

% Number of parameter combinations to add as points for each case
nPoints = [3, 5];


cols = hsv(5);
lbls = ["(a)", "(b)"];

h = figure(1);
h.Position = [      680         578        1035         420];
tiledlayout(1, nCases, "TileSpacing", "compact");

for iCase = 1:nCases
    a = a_vec(iCase);

    % Calculate parameter points to be plotted
    firstPoint = [1, 1+cumsum(nPoints(1:end-1))];
    [R0point, qpoint] = deal(zeros(1, nPoints(iCase)));
    for iPoint = 1:nPoints(iCase)
        par = getPar(firstPoint(iCase) + iPoint-1 );
        R0point(iPoint) = par.Beta/par.Gamma;
        qpoint(iPoint) = par.q;
    end
    
    % BDFEs
    Bstar = 0.5 * (1 + [1, -1]*sqrt(1-4/a));
    
    % Bifurcation values of q at each of the BDFEs
    qplus = (1./R0 + Bstar(1) - 1)./Bstar(1);
    qminus = (1./R0 + Bstar(2) - 1)./Bstar(2);
    
    % Find the value of R0 at the SNBs as a function of q numerically
    % Start at q=1 and work down (as it's easier to guess SNB location when q=1
    % as the S nullcline is vertical)
    qv = 1:-0.01:0;
    
    % First doe SNB2, which is the SNB with the lower value of B (and higher
    % value of R0), which exists for all b
    RSNB2 = zeros(size(qv));
    BSNB2 = zeros(size(qv));
    SSNB2 = zeros(size(qv));
    % Initial guess for R0 and B at SNB2 when q=1
    x0 = [1.5, 0.1];
    for ii = 1:length(qv)
        % Solve the system of two equtions to find the R0 and B where the nullclines intersect tangentially
        x = fsolve(@(x)myFn(x, qv(ii), a, b), x0);
        RSNB2(ii) = x(1);
        BSNB2(ii) = x(2);
        % Calculate the corresponding value of S
        SSNB2(ii) = 1/(RSNB2(ii)*(1+(qv(ii)-1)*BSNB2(ii)) );
        % Use this solution as the initial guess for the next value of q
        x0 = x;
    end
    
    % Look for SNB1 if a < 4
    if a < 4
        RSNB1 = zeros(size(qv));
        BSNB1 = zeros(size(qv));
        SSNB1 = zeros(size(qv));
        % Initial guess for R0 and B at SNB1 when q=1
        x0 = [1, 0.5];
        for ii = 1:length(qv)
            % Solve the system of two equtions to find the R0 and B where the nullclines intersect tangentially
            x = fsolve(@(x)myFn(x, qv(ii), a, b), x0);
            RSNB1(ii) = x(1);
            BSNB1(ii) = x(2);
            % Calculate the corresponding value of S
            SSNB1(ii) = 1/(RSNB1(ii)*(1+(qv(ii)-1)*BSNB1(ii)) );
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
    if a > 4
        % If a > 4 plot the transcritical bifurcations of BDFE and EE (and an empty SNB1 curve for the legend) 
        plot(R0(ind), qminus(ind), 'color', cols(2, :), 'LineWidth', 2, 'DisplayName', 'TC1')
        plot(R0(ind), qplus(ind), 'color', cols(1, :), 'LineWidth', 2, 'DisplayName', 'TC2')
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
    ylabel('q')




    % Add region labels
    if iCase == 1
        % Region 1
        text(0.35, 0.5, '1. NDFE')
        % Region 2
        text(1.1, 0.5, '2. EE low B')
        % Region 3
        text(1.55, 0.8, '3. EE+EE')
        % Region 4
        text(2.1, 0.9, '4. EE')
        text(2.15, 0.85, 'high B')
    else
        % Region 5
        text(0.2 ,0.5, '5. NDFE+BDFE')
        % Region 6
        text(1.03, 0.15, '6. EE+')
        text(1.05, 0.1, 'BDFE')
        % Region 7
        text(1.2, 0.48, '7. EE+')
        text(1.25, 0.43, 'BDFE')
        % Region 8
        text(1.9, 0.15, '8. BDFE')
        % Region 9
        text(1.15, 0.9, '9. EE+')
        text(1.2, 0.85, 'EE')
        % Region 10
        text(1.8, 0.9, '10. EE high B')


    end
    title(lbls(iCase) + " \tau/\alpha = " + sprintf('%.1f', a))

end

fName = figFolder + "bifurcation.png";
saveas(h, fName);
