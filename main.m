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


nCases = length(a_vec);

cols = hsv(5);
h = figure(1);
h.Position = [      680         578        1035         420];
tiledlayout(1, nCases, "TileSpacing", "compact");

for iCase = 1:nCases
    a = a_vec(iCase);
    
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
    if a > 4
        % If a > 4 plot the transcritical bifurcations of BDFE and EE  
        plot(R0(ind), qplus(ind), 'color', cols(1, :), 'LineWidth', 2)
        plot(R0(ind), qminus(ind), 'color', cols(2, :), 'LineWidth', 2)
    else
        % Otherwise plot SNB1
        plot(RSNB1, qv, 'color', cols(3, :), 'LineWidth', 2, 'LineStyle', '--')
    end
    % Plot SNB2
    plot(RSNB2, qv, 'color', cols(4, :), 'LineWidth', 2, 'LineStyle', '--')
    % Plot the transcritical bifurcation of DFE and EE
    xline(1, 'k-', 'LineWidth', 2)
    xlim([0, max(R0)])
    ylim([0 1])
    xlabel('R_0')
    ylabel('q')
    % Add region labels
    if iCase == 1
        text(0.4, 0.5, 'NDFE')
        text(1.3, 0.3, 'EE low B')
        text(1.7, 0.7, 'EE+EE')
        text(2.05, 0.9, 'EE high B')
    else
        text(0.2 ,0.5, 'NDFE+BDFE')
        text(1.15, 0.9, 'EE+EE')
        text(1.9, 0.7, 'EE high B')
        text(1.16, 0.48, 'EE+BDFE')
        text(1.9, 0.15, 'BDFE')
        text(1.04, 0.15, 'EE+')
        text(1.04, 0.1, 'BDFE')
    end
    title("\tau/\alpha = " + sprintf('%.1f', a) + ", \chi/\alpha = " + sprintf('%.1f', b) )

end

fName = figFolder + "bifurcation.png";
saveas(h, fName);
