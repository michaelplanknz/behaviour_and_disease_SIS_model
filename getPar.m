function par = getPar(iCase)

% Set ODE model parameters for each case

par.Gamma = 0.1;        % 0.2
par.Alpha = 0.025;             
par.Chi = 0.2 * par.Alpha;


if iCase == 1
    % Regime 2
    par.Tau = 3.8 * par.Alpha;
    par.Beta = 1.25 * par.Gamma;
    par.q = 0.3;
elseif iCase == 2
    % Regime 3
    par.Tau = 3.8 * par.Alpha;
    par.Beta = 1.75 * par.Gamma;
    par.q = 0.3;
elseif iCase == 3
    % Regime 4
    par.Tau = 3.8 * par.Alpha;
    par.Beta = 2.5 * par.Gamma;
    par.q = 0.3;
elseif iCase == 4
    % Regime 6
    par.Tau = 5 * par.Alpha;
    par.Beta = 1.4 * par.Gamma;
    par.q = 0.3;
elseif iCase == 5
    % Regime 7
    par.Tau = 5 * par.Alpha;
    par.Beta = 1.9 * par.Gamma;
    par.q = 0.3;
elseif iCase == 6
    % Regime 8
    par.Tau = 5 * par.Alpha;
    par.Beta = 1.15 * par.Gamma;
    par.q = 0.7;
elseif iCase == 7
    % Regime 9
    par.Tau = 5 * par.Alpha;
    par.Beta = 1.4 * par.Gamma;
    par.q = 0.7;
else
    % Regime 10
    par.Tau = 5 * par.Alpha;
    par.Beta = 1.9 * par.Gamma; 
    par.q = 0.7;
end



