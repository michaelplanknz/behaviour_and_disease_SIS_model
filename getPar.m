function par = getPar(iCase)

% Set ODe model parameters for each case

par.Gamma = 0.2;
par.Alpha = 0.025;              % 0.025 for slow, 0.1 for fast
par.Chi = 0.2 * par.Alpha;

if iCase == 1
    par.Tau = 3.8 * par.Alpha;
    par.Beta = 0.25;
    par.q = 0.7;
elseif iCase == 2
    par.Tau = 3.8 * par.Alpha;
    par.Beta = 0.42;
    par.q = 0.7;
elseif iCase == 3
    par.Tau = 3.8 * par.Alpha;
    par.Beta = 0.5;
    par.q = 0.7;
elseif iCase == 4
    par.Tau = 5 * par.Alpha;
    par.Beta = 0.28;
    par.q = 0.7;
elseif iCase == 5
    par.Tau = 5 * par.Alpha;
    par.Beta = 0.4;
    par.q = 0.7;
elseif iCase == 6
    par.Tau = 5 * par.Alpha;
    par.Beta = 0.24;
    par.q = 0.3;
elseif iCase == 7
    par.Tau = 5 * par.Alpha;
    par.Beta = 0.3;
    par.q = 0.3;
else
    par.Tau = 5 * par.Alpha;
    par.Beta = 0.4;
    par.q = 0.3;
end

