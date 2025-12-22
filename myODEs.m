function dydt = myODEs(t, y, par)

% Functiopn to evaluate the ODEs for the behaviour-dependent susceptibily
% version of the model

SN = y(1);
SB = y(2);
B = y(3);

S = SN+SB;

Omega = par.Tau*B^2 + par.Chi*(1-S);

dSNdt = -par.Beta*(1-S)*SN + par.Gamma*(1-B-SN) - Omega*SN + par.Alpha*SB; 
dSBdt = -par.Beta*par.q*(1-S)*SB + par.Gamma*(B-SB) + Omega*SN - par.Alpha*SB; 
dBdt = Omega*(1-B) - par.Alpha*B;

dydt = [dSNdt; dSBdt; dBdt];

