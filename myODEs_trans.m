function dydt = myODEs_trans(t, y, par)

% Functiopn to evaluate the ODEs for the behaviour-dependent transmission
% version of the model

SN = y(1);
SB = y(2);
B = y(3);

S = SN+SB;

Omega = par.Tau*B^2 + par.Chi*(1-S);

Lambda = par.Beta * (1-B-SN + (1-par.q)*(B-SB));

dSNdt = -Lambda*SN + par.Gamma*(1-B-SN) - Omega*SN + par.Alpha*SB; 
dSBdt = -Lambda*SB + par.Gamma*(B-SB) + Omega*SN - par.Alpha*SB; 
dBdt = Omega*(1-B) - par.Alpha*B;

dydt = [dSNdt; dSBdt; dBdt];

