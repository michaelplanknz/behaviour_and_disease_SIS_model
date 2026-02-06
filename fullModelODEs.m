function dydt = fullModelODEs(t, y, par)

% Function to evaluate the ODEs for S_N, S_B and B in the full model

SN = y(1);
SB = y(2);
B = y(3);

S = SN+SB;

Omega = par.Tau*B^2 + par.Chi*(1-S);

Lambda = par.Beta * (1-B-SN + (1-par.qc)*(B-SB));

dSNdt = -Lambda*SN + par.Gamma*(1-B-SN) - Omega*SN + par.Alpha*SB; 
dSBdt = -(1-par.qs)*Lambda*SB + par.Gamma*(B-SB) + Omega*SN - par.Alpha*SB; 
dBdt = Omega*(1-B) - par.Alpha*B;

dydt = [dSNdt; dSBdt; dBdt];

