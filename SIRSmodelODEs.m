function dydt = SIRSmodelODEs(t, y, par)

% Function to evaluate the ODEs for S_N, S_B, I_N, I_B and B in the SIRS version of the model

SN = y(1);
SB = y(2);
IN = y(3);
IB = y(4);
B = y(5);

S = SN+SB;
I = IN+IB;
RN = 1-B-SN-IN;
RB = B-SB-IB;

Omega = par.Tau*B^2 + par.Chi*I;

Lambda = par.Beta * (IN + (1-par.qc)*IB);

dSNdt = -Lambda*SN             + par.w*RN      - Omega*SN + par.Alpha*SB; 
dSBdt = -(1-par.qs)*Lambda*SB  + par.w*RB      + Omega*SN - par.Alpha*SB; 
dINdt = Lambda*SN              - par.Gamma*IN  - Omega*IN + par.Alpha*IB; 
dIBdt = (1-par.qs)*Lambda*SB   - par.Gamma*IB  + Omega*IN - par.Alpha*IB; 
dBdt = Omega*(1-B) - par.Alpha*B;

dydt = [dSNdt; dSBdt; dINdt; dIBdt; dBdt];

