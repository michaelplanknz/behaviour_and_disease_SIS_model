function dydt = myODEs_2D(t, y, par)

% Functiopn to evaluate the ODEs for the behaviour-dependent transmission
% version of the model

S = y(1);
B = y(2);


Omega = par.Tau*B^2 + par.Chi*(1-S);

Lambda = par.Beta * (1-S) * (1-B+par.q*B);

dSdt = -Lambda*S + par.Gamma*(1-S); 
dBdt = Omega*(1-B) - par.Alpha*B;

dydt = [dSdt; dBdt];

