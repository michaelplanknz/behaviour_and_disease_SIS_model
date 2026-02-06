function [dSdt, dBdt] = myVectorField(S, B, par)

% Functiopn to evaluate the vector field (dSdt and dBdt) for the
% transmission-mopdulated model

Omega = par.Tau*B.^2 + par.Chi*(1-S);

Lambda = par.Beta * (1-S) .* (1-par.qc*B);

dSdt = -Lambda.*S + par.Gamma*(1-S); 
dBdt = Omega.*(1-B) - par.Alpha*B;


