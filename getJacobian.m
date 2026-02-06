function J = getJacobian(x ,par)

% Calculate Jacobian matrix for the 2D model evaluate at x = [S; B] 

J = zeros(2, 2);

S = x(1);
B = x(2);

J(1, 1) = - par.Beta * (1-2*S) * (1-par.qc*B) - par.Gamma;
J(1, 2) = par.Beta*par.qc * S*(1-S);
J(2, 1) = -par.Chi*(1-B);
J(2, 2) = 2*par.Tau*B*(1-B) - (par.Tau*B^2 + par.Chi*(1-S)) - par.Alpha;


