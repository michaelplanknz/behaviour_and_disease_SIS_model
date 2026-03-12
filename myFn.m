function F = myFn(x, q, Tau, Chi)

R0 = x(1);
B = x(2);

F(1) = -B + (1-B) * (Tau*B^2 + Chi*(1 - 1/(R0*(1-q*B)) ) );
F(2) = -(1-2*q*B) + Tau*(2*B - 3*B^2 - 3*q*B^2 + 4*q*B^3) + Chi * (-q+2*q*B + 1/R0 - 1);




end