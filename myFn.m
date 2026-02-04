function F = myFn(x, q, a, b)

R0 = x(1);
B = x(2);

F(1) = -B + (1-B) * (a*B^2 + b*(1 - 1/(R0*(1-q*B)) ) );
F(2) = -(1-2*q*B) + a*(2*B - 3*B^2 - 3*q*B^2 + 4*q*B^3) + b * (-q+2*q*B + 1/R0 - 1);




end