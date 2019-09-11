function [y] = Beale(x)

x1 = x(1);
x2 = x(2);

term1 = (1.5 - x1 + x1*x2)^2;
term2 = (2.25 - x1 + x1*x2^2)^2;
term3 = (2.625 - x1 + x1*x2^3)^2;

y = term1 + term2 + term3;

end
