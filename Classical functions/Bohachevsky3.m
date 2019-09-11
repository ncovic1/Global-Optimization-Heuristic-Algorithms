function [y] = Bohachevsky3(x)

x1 = x(1);
x2 = x(2);

term1 = x1^2;
term2 = 2*x2^2;
term3 = -0.3 * cos(3*pi*x1 + 4*pi*x2);

y = term1 + term2 + term3 + 0.3;

end
