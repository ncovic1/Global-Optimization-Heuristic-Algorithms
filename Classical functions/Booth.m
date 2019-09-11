function [y] = Booth(x)

x1 = x(1);
x2 = x(2);

term1 = (x1 + 2*x2 - 7)^2;
term2 = (2*x1 + x2 - 5)^2;

y = term1 + term2;

end
