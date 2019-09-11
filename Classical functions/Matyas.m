function [y] = Matyas(x)

x1 = x(1);
x2 = x(2);

term1 = 0.26 * (x1^2 + x2^2);
term2 = -0.48*x1*x2;

y = term1 + term2;

end
