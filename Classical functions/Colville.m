function [y] = Colville(x)

x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);

term1 = 100 * (x1^2-x2)^2;
term2 = (x1-1)^2;
term3 = (x3-1)^2;
term4 = 90 * (x3^2-x4)^2;
term5 = 10.1 * ((x2-1)^2 + (x4-1)^2);
term6 = 19.8*(x2-1)*(x4-1);

y = term1 + term2 + term3 + term4 + term5 + term6;

end
