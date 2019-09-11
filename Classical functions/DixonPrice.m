function [y] = DixonPrice(x)

x1 = x(1);
d = length(x);
term1 = (x1-1)^2;

sum = 0;
for ii = 2:d
	xi = x(ii);
	xold = x(ii-1);
	new = ii * (2*xi^2 - xold)^2;
	sum = sum + new;
end

y = term1 + sum;

end
