function [y] = Michalewicz(x)

m = 10;

d = length(x);
sum = 0;

for ii = 1:d
	xi = x(ii);
	new = sin(xi) * (sin(ii*xi^2/pi))^(2*m);
	sum  = sum + new;
end

y = -sum;

end
