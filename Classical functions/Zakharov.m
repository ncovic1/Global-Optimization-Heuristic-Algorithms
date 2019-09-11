function [y] = Zakharov(x)


d = length(x);
sum1 = 0;
sum2 = 0;

for ii = 1:d
	xi = x(ii);
	sum1 = sum1 + xi^2;
	sum2 = sum2 + 0.5*ii*xi;
end

y = sum1 + sum2^2 + sum2^4;

end
