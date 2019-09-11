function [y] = Shubert(x)

x1 = x(1);
x2 = x(2);
sum1 = 0;
sum2 = 0;

for ii = 1:5
	new1 = ii * cos((ii+1)*x1+ii);
	new2 = ii * cos((ii+1)*x2+ii);
	sum1 = sum1 + new1;
	sum2 = sum2 + new2;
end

y = sum1 * sum2;

end
