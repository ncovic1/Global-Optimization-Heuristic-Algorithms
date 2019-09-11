function [y] = Griewank(x)
    d = length(x);
    sum = 0;
    prod = 1;

    for ii = 1:d
        xi = x(ii);
        sum = sum + xi^2/4000;
        prod = prod * cos(xi/sqrt(ii));
    end

    y = sum - prod + 1;

end
