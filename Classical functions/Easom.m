function z = Easom(x)
    d = length(x);
    z1 = 0;
    for i = 1:d
        z1 = z1 + (x(i)-pi)^2;
    end
    z = 1;
    for i = 1:d
        z = z*cos(x(i))*exp(-z1);
    end
    z = z*(-1)^(d+1);
end
