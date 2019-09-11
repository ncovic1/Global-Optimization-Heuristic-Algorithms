function z = Rastrigin(x)
    z = 10*length(x);
    for i=1:length(x)
        z = z + (x(i)^2-10*cos(2*pi*x(i)));
    end
end
