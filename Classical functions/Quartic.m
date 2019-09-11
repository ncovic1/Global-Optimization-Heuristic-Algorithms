function y = Quartic(x)
    j = 1:length(x);
    y = sum(j.*x.^4);
end