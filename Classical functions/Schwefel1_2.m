function y = Schwefel1_2(x)
    for i = 1:length(x)
        y(i) = sum(x(1:i))^2;
    end
    y = sum(y);
end