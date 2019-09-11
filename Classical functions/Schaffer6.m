function y = Schaffer6(x)
    y = 0;
    for i = 1:length(x)-1
    	y = y + 0.5 + ((sin(sqrt(x(i)^2+x(i+1)^2)))^2-0.5)/((1+0.001*(x(i)^2+x(i+1)^2))^2);
    end
end