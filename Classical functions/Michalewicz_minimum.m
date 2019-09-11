function [ fmin, xmin ] = Michalewicz_minimum( dim )

% michalewicz_min 
% input: problem dimension (starting at 1), output: function value at minimum and
% coordinates of minimum

    fmin = 0;
    xmin = zeros(1,dim);
    
    % compute the minimum for each dimension
    for d = 1:dim

        % compute the location of the peak, which is very close to the minimum
        n = round(0.25*d -0.5);
        fraction = sqrt((2*n+1)/(2*d));

        % if the fraction equals 0.5, the peak is located at the minimum.
        if fraction == 0.5
            fmin = fmin  + -1;
            xmin(d) = 0.5*pi;
            continue;
        end
        
        % determine the search domain for ternary search
        if (fraction < 0.5)
            x0 = fraction*pi;
            x3 = 0.5*pi; 
        else
            x0 = 0.5*pi; 
            x3 = fraction*pi;
        end
        
        % ternary search
        while(abs(x3-x0) > 1e-14)
            
            x1 = x0+(x3-x0)/3;
            f1 = michalewicz_term(x1,d);

            x2 = x3-(x3-x0)/3;
            f2 = michalewicz_term(x2,d);
            
            % update the search range
            if( f2 < f1 )
              x0 = x1;
            else
              x3 = x2;
            end
            
        end
        
        % update the values of the minimum
        xmin(d) = (x3+x0)/2;
        fmin = fmin + michalewicz_term(xmin(d),d);

        
    end
    
    

end

function f = michalewicz_term( x, dim ) 
% a single term of the michalewicz function.

    f = -sin(x).*(sin(dim*x.^2/pi).^20);

end
    
    