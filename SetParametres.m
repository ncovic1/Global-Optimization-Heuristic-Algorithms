function [f_real, constraints, fun] = SetParametres(fun_number, n)  
    global x_real;
    if fun_number == 1
        fun = @Step;
        x_real = -0.5*ones(1, n);
        f_real = fun(x_real(1,:));
        constraints = 5.12*[-ones(n,1), ones(n,1)];
    elseif fun_number == 2
        fun = @Sphere;
        x_real = zeros(1, n);
        f_real = fun(x_real(1,:)); 
        constraints = 100*[-ones(n,1), ones(n,1)]; 
    elseif fun_number == 3
        fun = @SumSquares;
        x_real = zeros(1, n);
        f_real = fun(x_real(1,:)); 
        constraints = 10*[-ones(n,1), ones(n,1)]; 
    elseif fun_number == 4
        fun = @Quartic;
        x_real = zeros(1, n);
        f_real = fun(x_real(1,:)); 
        constraints = 1.28*[-ones(n,1), ones(n,1)];
    elseif fun_number == 5
        fun = @Beale;
        x_real = [3, 0.5];
        f_real = fun(x_real(1,:));
        constraints = 4.5*[-ones(n,1), ones(n,1)];
    elseif fun_number == 6
        fun = @Easom;
        x_real = pi*ones(1, n);
        f_real = fun(x_real(1,:));
        constraints = 100*[-ones(n,1), ones(n,1)];
    elseif fun_number == 7
        fun = @Matyas;
        x_real = zeros(1, n);
        f_real = fun(x_real(1,:));
        constraints = 10*[-ones(n,1), ones(n,1)];
    elseif fun_number == 8
        fun = @Colville;
        x_real = ones(1, n);
        f_real = fun(x_real(1,:));
        constraints = 10*[-ones(n,1), ones(n,1)];
    elseif fun_number == 9
        fun = @Zakharov;
        x_real = zeros(1, n);
        f_real = fun(x_real(1,:));
        constraints = [-5*ones(n,1), 10*ones(n,1)];
    elseif fun_number == 10
        fun = @Schwefel2_22;
        x_real = zeros(1, n);
        f_real = fun(x_real(1,:));
        constraints = 10*[-ones(n,1), ones(n,1)];
    elseif fun_number == 11
        fun = @Schwefel1_2;
        x_real = zeros(1, n);
        f_real = fun(x_real(1,:));
        constraints = 100*[-ones(n,1), ones(n,1)];
    elseif fun_number == 12
        fun = @DixonPrice;
        i = 1:n;
        x_real(i) = 2.^(-(2.^i-2)./2.^i);
        f_real = 0;
        constraints = 10*[-ones(n,1), ones(n,1)];
    elseif fun_number == 13
        fun = @Bohachevsky1;
        x_real = zeros(1, n);
        f_real = fun(x_real(1,:));
        constraints = 100*[-ones(n,1), ones(n,1)];
    elseif fun_number == 14
        fun = @Booth;
        x_real = [1 3];
        f_real = fun(x_real(1,:));
        constraints = 10*[-ones(n,1), ones(n,1)];
    elseif fun_number == 15
        fun = @HolderTable;
        x_real = [8.055023472141116, 9.664590028909654; 8.055023472141116, -9.664590028909654; 
            -8.055023472141116, 9.664590028909654; -8.055023472141116, -9.664590028909654];
        f_real = fun(x_real(1,:));
        constraints = 10*[-ones(n,1), ones(n,1)];
    elseif fun_number == 16
        fun = @Michalewicz;
        [f_real, x_real] = Michalewicz_minimum(2);
        constraints = [zeros(n,1), pi*ones(n,1)];
    elseif fun_number == 17
        fun = @Michalewicz;
        [f_real, x_real] = Michalewicz_minimum(5);
        constraints = [zeros(n,1), pi*ones(n,1)];
    elseif fun_number == 18
        fun = @Michalewicz;
        [f_real, x_real] = Michalewicz_minimum(10);
        constraints = [zeros(n,1), pi*ones(n,1)];
    elseif fun_number == 19
        fun = @Rastrigin;
        x_real = zeros(1, n);
        f_real = fun(x_real(1,:));
        constraints = 5.12*[-ones(n,1), ones(n,1)];
    elseif fun_number == 20
        fun = @Schaffer2; 
        x_real = [0 0];
        f_real = fun(x_real(1,:));
        constraints = 100*[-ones(n,1), ones(n,1)];
    elseif fun_number == 21
        fun = @Schaffer4; 
        x_real = [0, 1.25313; 0, -1.25313; 1.25313, 0; -1.25313, 0];
        f_real = fun(x_real(1,:));
        constraints = 100*[-ones(n,1), ones(n,1)];
    elseif fun_number == 22
        fun = @Schaffer6; 
        x_real = zeros(1, n);
        f_real = fun(x_real(1,:));
        constraints = 100*[-ones(n,1), ones(n,1)];
    elseif fun_number == 23
        fun = @SixHumpCamelBack; 
        x_real = [0.0898, -0.7126; -0.0898, 0.7126];
        f_real = -1.031628453489836; 
        constraints = 5*[-ones(n,1), ones(n,1)];
    elseif fun_number == 24
        fun = @Bohachevsky2;
        x_real = zeros(1, n);
        f_real = fun(x_real(1,:));
        constraints = 100*[-ones(n,1), ones(n,1)];
    elseif fun_number == 25
        fun = @Bohachevsky3;
        x_real = zeros(1, n);
        f_real = fun(x_real(1,:));
        constraints = 100*[-ones(n,1), ones(n,1)];
    elseif fun_number == 26
        fun = @Shubert;
        x_real = [-7.708303 -0.800343]; % one of 18 global optima
        f_real = -186.7309;
        constraints = 10*[-ones(n,1), ones(n,1)];
    elseif fun_number == 27
        fun = @DropWave;
        x_real = zeros(1, n);
        f_real = fun(x_real(1,:));
        constraints = 5.12*[-ones(n,1), ones(n,1)];
    elseif fun_number == 28
        fun = @Rosenbrock;
        x_real = ones(1, n);
        f_real = fun(x_real(1,:));
        constraints = 30*[-ones(n,1), ones(n,1)];
    elseif fun_number == 29
        fun = @Griewank;
        x_real = zeros(1, n);
        f_real = fun(x_real(1,:));
        constraints = 600*[-ones(n,1), ones(n,1)];  
    elseif fun_number == 30
        fun = @Ackley;
        x_real = zeros(1, n);
        f_real = fun(x_real(1,:));
        constraints = 32*[-ones(n,1), ones(n,1)];
    end   
       
end