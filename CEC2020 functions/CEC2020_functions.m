function y = CEC2020_functions(x)
    global fun_num;
    y = feval(@cec20_func, x', fun_num);
end