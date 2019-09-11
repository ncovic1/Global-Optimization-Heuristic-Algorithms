function y = CEC2014_functions(x)
    global fun_num;
    y = feval(@cec14_func, x', fun_num);
end