function result = correlation_int(alpha)

    fun = @(x) x./(exp(x)-1);
    result = 1-4/alpha*(integral(fun,0,alpha)+alpha/2-1);

end