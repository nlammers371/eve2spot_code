function mean_rate = mean_rate_fun(kon, koff, r, T)

% formula from choubey 2015
mean_rate = r.*kon./(kon+koff).*T;