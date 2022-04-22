function IN = intrinsic_noise_fun(kon, koff, r, T)

% formula from choubey 2015
IN = r.*kon./(kon+koff).*T.*(1 + (2.*r.*koff)./(kon+koff).^2 + (2.*r.*koff)./(kon + koff).^3 .* (exp(-T.*(kon+koff))-1)./T);