function Fano = fano_fun(kon, koff, r, T)

% formula from choubey 2015
Fano = 1 + (2.*r.*koff)./(kon+koff).^2 + (2.*r.*koff)./(kon + koff).^3 .* (exp(-T.*(kon+koff))-1)./T;