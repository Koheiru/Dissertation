function Y = heaviside_restricted(X)
%HEAVISIDE_RESTRICTED    Step function.
%    HEAVISIDE(X) is 0 for X < 0, 1 for X >= 0.

Y = zeros(size(X),'like',X);
Y(X >= 0) = 1;
Y(isnan(X)) = NaN;
