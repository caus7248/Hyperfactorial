function [Z,Y] = HZ(s,aa)
% _______________________________________________________________________ %
% Evaluate the complex Hurwitz zeta function HZ(s,a), and its derivative.
% ________________________________Inputs_________________________________ %
% s :  a real or complex scalar
% aa:  a real or complex array
% _______________________________Outputs_________________________________ %
% Z:  The Hurwitz zeta function at (s,aa).
% Y:  The derivative of Z with resepect to s, at (s,aa).
% _______________________________Algorithm_______________________________ %
% We begin with the definition
%                          Z = sum_{n=1}^infty  1 /[ (n+a)^s ],
% and take Y to be the derivative of Z with respect to s.
% The first p terms of the sum are aggregated directly. Then the remainder
% is replaced by one of two asymptotic representations based on the
% Euler-Maclaurin expansion, based on the value of s. The latter series
% is in general divergent, and must be trunated after P terms (see below).
% _____________________________Miscellaneous_____________________________ %
% External dependencies:  gampsi.m, get_params.m
% Built-in Matlab functions:  bernoulli
% _______________________________________________________________________ %

% _____________________________Initialization____________________________ %
a   = aa(:);    % flatten into a vector
Z   = 0*a;
Y   = 0*a;
p   = 10;       % use argument reduction to move x -> x+p
P   = 6;        % number of Bernoulli coefficients in the asymptotics
t   = -s;
[c,d,q] = get_params(t,P);

% _______________________________Evaluation______________________________ %
for k = 1:length(a)
    ak = a(k);
    Z_sum = 0;
    Y_sum = 0;
    while real(ak)<p        % directly sum first few terms
        Z_sum = Z_sum + ak^t;
        Y_sum = Y_sum - ak^t*log(ak);
        ak = ak+1;
    end
    % Add the direct sum to the Euler_Maclaurin expansion
    Z(k) = ak^t/2         - ak^(t+1)/(t+1)   + sum(c.*ak.^q);
    Y(k) =- Z(k).*log(ak) - ak^(1+t)/(t+1)^2 + sum(d.*ak.^q);
    Z(k) = Z(k) + Z_sum;
    Y(k) = Y(k) + Y_sum;
end
Z	= reshape(Z,size(aa));
Y	= reshape(Y,size(aa));

end