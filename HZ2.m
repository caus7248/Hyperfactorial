function [Z,Y] = HZ2(ss,a)
% _______________________________________________________________________ %
% Evaluate the complex Hurwitz zeta function HZ(s,a), and its derivative.
% ________________________________Inputs_________________________________ %
% ss:  a real or complex array
% a :  a real or complex scalar
% _______________________________Outputs_________________________________ %
% Z:  The Hurwitz zeta function at (ss,a).
% Y:  The derivative of Z with resepect to s, at (ss,a).
% _______________________________Algorithm_______________________________ %
% We begin with the definition
%                          Z = sum_{n=1}^infty  1 /[ (n+a)^s ],
% and take Y to be the derivative of Z with respect to s.
% The first p terms of the sum are aggregated directly. Then the remainder
% is replaced by one of two asymptotic representations based on the
% Euler-Maclaurin expansion, based on the value of s. The latter series
% is in general divergent, and must be trunated after P terms (see below).
% _____________________________Miscellaneous_____________________________ %
% External dependency:  gampsi.m, get_params.m
% Built-in Matlab functions:  bernoulli
% _______________________________________________________________________ %

% _____________________________Initialization____________________________ %
s   = ss(:);    % flatten into a vector
Z   = 0*s;
Y   = 0*s;
Z_sum = 0*s;
Y_sum = 0*s;
p   = 10;
P   = 6;

% _______________________________Evaluation______________________________ %
% First get the direct sum
while real(a)<p
    Z_sum   = Z_sum + a.^(-s);
    Y_sum   = Y_sum - a.^(-s)*log(a);
    a   = a+1;
end
% Next compute the asymptotic series, looping over s
for k = 1:length(s)
    t   = -s(k);
    [c,d,q] = get_params(t,P);
    Z(k) = a^t/2        - a^(t+1)/(t+1)   + sum(c.*a.^q);
    Y(k) =-Z(k)*log(a)  - a^(t+1)/(t+1)^2 + sum(d.*a.^q);
end
Z   = Z+Z_sum;
Y   = Y+Y_sum;
Z	= reshape(Z,size(ss));
Y	= reshape(Y,size(ss));

end