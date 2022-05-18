function [c,d,q] = get_params(t,P)
% _______________________________________________________________________ %
% Parameters for the asymptotic series in the Hurwitz zeta function Z(s,a)
% and its derivative Y(s,a).
% _________________________________Inputs________________________________ %
% t:  a real or complex scalar or array
% P:  number of terms in the asymptotic series
% ________________________________Outputs________________________________ %
% c:  asymptotic coefficients in Z(s,a)
% d:  asymptotic coefficients in Y(s,a)
% q:  asymptotic exponents in Z(s,a) and Y(s,a)
% _______________________________Algorithm_______________________________ %
% We work with L(t) = log(G(t)), and note that S(t) = L'(t).
% The algorithm has two major steps:
% 1) use the recursion L(t) = L(t+1) - log(t), until real(t+M) >= 10
% 2) evaluate L(t+M) and its derivative using asymptotic series
% _______________________________________________________________________ %

m   = fix(real(t));
if m<0  % Real(s) > 0:  use the standard asymptotic series 
    q1  = 2:2:2*P;
    q   = t+1-q1;
    [gq,sq] = gampsi(-q);
    [gt,st] = gampsi(-t);
    c   = bernoulli(q1).*gq./gt./factorial(q1);
    d   = c.*(sq-st);
else    % Real(s) <= 0:  use analytic continuation, reflection formula
    q1  = 2:2:m+1;
    q2  = 2*fix( (m+1)/2 )+(2:2:2*P);
    t1  = t+1-q1;                           % these exponents are positive
    t2  = t+1-q2;                           % these exponents are negative
    q   = [t1 t2];
    [gt,st] = gampsi(t+1);
    [g1,s1] = gampsi(t1+1);
    [g2,s2] = gampsi(-t2);
    c1  = bernoulli(q1)./factorial(q1).*( gt./g1 );
    c2  = bernoulli(q2)./factorial(q2).*( gt.*g2 );     % reflection
    d1  = c1.*(s1-st);
    d2  = c2.*((s2-st)*sin(pi*t)/pi-cos(pi*t));         % reflection
    d   = -[d1 d2];
    c2  = c2*sin(pi*t)/pi;                              % reflection
    c   = -[c1 c2];
end

end