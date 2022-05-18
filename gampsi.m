function [G,S] = gampsi(t)
% _______________________________________________________________________ %
% Gamma function G = gamma(t), and its logarithmic deriative S = psi(t)
% _________________________________Inputs________________________________ %
% t:  a real or complex scalar or array
% ________________________________Outputs________________________________ %
% G:  gamma function
% S:  di-gamma function
% _______________________________Algorithm_______________________________ %
% We work with L(t) = log(G(t)), and note that S(t) = L'(t).
% The algorithm has two major steps:
% 1) use the recursion L(t) = L(t+1) - log(t), until real(t+M) >= 10
% 2) evaluate L(t+M) and its derivative using asymptotic series
% _______________________________________________________________________ %

% _____________________________Initialization____________________________ %
L   = 0*t;
S   = 0*t;
b   = [1/6 -1/30 1/42 -1/30 5/66 -691/2730];    % Bernoulli numbers
N   = length(b);
a   = 2*(1:N);
c   = b./( a.*(a-1) );                          % coefficients for L(t)
d   =-b./( a );                                 % coefficients for S(t)

% ___________________________Recursion:  t->t+M__________________________ %
while min(real(t))<10  % note: M is implicit. e.g. if t = -4.5+2i, M = 14.
    L = L-log(t);
    S = S-1./t;
    t = t+1;
end

% __________________________Evaluation:  L(t+M)__________________________ %
L   = L+(t-1/2).*log(t)-t+log(2*pi)/2;          % dominant terms
S   = S+log(t)-1./(2*t);                        % dominant terms
for n = 1:N
    L = L+c(n)./t.^(2*n-1);                     % subdominant terms
    S = S+d(n)./t.^(2*n  );                     % subdominant terms
end
G   = exp(L);
end