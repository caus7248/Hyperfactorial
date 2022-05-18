function [H,L] = HF(m,z,opt)
% _______________________________________________________________________ %
% Compute the hyperfactorial HF(m,z) for complex arguments m, and z.
% ________________________________Inputs_________________________________ %
% m  :  a real or complex scalar/array for the order
% z  :  a real or complex scalar/array for the argument
% opt:  1 for vectorized z, 2 for vectorized m, 3 for Matlab built-in
% _______________________________Outputs_________________________________ %
% H:  The hyperfactorial at (m,z)
% L:  The logarithm of H
% _______________________________Algorithm_______________________________ %
% When z=n+1 is an integer, H is defined by a product
%                          H(m,n+1) = prod_{k=1}^n k^{k^m},
% which reduces to the factorial when m = 0.
% For complex m, we use the continuation
%                         L(m,z) = Y(-m,z) - Y(-m,1),
% where Y(s,z) is the derivative of the Hurwitz zeta function.
% _____________________________Miscellaneous_____________________________ %
% Matlab has only recently added an internal Hurwitz Zeta function. But it
% is very slow. Instead, we use opt = 1 if z is an array, or opt = 2 if m
% is an array.
% _________________________________Usage_________________________________ %
% >>HF(2,4)
% 
% ans = 
%           314928
% _______________________________________________________________________ %
if nargin<3, opt = 1; end
switch opt
    case 1  % use this when z is an array and m is a scalar
        [~,Y1] = HZ(-m,z);
        [~,Y2] = HZ(-m,1);
        L   = Y1-Y2;
    case 2  % use this when m is an array and z is a scalar
        [~,Y1] = HZ2(-m,z);
        [~,Y2] = HZ2(-m,1);
        L   = Y1-Y2;
    case 3  % use the Matlab built-in function, but is very slow
        L   = hurwitzZeta(1,-m,z)-hurwitzZeta(1,-m,1);
end
H   = exp(L);

% __________________Clean up instances of whole numbers__________________ %
if m==round(m)
    wholes = find( z==round(z) );
    H(wholes) = round(H(wholes));
end

end