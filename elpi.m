
% ELPI Compute kernel/transfer function for ELPI+.
%  Uses the fits from Järvinen et al. (2014).
%  
%  K = elpi(DA) computes the kernel, K, for the aerodynamic
%  dimaeters DA in nanometers. 
%  
%  AUTHOR: Timothy Sipkens, 2021-04-05

function [K, d50] = elpi(da, prop)

if ~exist('prop', 'var'); prop = []; end
if isempty(prop); prop = prop_elpi(); end

d50 = prop.d50;

% If data is transposed from expected.
if size(da, 2) == 1; da = da'; end

En = @(da, d50, s50) 1 ./ (1 + (d50 ./ da) .^ (2 .* s50));


% Generate kernel.
K = zeros(length(prop.d50), length(da));
B = zeros(1, length(da));
for ii=length(prop.d50):-1:1
    t0 = En(da, prop.d50(ii), prop.s50(ii));
    
    K(ii, :) = t0 - B;
    
    B = 1 - (1 - B) .* (1 - t0);
end


end

