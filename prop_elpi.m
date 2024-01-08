
% PROP_ELPI  Get properties of impactors, such as the ELPI. 
%  AUTHOR: Timothy Sipkens, 2024-01-08.

function prop = prop_elpi(spec)

if ~exist('spec', 'var'); spec = []; end
if isempty(spec); spec = 'default'; end

switch spec
    case {'default', 'elpi'}
        % Values from Järvinen et al. (2014).
        prop.d50 = [15.7, 30.4, 54.1, 94.3, 154., ...
               254., 380., 600., 943., 1620, ...
               2460, 3640, 5340];
        prop.s50 = [3.32, 3.65, 3.89, 3.05, 3.62, ...
               6.30, 8.43, 7.16, 6.21, 5.32, ...
               5.33, 4.14, 3.66];
end

end
