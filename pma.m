
% PMA  Bridging function used to evaluate particle mass analyer (PMA) transfer function.
%  AUTHOR: Timothy Sipkens, 2024-01-05


function [Lambda_ii, prop, f_z, qbar] = pma(sp, m, d, z, prop, opt)

% Add mat-tfer-pma package to MATLAB path.
fd = fileparts(mfilename('fullpath'));
addpath([fd, filesep, 'tfer_pma']);


%-- Parse inputs ---------------------------------------------------------%
if ~exist('opt','var'); opt = []; end

% By default, use Taylor series solution baout rc (Case 1C) without diffusion.
% See Sipkens et al., Aerosol Sci. Technol. (2019) for more information.
if isempty(opt); opt = '1C'; end

% If not given, import default properties of PMA, 
% as selected by prop_pma function.
if ~exist('prop','var'); prop = []; end
if isempty(prop); prop = kernel.prop_pma; end

if ~exist('z','var'); z = []; end
if isempty(z); z = (1:4)'; end
%-------------------------------------------------------------------------%


% Assign tfer_pma function to use.
fun = str2func(['tfer_',opt]); % call relevant function from submodule

% For first charge state.
Lambda_ii = zeros([length(sp), length(m), length(z)]);

% Add additional charge states.
for ii=1:length(z)
    Lambda_ii(:, :, ii) = fun(sp, m' .* 1e-18, d' .* 1e-9, z(ii), prop);
end

end

