
% PROP_PMA  Generates the prop struct used to summarize CPMA parameters.
%  
%  PROP = kernel.prop_pma() creates a default PMA properties structure for 
%  use in evaluating transfer function. This is equiavlent to 
%  prop_pma('cpma').
%  
%  PROP = kernel.prop_pma(SPEC) add a string specifying parameter set. 
%  
%  AUTHOR: Timothy Sipkens, 2019-06-26

function [prop] = prop_pma(spec)

if ~exist('spec', 'var'); spec = []; end
if isempty(spec); spec = 'cpma'; end
if strcmp(spec, 'olfert'); spec = 'cpma'; end  % for backward compatibility

% Load from config file.
prop = load_spec(spec);

% --If Dm not set, assign one now ----------------------------------------%
if ~isfield(prop, 'Dm')
    % Soot.
    % prop.rho0 = 0.0612;
    % prop.Dm = 2.48;
    
    % Water spheres. 
    prop.rho0 = pi*1000/6; % ~524;
    prop.Dm = 3;
end


%-- Parameters related to CPMA geometry ----------------------------------%
prop.rc = (prop.r2 + prop.r1)/2;
prop.r_hat = prop.r1 / prop.r2;
prop.del = (prop.r2 - prop.r1)/2; % half gap width

prop.A = pi*(prop.r2^2 - prop.r1^2); % cross sectional area of APM
prop.v_bar = prop.Q/prop.A; % average flow velocity


%-- For diffusion --------------------------------------------------------%
kB = 1.3806488e-23; % Boltzmann's constant
prop.D = @(B) kB .* prop.T .* B; % diffusion coefficient


% Fill mass-mobility relation equivalents.
% First, add autils package, which should be a folder in upper directory. 
addpath('autils');
prop = massmob.init(prop);  % mass-mobility relationship

end

