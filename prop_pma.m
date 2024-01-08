
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

if ~exist('spec','var'); spec = []; end
if isempty(spec); spec = 'cpma'; end
if strcmp(spec,'olfert'); spec = 'cpma'; end  % for backward compatibility

% Check if code is in +tfer package.
if isfile('+tfer/load_spec.m'); load = @tfer.load_spec;
else; load = @load_spec; end

prop = load(spec);  % load from config file

% --If Dm not set, assign one now ----------------------------------------%
if ~isfield(prop, 'Dm')
    prop.rho0 = pi*1000/6; % ~524;
    prop.Dm = 3;
    
    %-- For soot --%
    % prop.rho0 = 0.0612;
    % prop.Dm = 2.48;
end


%-- Parameters related to CPMA geometry ----------------------------------%
prop.rc = (prop.r2 + prop.r1)/2;
prop.r_hat = prop.r1 / prop.r2;
prop.del = (prop.r2 - prop.r1)/2; % half gap width

prop.A = pi*(prop.r2^2 - prop.r1^2); % cross sectional area of APM
prop.v_bar = prop.Q/prop.A; % average flow velocity


%-- For diffusion --------------------------------------------------------%
kB = 1.3806488e-23; % Boltzmann's constant
prop.D = @(B) kB.*prop.T.*B; % diffusion coefficient


% Fill mass-mobility relation equivalents.
% First, add mat-tfer-pma package to MATLAB path.
fd = fileparts(mfilename('fullpath'));
addpath([fd, filesep, 'tfer-pma']);
prop = prop_massmob(prop);  % mass-mobility relationship

end

