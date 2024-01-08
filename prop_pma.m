
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


% Check if load configuration file.
prop = check_spec(spec);

if ~isstruct(prop)  % otherwise load hard coded examples
    switch spec
        
        %-- CPMA parameters from Olfert lab ----------------------------------%
        case {'olfert','cpma',' CPMA'} % final entry associated with data import
            prop.r1 = 0.06; % inner electrode radius [m]
            prop.r2 = 0.061; % outer electrode radius [m]
            prop.L = 0.2; % length of chamber [m]
            prop.p = 1; % pressure [atm]
            prop.T = 293; % system temperature [K]
            prop.Q = 3/1000/60;%0.3/1000/60;%1.5/1000/60; % volume flow rate (m^3/s) (prev: ~1 lpm)
            prop.omega_hat = 32/33; % ratio of angular speeds
    
        %-- CPMA/APM parameters from Buckley et al. --------------------------%
        case 'buckley'
            prop.r2 = 0.025; % outer electrode radius [m]
            prop.r1 = 0.024; % inner electrode radius [m]
            prop.L = 0.1;    % length of APM [m]
            RPM = 13350; % rotational speed [rpm]
            prop.omega = RPM*2*pi/60; % rotational speed [rad/s]
            prop.omega_hat = 1; % APM, so rotational speed is the same
            prop.Q = 1.02e-3/60; % aerosol flowrate [m^3/s]
            prop.T = 298; % system temperature [K]
            prop.p = 1; % system pressure [atm]
        
        %-- CPMA parameters from Olfert lab ----------------------------------%
        case {'flarenet18','fn18'}
            prop.r1 = 0.06; % inner electrode radius [m]
            prop.r2 = 0.061; % outer electrode radius [m]
            prop.L = 0.2; % length of chamber [m]
            prop.p = 1; % pressure [atm]
            prop.T = 293; % system temperature [K]
            prop.Q = 0.3/1000/60; % volume flow rate (m^3/s) (prev: ~1 lpm)
            prop.omega_hat = 32/33; % ratio of angular speed
            prop.rho0 = pi*1000/6; % ~524;
            prop.Dm = 3;
            
        case {'soot-salt'}
            prop.r1 = 0.06; % inner electrode radius [m]
            prop.r2 = 0.061; % outer electrode radius [m]
            prop.L = 0.2; % length of chamber [m]
            prop.p = 1; % pressure [atm]
            prop.T = 293; % system temperature [K]
            prop.Q = 0.3/1000/60; % volume flow rate (m^3/s) (prev: ~1 lpm)
            prop.omega_hat = 32/33; % ratio of angular speed
            prop.rho0 = 0.0612; % ~524;
            prop.Dm = 2.48;
            
        %-- APM parameters from Ehara et al. -------------%
        case 'ehara'
            prop.r2 = 0.103; % outer electrode radius [m]
            prop.r1 = 0.1; % inner electrode radius [m]
            prop.L = 0.2;    % length of APM [m]
            prop.omega_hat = 1; % APM, so rotational speed is the same
            prop.Q = 0.5/1000/60; % aerosol flowrate [m^3/s], assumed
            prop.T = 298; % system temperature [K]
            prop.p = 1; % system pressure [atm]
        
        %-- Parameters from Olfert and Collings -------------%
        %   Nearly identical to the Ehara et al. case
        case 'olfert-collings'
            prop.r2 = 0.103; % outer electrode radius [m]
            prop.r1 = 0.1; % inner electrode radius [m]
            prop.L = 0.2;
            prop.omega_hat = 0.945;
            prop.Q = 0.5/1000/60; % aerosol flowrate [m^3/s]
            prop.T = 295; % system temperature [K]
            prop.p = 1; % system pressure [atm]
        
        %-- Parameters from Kuwata --------------------------%
        case 'kuwata'
            prop.r2 = 0.052; % outer electrode radius [m]
            prop.r1 = 0.05; % inner electrode radius [m]
            prop.L = 0.25;
            prop.omega_hat = 1;
            prop.Q = 1.67e-5; % aerosol flowrate [m^3/s]
            prop.T = 295; % system temperature [K]
            prop.p = 1; % system pressure [atm]
            
	    case 'santavac'
            prop.r1 = 0.06; % inner electrode radius [m]
            prop.r2 = 0.061; % outer electrode radius [m]
            prop.L = 0.2; % length of chamber [m]
            prop.p = 1; % pressure [atm]
            prop.T = 293; % system temperature [K]
            prop.Q = 0.3/1000/60; % volume flow rate (m^3/s) (prev: ~1 lpm)
            prop.omega_hat = 32/33; % ratio of angular speed
            prop.rho0 = pi*1000/6; % ~524;
            prop.Dm = 3;
    
    end
end


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
addpath tfer_pma;
prop = prop_massmob(prop);

end

