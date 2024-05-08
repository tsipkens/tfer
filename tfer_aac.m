
% TFER_AAC  Transfer function of the AAC. 
%  Based on limited trajectory model. 
%  See also Johnson et al., Aerosol Sci. Technol (2021) and 
%  Tavakoli and Olfert (2013). 
%  
%  LAMBDA = tfer_aac(D_STAR, D, PROP) computes the AAC transfer function at 
%  aerodynamic diameters D [nm] for a specified da septoint, D_STAR [nm],  
%  and AAC properties, PROP. 
%  
%  LAMBDA = tfer_aac(OMEGA, D, PROP, OPTS) instead uses an angular speed,
%  OMEGA, as the input to constrain the setpoint, requiring that 
%  OPTS.input = 'omega'.
%  
%  [LAMBDA,DA_STAR,PROP] = tfer_aac(...) adds outputs for setpoint
%  information, DA_STAR, and additional characteristics of the classifier,
%  PROP.
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2024-01-10

function [Lambda, da_star, prop] = tfer_aac(da_star, da, prop, opts, varargin)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('opts', 'var'); opts = []; end
if isempty(opts); opts = struct(); end
if ~isfield(opts, 'input'); opts.input = 'd_star'; end  % by default, assume input is d_star
if ~isfield(opts, 'scan'); opts.scan = 0; end  % by default, use steady state (not scanning)
if ~isfield(opts, 'model'); opts.model = 'ps'; end  % by default, limited trajectory
if ~isfield(opts, 'ideal'); opts.ideal = false; end
if ~isfield(opts, 'diffusion'); opts.diffusion = true; end  % whether to include diffusion

% Convert from nm to m for calculations.
% Conversion of da_star is below, depending on opts structure. 
da = da .* 1e-9;
%-------------------------------------------------------------------------%

% Computer other classifier properties.
prop.del = (prop.Qs - prop.Qa) ./ (prop.Qs + prop.Qa);
prop.bet = (prop.Qs + prop.Qa) ./ (prop.Qsh + prop.Qexh);

% Condense if not unique values.
if all(prop.del == prop.del(1)); prop.del = prop.del(1); end
if all(prop.bet == prop.bet(1)); prop.bet = prop.bet(1); end

% Compute resolution.
prop.Rt = 1 ./ prop.bet;

% Get some basic quanities. 
mu = 1.82e-5;  % gas viscosity [Pa*s]
rc = prop.r1 ./ prop.r2;  % ratio of radii
tf = pi .* prop.L .* (prop.r2.^2 - prop.r1.^2) ./ ...
    (prop.Qa + prop.Qsh);  %  axial velovity in the classifier
gam = (prop.Qa + prop.Qsh - prop.Qs .* (1 - rc .^ 2)) ./ ...
    (prop.Qa + prop.Qsh - prop.Qsh .* (1 - rc .^ 2));


% Interpret setpoint.
if strcmp(opts.input, 'd_star')  % d_star given, so compute omega
    da_star = da_star .* 1e-9;  % convert input from nm to m
    tau_star = Cc(da_star) .* 1e3 .* da_star.^2 ./ (18 .* mu);
    omega = sqrt((prop.Qsh + prop.Qexh) ./ ...
        (pi .* (prop.r2 + prop.r1).^2 .* prop.L .* tau_star));

elseif strcmp(opts.input, 'omega')  % omega given, so compute d_star
    omega = da_star;
    tau_star = (prop.Qsh + prop.Qexh) ./ ...
        (pi .* (prop.r2 + prop.r1).^2 .* prop.L .* omega.^2);
    da_star = cacl_d_star(tau_star, mu);

end

% Relaxation time for da input.
tau = Cc(da) .* 1e3 .* da.^2 ./ (18 .* mu);

%-- Constants for scanning -----------------------------------------------%
%   prop.omega_s is initial angular speed
%   prop.omega_e is final angular speed
%   prop.tsc is scan time
if opts.scan
    tau_sc = prop.tsc ./ (2 .* log(prop.omega_e ./ prop.omega_s));
    tf = pi .* prop.L .* ...
        (prop.r2.^2 - prop.r1.^2) ./ (prop.Qa + prop.Qsh);  % Johsnon et al., Eq. 13

    c_sc = prop.omega_s.^2 .* tau_sc .* (1 - exp(-tf ./ tau_sc));
end

%-- Continue with transfer function calculation. -------------------------%
switch opts.model

    %-- Limited trajectory approach --------------------------------------%
    case {'lt', 'limited-trajectory'}
        if ~opts.scan  % steady state
            K = omega .^ 2 .* tf;  % for stepping
            
        else  % scanning
            c_tau_star = 1 ./ (2 .* c_sc) .* (log(1./rc) + 1/2 .* log(gam));
            
            if strcmp(opts.input, 'time')  % OPTION 1: time to d_star
                tm = d_star;  % classifier exit time
                tau_star = c_tau_star .* exp(-tm ./ tau_sc);
                cacl_d_star(tau_star, mu);
                
            else  % OPTION2: d_star to time
                tm = -tau_sc .* log(tau_star ./ c_tau_star);
            end
            
            K = c_sc .* exp(tm ./ tau_sc);
        end
        
        %{
        tau_max = 1./ K .* log(prop.r2 ./ prop.r1);
        tau_min = 1 ./ (2 .* K) .* log(gam);
        
        % Setpoint check calculation.
        tau_star2 = (tau_min + tau_max) ./ 2;
        d_star2 = fminsearch(@(d) ...
            1e9 .* norm(Cc(d) .* 1e3 .* d.^2 ./ (18 .* mu) - tau_star2), 100e-9);
        omega2 = sqrt((prop.Qsh + prop.Qexh) ./ ...
             (pi .* (prop.r2 + prop.r1).^2 .* prop.L .* tau_star2));
        %}
        
        f1 = (prop.Qa + prop.Qsh .* rc.^2 - ...
            exp(-2 .* tau .* K) .* (prop.Qa + prop.Qsh - prop.Qs .* (1 - rc .^ 2))) ./ ...
            (prop.Qa .* (1 - rc .^ 2));
        
        f2 = (prop.Qa + prop.Qsh) ./ prop.Qa .* ...
            (exp(-2 .* tau .* K) - rc .^ 2) ./ (1 - rc .^ 2);
        
        f3 = prop.Qs ./ prop.Qa .* ones(size(f1));
        
        Lambda = max(zeros(size(f1)), min(min(min(f1, f2), f3), ones(size(f1))));


    %-- Particle streamline approach -------------------------------------%
    %   Diffusional Transmission Efficiency from Johnson (2018). 
    %   Orginally defined by Karlsson (2003).
    %   Code adapted from Nishan Sapkota. 
    case {'ps', 'particle-streamline'}
        addpath autils;
        
        % Transmission constants (Johnson, 2018)
        prop.lambda_e = 0.8; % classifier entrance/exit transmission efficiency 
        prop.L_eff = 46;  % length of a circular tube with the same diffusion deposition as the classifier 
        
        % Diffusion properties.
        % Invoke mass-mobility relationship for dm calculation.
        if opts.diffusion
            if ~isfield(prop, 'm0')
                warning('Mass-mobility information not given for AAC transfer function. Assuming dm = da.')
                dm = da;  % assume dm = da, otherwise must add massmob-info, dm2da(da, prop);
            else
                dm = dm2da(da, prop);  % invoke autils function for conversion
            end
            B = dm2zp(dm, [], prop.T, prop.p);  % compute mobility
            
            kB = 1.3806488e-23;  % Boltzmann's constant
            D = kB .* prop.T .* B;  % diffusion coefficient
        else  % neglect diffusion
            D = 0;
        end
        
        % Transfer width using empirical correlation.
        if ~opts.ideal
            % Incorporate deposition parameter.
            del_dep = (prop.L_eff .* D) ./ prop.Qa;
            fl = del_dep >= 0.007;
            
            lambda_d = zeros(size(del_dep));
            lambda_d(fl) = 0.891 .* exp(-11.5.*del_dep(fl)) + ...
                0.0975 .* exp(-70.1.*del_dep(fl)) + ...
                0.0325 .* exp(-179.*del_dep(fl));
            lambda_d(~fl) = 1 - 5.50 .* (del_dep(~fl).^(2/3)) + ...
                3.77 .* del_dep(~fl) + 0.814 .* (del_dep(~fl).^(4/3));
            
            lambda_tf = prop.lambda_e .* lambda_d;  % total transmission efficiency

            a = -1.73;  b = -0.0316;  c = 1.999;
            mu_tf = a .* (da.*1e9).^b + c;
        else
            lambda_tf = 1;
            mu_tf = 1;
        end
    
        % User defined transfer function variables for transmission efficiency.
        A0 = (lambda_tf .* (mu_tf.^2)) ./ (2.*prop.bet);
        B0 = prop.bet ./ mu_tf;
        
        tau_norm = tau ./ tau_star;

        Lambda = A0 .* (abs(tau_norm - (1 + B0)) + ...
            abs(tau_norm - (1 - B0)) - ...
            abs(tau_norm - (1 + B0 .* prop.del)) - ...
            abs(tau_norm - (1 - B0 .* prop.del)));
end


end



function d_star = cacl_d_star(tau_star, mu)
% A function to compute d_star when note directly supplied.

d_star = zeros(size(tau_star));

for ii=1:length(tau_star)  % loop through setpoints
    d_star(ii) = fminsearch(@(d) ...  % only used for output
        1e9 .* norm(Cc(d) .* 1e3 .* d.^2 ./ (18 .* mu) - ...
        tau_star(ii)), 100e-9);
end

end
