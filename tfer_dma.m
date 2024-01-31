
% TFER_DMA  Evaluates the transfer function of a differential mobility analyzer.
% 
%  OMEGA = kernel.tfer_dma(D_STAR, D, Z) uses the mobility diameter set points
%  specified by D_STAR [nm] and evalautes the DMA transfer function at D
%  [nm] for an integer charge state of Z (a scalar, integer). Uses default
%  properties specified by kernel.prop_dma. Explicitly stating prop_dma is
%  preferred. Output is transfer function, OMEGA. 
%  
%  NOTE: D should be n2 x 1. D_STAR should be 1 x n1. If any entries of
%  prop are vectors,they should have the same dimensions as D_STAR.
%  
%  OMEGA = kernel.tfer_dma(D_STAR, D, Z, PROP_DMA) explicitly specified the
%  properties of the DMA (e.g., the radii) as a data structure. This is the
%  preferred usage to the previous call. For structure of PROP_DMA, see
%  kernel.prop_dma(...). 
%  
%  OMEGA = kernel.tfer_dma(...,OPTS) adds an options structure with the
%  field specified below: 
%   OPTS.diffusion  Indicates whether to include diffusion, boolean (optional)
%       .solver     Indicates the method by which diffusion is calculated (optional)
%       .param      String indicated which parameter set to use (see prop_DMA.m)
%  
%  [OMEGA, ZP_TILDE] = kernel.tfer_dma(...) adds an output containing the
%  non-dimensional electrical mobility as a vector.
%  
%  [OMEGA, ZP_TILDE, PROP] = kernel.tfer_dma(...) adds output for updated
%  property struct, containing transfer function specific information
%  (e.g., DMA resolution). 
%  
%  [..., SP] = kernel.tfer_dma(...) adds an output containing setpoint
%  information. 
%  
%  ------------------------------------------------------------------------
% 
%  AUTHOR: Timothy Sipkens, 2018-12-27
%  Some code adapted from Buckley et al. (2017) and Olfert group.

function [Omega, Zp_tilde, prop, sp] = tfer_dma(d_star, d, z, prop, opts)


%-- Parse inputs ---------------------------------------------------------%
% Add mat-tfer-pma package to MATLAB path.
% Used for size conversions calculations.
fd = fileparts(mfilename('fullpath'));
addpath([fd, filesep, 'autils']);

% Set defaults, if opts not given.
if ~exist('opts', 'var'); opts = struct(); end  % initialize options struct
if ~isfield(opts, 'solver'); opts.solver = 'fullydeveloped'; end
if ~isfield(opts, 'diffusion'); opts.diffusion = 1; end
if ~isfield(opts, 'type'); opts.type = 'd_star'; end

% Get prop, if prop not given. Use opts struct to build.
if ~exist('prop', 'var'); prop = []; end
if isempty(prop); prop = prop_dma(opts); end

d = d .* 1e-9;  % convert from nm to m for calculations
%-------------------------------------------------------------------------%


%-- Compute additional classifier properties -----------------------------%
prop.bet = (prop.Qs + prop.Qa) ./ (prop.Qc + prop.Qm);
prop.del = (prop.Qs - prop.Qa) ./ (prop.Qs + prop.Qa);

% Check if bet and del change. If not, reduce to scalar for calcs. (faster).
if and(all(prop.bet == prop.bet(1)), all(prop.del == prop.del(1)))
    prop.bet = prop.bet(1);
    prop.del = prop.del(1);
end

prop.Rd = 1 ./ prop.bet;  % compute theoretical resolution

gam = (prop.R1 ./ prop.R2) .^ 2;
kap = prop.L ./ prop.R2;  % Stolzenburg Manuscript, Eq. 9
% kap = prop.L .* prop.R2 ./ (prop.R2.^2 - prop.R1.^2);  % Buckley et al.


%-- Parse setpoint input -------------------------------------------------%
switch opts.type
    case 'd_star'  % default
        d_star = d_star .* 1e-9;  % convert from nm to m for calculations

        %-- Evaluate particle mobility -----------------------------------%
        if strcmp(opts.solver, 'buckley')
            [~, Zp_star] = dm2zp(d_star, 1);  % evaluate electrical mobility (Davies)
        else
            [~, Zp_star] = dm2zp(d_star, 1, prop.T, prop.p);  % evaluate electrical mobility (Kim et al.)
        end
        
        % Classifier voltage (TSI DMA 3080 Manual Equation B-5).
        V = (prop.Qc ./ ((2*pi) .* Zp_star .* prop.L)) .* ...
            log(prop.R2 ./ prop.R1);

    case 'V'
        V = d_star;  % input is actually voltage, copy over

        % Classifier voltage to mobility.
        Zp_star = (prop.Qc ./ ((2*pi) .* V .* prop.L)) .* ...
            log(prop.R2 ./ prop.R1);
        
        %-- Evaluate particle mobility -----------------------------------%
        if strcmp(opts.solver, 'buckley')
            d_star = zp2dm(Zp_star, 1);  % evaluate electrical mobility (Davies)
        else
            d_star = zp2dm(Zp_star, 1, prop.T, prop.p);  % evaluate electrical mobility (Kim et al.)
        end
end

% Generate setpoint structure (for output).
sp(length(d_star)) = struct();
for ii=1:length(d_star)
    sp(ii).d_star = d_star(ii);
    sp(ii).Zp_star = Zp_star(ii);
    sp(ii).V = V(ii);
end


%-- Calculate G_DMA ------------------------------------------------------%
%   Pre-computed before looping through charge states below.
switch opts.solver
    case 'buckley' % evaluation from Buckley et al.
        I_gamma = (0.25*(1-gam^2)*(1-gam)^2+(5/18)*(1-gam^3)*(1-gam)*log(gam) + ...
            (1/12)*(1-gam^4)*(log(gam))^2)/((1-gam)*(-0.5*(1+gam)*log(gam)-(1-gam))^2);
        prop.G_DMA = 4 .* (1 + prop.bet).^2 ./ (1-gam) .* ...
            (I_gamma + (2 .* (1 + prop.bet) .* kap) .^ (-2));
        
    case 'plug' % Plug flow, Stolzenburg, 2018
        omega_a = 1 - prop.bet .* (1 - prop.del) ./ (2 .* (1 + prop.bet)) .* (1 - gam);
        omega_s = gam + prop.bet .* (1 + prop.del) ./ (2 .* (1 + prop.bet)) .* (1 - gam);
        
        I_gamma = @(omega) (1 - omega.^2) ./ (2 .* (1 - gam));
        G_o_corr = (4 .* (1 + prop.bet).^2 ./ (1 - gam) .* ...
            (I_gamma(omega_s) - I_gamma(omega_a)) + ...
            (omega_a - omega_s) ./ kap.^2);
        prop.G_DMA = G_o_corr .* (-log(gam) ./ 2);

    case {'olfert', 'fullydeveloped'} % Olfert laboratory; Fully developed flow, Stolzenburg, 2018
        if or(any(size(prop.bet) > 1), any(size(prop.del) > 1))
            error('Cannot currently compute transfer function for when del/bet change.');
        end
        
        omega_a = zeros(size(prop.bet));  omega_s = omega_a;  % initialize omegas
        for ii=1:length(prop.bet)  % loop over bet/del (as they may be vectors)
            fun = @(omega) ((1-omega).^2 .* log(gam)./2 ./ (1 - gam) + ...
                omega .* log(omega) + (1-omega)) ./ ((1+gam) .* log(gam)./2 + ...
                (1-gam)) - prop.bet(ii) .* (1-prop.del(ii))./2 ./ (1+prop.bet(ii));
            omega_a(ii) = fzero(fun, 0.9);
            
            fun = @(omega) (1 - ((1-omega).^2 .* log(gam)./2 ./ (1-gam) + ...
                omega .* log(omega) + (1-omega)) ./ ((1+gam) .* log(gam)./2 + ...
                (1-gam))) - prop.bet(ii) .* (1+prop.del(ii))./2 ./ (1+prop.bet(ii));
            omega_s(ii) = fzero(fun, 0.3);
        end

        A = (-1/2*(1+gam) .* log(gam) - (1-gam)).^(-1);
        I_gamma = @(omega) A.^2 ./ (1-gam) .* ...
            (-omega^2 .* ((1-gam) .* log(omega) - (1-omega) .* log(gam)).^2./2 + ...
            (omega.^2 .* (1-gam)./2 + omega.^3 .* log(gam)./3) .* ((1-gam) .* log(omega) - ...
            (1-omega) .* log(gam)) + (1-omega.^2) .* (1-gam).^2./4 + ...
            5 .* (1-omega.^3) .* (1-gam) .* log(gam)./18 + (1-omega.^4) .* (log(gam)).^2./12);
        
        G_o_corr = (4 .* (1+prop.bet).^2) .* ...
            (I_gamma(omega_s) - I_gamma(omega_a)) ./ (1-gam) +...
            (omega_a - omega_s) ./ kap.^2;  % Stolzenburg Manuscript Equation 8 and 12
        prop.G_DMA = G_o_corr .* log(prop.R2/prop.R1);  % Stolzenburg Manuscript Equation 23

    otherwise
        disp('Invalid solver specified.');
        return;
end


%-- Loop through charge states and evaluate transfer function ------------%
Omega = zeros(length(d_star), length(d), length(z));
Zp_tilde = zeros(length(d), length(d_star), length(z));
for ii=1:length(z)
    if z(ii) == 0  % then hard code that no particles transmit
        Omega(:,:,ii) = zeros(length(d_star), length(d));
    else
        [Omega(:,:,ii), Zp_tilde(:,:,ii)] = ...
            tfer_dma0(Zp_star, V, d, z(ii), prop, opts);
    end
end

end



%== TFER_DMA0 ============================================================%
function [Omega, Zp_tilde] = tfer_dma0(Zp_star, V, d, z, prop, opts)
% TFER_DMA0  Evalutes transfer function for a single charge state.
% Loop through charge states and pre-compute values in previous function.

%-- Physical constants ---------------------------------------------------%
kb = 1.38064852e-23; % Boltzmann constant [m^2 kg s^-2 K^-1]
e = 1.6022E-19; % electron charge [C]

%-- Evaluate particle mobility -------------------------------------------%
if strcmp(opts.solver, 'buckley')
    [B, Zp] = dm2zp(d, z);  % evaluate electrical mobility (Davies)
else
    [B, Zp] = dm2zp(d, z, prop.T, prop.p);  % evaluate electrical mobility (Kim et al.)
end
Zp_tilde = Zp ./ Zp_star;  % array of non-dimensional mobilities


%-- Evaluate transfer function -------------------------------------------%
if opts.diffusion
    switch opts.solver
        case 'buckley'  % evaluation from Buckley et al.
            D = prop.D(B) .* z;  % diffusion
            sigma = sqrt(prop.G_DMA .* (2*pi) .* prop.L .* D ./ prop.Qc);
            
        case 'plug'  % plug flow, Stolzenburg, 2018
            sigma_star = sqrt(((kb*prop.T) / (z*e*V)) * prop.G_DMA); % Stolzenburg Manuscript Equation 20
            sigma = sqrt(sigma_star.^2 .* Zp_tilde); % Stolzenburg Manuscript Equation 19
            
        case {'olfert', 'fullydeveloped'}  % Olfert laboratory; Fully developed flow, Stolzenburg, 2018
            sigma_star = sqrt(((kb*prop.T) ./ (z .* e .* V)) .* prop.G_DMA); % Stolzenburg Manuscript Equation 20
            sigma = sqrt(sigma_star.^2 .* Zp_tilde); % Stolzenburg Manuscript Equation 19
            
        otherwise
            disp('Invalid solver specified.');
            return;
            
    end
    
    epsilon = @(x) x .* erf(x) + exp(-x.^2) ./ (pi^0.5);
    
    %-- Standard DMA transfer function for diffusion --------%
    %	Computes the diffusive transfer function for the DMA
    %	based on Stolzenberg's 1988 analysis. 
    Omega = sigma ./ (2^0.5*prop.bet*(1 - prop.del)).*(...
        epsilon((Zp_tilde - 1 - prop.bet)./(2^0.5*sigma))+...
        epsilon((Zp_tilde - 1 + prop.bet)./(2^0.5*sigma))-...
        epsilon((Zp_tilde - 1 - prop.bet*prop.del)./(2^0.5*sigma))-...
        epsilon((Zp_tilde - 1 + prop.bet*prop.del)./(2^0.5*sigma)));
    
else  % simpler evaluation for the case of exluding diffusion
    Omega = 1 / (2 .* prop.bet .* (1 - prop.del)).*(abs(Zp_tilde - 1 + prop.bet)+...
        abs(Zp_tilde - 1 - prop.bet) - ...
        abs(Zp_tilde - 1 + prop.bet .* prop.del) - ...
        abs(Zp_tilde - 1 - prop.bet .* prop.del));
end
%-------------------------------------------------------------------------%


% Remove negative and small values. 
negTransform = Omega > 0;
Omega = negTransform.*Omega;
smallTransform = Omega > 1e-14;
Omega = smallTransform .* Omega;

% Remove numerical noise in kernel.
Omega(Omega < (1e-7 .* max(max(Omega)))) = 0;

Omega = Omega'; % transpose data to output desired format

end

