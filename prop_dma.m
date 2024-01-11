
% PROP_DMA  Generates the prop struct used to summarize DMA parameters.
% Author:   Timothy Sipkens, 2019-07-17
%=========================================================================%

function [prop] = prop_dma(varargin)

%-- Parse inputs ---------------------------------------------------------%
if isempty(varargin); opts = struct(); opts.param = '';
elseif length(varargin) == 1
    opts = varargin{1};
    if ~isstruct(opts)
        t0 = opts;
        opts = struct();
        opts.param = t0;  % then assign supplied input to 'param' field
    end
elseif length(varargin) == 2
    opts.param = '';
    opts.model = varargin{1};
    opts.flow = varargin{2};
else; error('Too many inputs.')
end

% Parameter specs. Default is first option.
if ~isfield(opts, 'param'); opts.param = ''; end
switch opts.param
    case '' % Do nothing.
    case 'olfert'
        opts.model = '3018';
        opts.flow = 'low';
        
    case 'buckley'
        opts.model = '3081';
        opts.flow = 'low';  % custom flow rates in original
        
    otherwise
        if ~isempty(opts.param)
            opts.model = opts.param;  % by default, assume param is just the model number
        end
end

% Otherwise assign defaults.
if ~isfield(opts, 'model'); opts.model = ''; end
if isempty(opts.model); opts.model = '3018'; end  % default is 3018
prop.model = opts.model;

if ~isfield(opts, 'flow'); opts.flow = ''; end
if isempty(opts.flow); opts.flow = 'low'; end  % default is low flow
prop.flow = opts.flow;

if ~isfield(opts,'solver'); opts.solver = 'fullydeveloped'; end % if order not specified
prop.solver = opts.solver;

% Default temperature and pressure. 
prop.T = 293; % default temperature [K]
prop.p = 1;


%-- Assign classifier dimensions -----------------------------------------%
switch prop.model
    case {'3018', '3081', '3080'}  % DMA properties from Olfert/Buckley
        prop.L = 0.44369; % length, m
        prop.R1 = 0.00937; % inner radius, m
        prop.R2 = 0.01961; % outer radius, m

    case 'custom'
        prop = opts.prop; % read in from opts
end


%-- Assign flow rates ----------------------------------------------------%
switch prop.flow
    case 'low'
        prop.Q_s = 0.3/60/1000; % Sample flow [m^3/s]
        prop.Q_a = 0.3/60/1000; % Aerosol flow [m^3/s]
        prop.Q_c = 3/60/1000; % Sheath flow [m^3/s]
        prop.Q_m = 3/60/1000; % Exhaust flow [m^3/s]

    case 'high'
        prop.Q_s = 1.5/60/1000;
        prop.Q_a = 1.5/60/1000;
        prop.Q_c = 15/60/1000;
        prop.Q_m = 15/60/1000;

    case 'buckley'
        prop.Q_c = 4.89E-3/60; % sheath flow [m^3/s]
        prop.Q_a = 1.02E-3/60; % aerosol flow [m^3/s]
        prop.Q_s = Q_a; % sample flow [m^3/s]
        prop.Q_m = Q_c; % exhaust flow [m^3/s]
end


prop.bet = (prop.Q_s + prop.Q_a)  /(prop.Q_c + prop.Q_m);
prop.del = (prop.Q_s - prop.Q_a) / (prop.Q_s + prop.Q_a);

prop.Rd = 1 / prop.bet;  % compute theoretical resolution

gam = (prop.R1/prop.R2)^2;
kap = prop.L/prop.R2; % Stolzenburg Manuscript Eqn 9
% kap = L*R2/(R2^2-R1^2); % Buckey et al.


%-- Calculate G_DMA ------------------------------------------------------%
switch opts.solver
    case 'Buckley' % evaluation from Buckley et al.
        I_gamma = (0.25*(1-gam^2)*(1-gam)^2+(5/18)*(1-gam^3)*(1-gam)*log(gam) + ...
            (1/12)*(1-gam^4)*(log(gam))^2)/((1-gam)*(-0.5*(1+gam)*log(gam)-(1-gam))^2);
        prop.G_DMA = 4*(1+prop.bet)^2/(1-gam)*(I_gamma+(2*(1+prop.bet)*kap)^(-2));

    case 'plug' % Plug flow, Stolzenburg, 2018
        omega_a = 1-prop.bet*(1-prop.del)/(2*(1+prop.bet))*(1-gam);
        omega_s = gam+prop.bet*(1+prop.del)/(2*(1+prop.bet))*(1-gam);

        I_gamma = @(omega) (1-omega^2)/(2*(1-gam));
        G_o_corr = (4*(1+prop.bet)^2/(1-gam)*(I_gamma(omega_s)-I_gamma(omega_a))+...
            (omega_a-omega_s)/kap^2);
        prop.G_DMA = G_o_corr*(-log(gam)/2);

    case {'Olfert','fullydeveloped'} % Olfert laboratory; Fully developed flow, Stolzenburg, 2018
        fun = @(omega) ((1-omega)^2*log(gam)/2/(1-gam)+omega*log(omega)+(1-omega))/((1+gam)*log(gam)/2+(1-gam))-prop.bet*(1-prop.del)/2/(1+prop.bet);
        guess = 0.9;
        omega_a = fzero(fun,guess);

        fun = @(omega) (1-((1-omega)^2*log(gam)/2/(1-gam)+omega*log(omega)+(1-omega))/((1+gam)*log(gam)/2+(1-gam)))-prop.bet*(1+prop.del)/2/(1+prop.bet);
        guess = 0.3;
        omega_s = fzero(fun,guess);

        A = (-1/2*(1+gam)*log(gam)-(1-gam))^(-1);
        I_gamma = @(omega) A^2/(1-gam)*...
            (-omega^2*((1-gam)*log(omega)-(1-omega)*log(gam))^2/2+...
            (omega^2*(1-gam)/2+omega^3*log(gam)/3)*((1-gam)*log(omega)-...
            (1-omega)*log(gam))+(1-omega^2)*(1-gam)^2/4+...
            5*(1-omega^3)*(1-gam)*log(gam)/18+(1-omega^4)*(log(gam))^2/12);

        G_o_corr = (4*(1+prop.bet)^2)*(I_gamma(omega_s)-I_gamma(omega_a))/(1-gam)...
            +(omega_a-omega_s)/kap^2; % Stolzenburg Manuscript Equation 8 and 12
        prop.G_DMA = G_o_corr*log(prop.R2/prop.R1); % Stolzenburg Manuscript Equation 23

    otherwise
        disp('Invalid solver specified.');
        return;

end

end

