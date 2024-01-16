
% PROP_DMA  Generates the prop struct used to summarize DMA parameters.
%  AUTHOR: Timothy Sipkens, 2019-07-17

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
%-------------------------------------------------------------------------%


% Parameter specs. Default is first option.
if ~isfield(opts, 'param'); opts.param = ''; end
switch opts.param
    case ''  % do nothing
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
        prop.Qs = 0.3/60/1000; % Sample flow [m^3/s]
        prop.Qa = 0.3/60/1000; % Aerosol flow [m^3/s]
        prop.Qc = 3/60/1000; % Sheath flow [m^3/s]
        prop.Qm = 3/60/1000; % Exhaust flow [m^3/s]

    case 'high'
        prop.Qs = 1.5/60/1000;
        prop.Qa = 1.5/60/1000;
        prop.Qc = 15/60/1000;
        prop.Qm = 15/60/1000;

    case 'buckley'
        prop.Qc = 4.89E-3/60; % sheath flow [m^3/s]
        prop.Qa = 1.02E-3/60; % aerosol flow [m^3/s]
        prop.Qs = Q_a; % sample flow [m^3/s]
        prop.Qm = Q_c; % exhaust flow [m^3/s]
end



end

