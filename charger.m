
% CHARGER  Calculates the fraction of particles with a specific integer charge.
%  Based on code from Buckley et al., J. Aerosol Sci. (2008) and Olfert
%  laboratory at the University of Alberta.
% 
%  INPUTS:
%   d      Particle diameter [nm]
%   z      Integer particle charge state (optional, default = 0:6)
%   T      Temperature [K] (optional, default = 298 K)
%   model  String that specifies which model to use (optional, default = 'hybrid')
% 
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2018-12-27

function [fn, qbar, model] = charger(d, z, T, model, opt)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('T','var')
    T = 298;
elseif isempty(T)
    T = 298;
end

if ~exist('z','var') % if z is not specified, output for states 0 to 6
    z = [];
elseif isempty(z)
    z = (0:6);
end
if size(z,2) > size(z,1)  % transpose depending on the input
    z = z';
end

if ~exist('model','var'); model = []; end
if strcmpi(model, 'fuchs'); model = 'f'; end
if strcmpi(model, 'li'); model = 'l'; end
if isempty(model)
    model = repmat('w', [1,length(z)]);
    model(abs(z) < 3) = 'g';
end

if ~exist('opt', 'var'); opt = []; end
if isempty(opt); opt = struct(); end

d = d .* 1e-9;  % convert from nm to m for calculations
%-------------------------------------------------------------------------%


e = 1.602177e-19; % elementary charge
epi = 8.85418e-12; % dielectric constant (for air) [F/m]
kB = 1.38065e-23; % Boltzmann's constant
Z_Z = 0.875; % ion mobility ratio (Wiedensohler, 1988)

[vec_d,vec_z] = ndgrid(d,z); % used in boolean expressions below
fn = zeros(size(vec_d));


% Wiedensohler.
if any(model == 'w')
    indw = (model == 'w')';
    if any(abs(z) < 3)  % if charge state less than 3
        ind = and(indw, abs(z) < 3);

        a = [-26.3328,-2.3197,-0.0003,-2.3484,-44.4756;
            35.9044,0.6175,-0.1014,0.6044,79.3772;
            -21.4608,0.6201,0.3073,0.4800,-62.8900;
            7.0867,-0.1105,-0.3372,0.0013,26.4492;
            -1.3088,-0.1260,0.1023,-0.1553,-5.7480;
            0.1051,0.0297,-0.0105,0.0320,0.5049];
                % coefficients for fit

        exponent = zeros(length(d),sum(ind));
        for ii = 1:6
            exponent = exponent + a(ii,z(ind)+3).*log10(d.*1e9).^(ii-1);
        end
        fn(:,ind) = 10.^exponent;

        fn(and(vec_d<20e-9, abs(vec_z)==2)) = 0;
    end

    if any(abs(z) >= 3) % if charge state is 3 or more
        ind = and(indw, abs(z) >= 3);

        fn(:,ind) = e./sqrt(4*pi*pi*epi*kB*T.*vec_d(:,ind)) .* ...
            exp(-(vec_z(:,ind) - 2*pi*epi*kB*T*log(Z_Z).*vec_d(:,ind)./e^2).^2 ./ ...
            (4*pi*epi*kB*T.*vec_d(:,ind)./e^2));

        fn(and(vec_d<69.78e-9, abs(vec_z)>=3)) = 0;
        fn(fn<=6e-5) = 0; % remove unnecessary small values
    end
end


% Gopalakrishnan.
if any(model == 'g')
    indg = (model == 'g')';
    if any(z<3)
        ind = and(indg, abs(z) < 3);
        
        if ~isfield(opt, 'conduction'); cond = 1;
        else; cond = opt.conduction; end
        a = get_a(cond);

        exponent = zeros(length(d),sum(ind));
        for ii = 1:4
            exponent = exponent + ...
                a(ii,z(ind)+3) .* log(d.*1e9) .^ (ii-1);
        end
        fn(:,ind) = exp(exponent);
    end
end
qbar = sum(fn .* z', 2) ./ sum(fn, 2);


if any(model == 'f')
    % Dielectric constant (default from T. Johnson).
    if ~isfield(opt, 'eps'); eps = 13.5;
    else; eps = opt.eps; end
    
    % ion·s/m3
    if ~isfield(opt, 'nit'); nit = 1.01e13;
    else; nit = opt.nit; end
    
    % If many diameters, create textbar.
    if length(d)>5
        disp(' Running Fuchs charging model:');
        tools.textbar([0,length(d)]);
    end
    
    % Main loop over diameters.
    qbar = zeros(length(d), 1);  % initialize
    for ii=1:length(d)
        [p0, pz, qbar(ii)] = fuchs(d(ii), max(80,round(max(z) .* 1.2)), T, 1, nit, eps);
        if z(1)~=0  % account for if z = 0 is present
            fn(ii,:) = pz(z);
        else
            fn(ii,2:end) = pz(z(2:end));
            fn(ii,1) = p0;
        end
        
        if length(d)>5; tools.textbar([ii,length(d)]); end
    end
    if length(d)>5; disp(' '); end
end


if any(model == 'l')
    % Dielectric constant (default from T. Johnson).
    if ~isfield(opt, 'eps'); eps = 13.5;
    else; eps = opt.eps; end
    
    % ion·s/m3
    if ~isfield(opt, 'nit'); nit = 1.01e13;
    else; nit = opt.nit; end

    % If many diameters, create textbar.
    if length(d)>5
        disp(' Running Li et al. charging model (interpolated):');
        tools.textbar([0,length(d)]);
    end

    % Unpack stored collision kernels.
    colls = [];
    fd = fileparts(mfilename('fullpath'));
    in = load([fd, filesep, 'li_v4_collkernel.mat']);
    for ii=1:length(in.dvec)
        if size(colls, 2) ~= length(in.collkernel0{ii})
            colls = [colls, zeros(size(colls, 1), ...
                length(in.collkernel0{ii}) - size(colls, 2))];
        end
        colls = cat(1, colls, in.collkernel0{ii});
    end

    % Main loop over diameters.
    qbar = zeros(length(d), 1);  % initialize
    for ii=1:length(d)
        idx = find(d(ii) * 1e9 < in.dvec, 1);

        if isempty(idx)  % input data is not available
            fn(ii, :) = 0;  % or NaN
            continue;
        end

        coll = exp((log(colls(idx,:)) - log(colls(idx - 1,:))) .* ...
            (log(d(ii) * 1e9) - log(in.dvec(idx - 1))) ./ ...
            (log(in.dvec(idx)) - log(in.dvec(idx - 1))) + ...
            log(colls(idx - 1,:)));
        
        [qbar(ii), pz] = collkernel2charge(coll, nit);
        
        z_int = intersect(0:size(colls, 2) - 1, z);
        z_xor = setxor(0:size(colls, 2) - 1, z);
        
        fn(ii, z_int + 1) = pz;
        fn(ii, z_xor + 1) = 0;  % or NaN
        
        if length(d)>5; tools.textbar([ii,length(d)]); end
    end
    if length(d)>5; disp(' '); end
end



fn = fn';

end


function [a] = get_a(cond) % coefficients from Gopalakrishnan et al.

if cond == 1 % conducting values
    a = [-45.405, -7.8696, -0.3880, -8.0157, -40.714;...
          20.049,  3.1036,  0.4545,  3.2536,  17.487;...
         -3.0570, -0.4557, -0.1634, -0.5018,  -2.6146;...
          0.1534,  0.0187,  0.0091,  0.0223,   0.1282];

else % non-conducting values
    a = [-63.185, -16.801,  -1.212,  -16.704,  -71.051;...
          26.833,   7.5947,  1.1068,   7.5438,  31.209;...
          -3.8723, -1.1975, -0.2934,  -1.1938,  -4.6696;...
           0.1835,  0.0590,  0.0169,   0.0589,   0.2301];
end

end


