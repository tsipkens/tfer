
% PROP_AAC  Set up property structure for an AAC.
%  AUTHOR: Timothy Sipkens, 2024-01-10

function prop = prop_aac()

% Default temperature and pressure. 
prop.T = 293; % default temperature [K]
prop.p = 1;

% Classifier dimensions.
prop.r1 = 0.056; % 0.043;
prop.r2 = 0.060; % 0.045;
prop.L = 0.206; % 0.21;

%-{
% Low flow.
prop.Qs = 0.3/60/1000;  % sample flow (leaving classifier) [m^3/s]
prop.Qa = prop.Qs;  % aerosol flow (entering classifier) [m^3/s]
prop.Qsh = 3/60/1000;  % sheath flow [m^3/s]
prop.Qexh = prop.Qsh + prop.Qa - prop.Qs;  % outlet sheath [m^3/s]
%}

%{
% High flow.
prop.Qs = 1.5/60/1000;  % sample flow (leaving classifier) [m^3/s]
prop.Qa = prop.Qs;  % aerosol flow (entering classifier) [m^3/s]
prop.Qsh = 15/60/1000;  % sheath flow [m^3/s]
prop.Qexh = prop.Qsh + prop.Qa - prop.Qs;  % outlet sheath [m^3/s]
%}

end