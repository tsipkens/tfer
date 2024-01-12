
% TFER_1C_PB    Evaluates the transfer function for a PMA in Case B (w/ parabolic flow).
% Author:       Timothy Sipkens, 2019-03-21
% 
% Inputs:
%   sp          Structure defining various setpoint parameters 
%               (e.g. m_star, V). Use 'get_setpoint' method to generate 
%               this structure.
%   m           Particle mass
%   d           Particle mobility diameter
%   z           Integer charge state
%   prop        Device properties (e.g. classifier length)
%
% Outputs:
%   Lambda      Transfer function
%   G0          Function mapping final to initial radial position
%=========================================================================%

function [Lambda, G0] = tfer_1C_pb(sp, m, d, z, prop)

[tau,C0,~,rs] = parse_inputs(sp,m,d,z,prop);
        % parse inputs for common parameters

%-- Taylor series expansion constants ------------------------------------%
C3 = tau .* ([sp.alpha]' .^ 2 .* prop.rc + ...
    2 .* [sp.alpha]' .* [sp.beta]' ./ prop.rc + ...
    [sp.beta]' .^ 2 ./ (prop.rc .^ 3) - C0 ./ (m .* prop.rc));
C4 = tau .* ([sp.alpha]' .^ 2 - ...
    2 .* [sp.alpha]' .* [sp.beta]' ./ (prop.rc .^ 2) - ...
    3 .* [sp.beta]' .^ 2 ./ (prop.rc .^ 4) + ...
    C0 ./ (m .* (prop.rc .^ 2)));

A1 = -3 * prop.v_bar ./ (4 .* C4 .^ 3 .* prop.del ^ 2);
A2 = 2 .* (C3 .^ 2 - C4 .^ 2 .* prop.del ^ 2);

% Loop over setpoints.
Lambda = zeros(size(A1));
for jj=1:size(Lambda,1)
    A3 = @(r,ii) C4(jj,ii) .^ 2 .* (r-prop.rc) .^ 2 - 2 .* C3(jj,ii) .* C4(jj,ii) .* (r-prop.rc);
    
    %-- Set up F function for minimization -------------------------------%
    F = @(r,ii) A1(jj,ii) .* (A2(jj,ii) .* ...
        log(abs(C4(jj,ii) .* (r - prop.rc) + C3(jj,ii))) + A3(r,ii));
    min_fun = @(rL,r0,ii) F(rL,ii) - F(r0,ii) - prop.L;
    
    %-- Evaluate G0 and transfer function --------------------------------%
    G0 = @(r) G_fun(min_fun, r, rs(jj,:), ...
        prop.r1, prop.r2, sp(jj).alpha, sp(jj).beta);
    
    ra = min(prop.r2, max(prop.r1, G0(prop.r1)));
    rb = min(prop.r2, max(prop.r1, G0(prop.r2)));
    
    Lambda(jj,:) = 3/4 .* (rb - ra) ./ prop.del - ...
        1/4 .* ((rb - prop.rc) ./ prop.del) .^ 3 + ...
        1/4 .* ((ra - prop.rc) ./ prop.del) .^ 3;
end

end

