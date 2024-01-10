
% TFER_TRI Evaluates a triangular transfer function from a setpoint and resolution.
%  
%  LAMBDA = tfer_tri(S_STAR, S, R) computes the transfer function using the
%  set of setpoint in S_STAR and resolutions in R on the array of sizes in
%  S. 
%  
%  LAMBDA = tfer_tri(SP, S, ZET) extracts setpoint information and 
%  resolution from a PMA setpoint structure and uses the mass-mobility
%  exponents ZET to compute transfer function. In this case, one can also 
%  add multiple charging by adding Z argument.
%  
%  ------------------------------------------------------------------------
% 
%  AUTHOR: Timothy Sipkens, 2018-12-27, 2024-01-10

function Lambda = tfer_tri(s_star0, s, R0, z)

% Simple evaluation at a single charge state.
if ~isstruct(s_star0)
    % Check if resolution is a scalar and extend. 
    if all(size(R0) == 1)
        R0 = R0 .* ones(size(s_star0));
    end

    Lambda = tfer_tri0(s_star0, s, R0);  % simple evaluation
    

% Handle if PMA setpoint is given(allows for multiple charging). 
else
    zet = R0;
    R0 = [s_star0.Rm];
    s_star0 = [s_star0.m_star];

    % If charge states not given. 
    if ~exist('z', 'var'); z = []; end
    if isempty(z); z = 1; end  % only evaluate a single charge peak
    
    % Continue with calculation.
    Lambda = zeros(length(s), length(s_star0), length(z));
    for ii=1:length(z)  % loop over charge states (only used for z)
    
        s_star = s_star0 .* z(ii);  % incorporate multiple charging (= s_star for single charge)
        R = R0 .* z(ii) .^ (1./zet);

        %-- Actucal transfer function evaluation -----------------------------%
        Lambda(:,:,ii) = tfer_tri0(s_star, s, R);
    end
end

end



function Lambda = tfer_tri0(s_star, s, R)
% Simple transfer functione evaluation.

s_max = s_star .* (1 ./ R + 1);  % calculate maximum size for triangle (for this charge state)

s_del = s_max - s_star; % FWHM of the transfer function (related to resolution)
s_min = 2 .* s_star - s_max; % lower end of the transfer function

Lambda = ...
        (s <= s_star) .* (s > s_min) .* (s - s_min) ./ s_del + ...
        (s > s_star) .* (s < s_max) .* ((s_star - s) ./ s_del + 1);

end


