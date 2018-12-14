%-------------------------------------------------------------------------%
% nlink_DAE - Jeremy Turner
% 
% Right-hand side function called by ode45 in npend.m to simulate system
% using the symbolically-derived ODEs using the DAE method. Includes
% shaking end.
%
% Input: z - State vector 2nx1 (theta1; ...; thetan; omega1; ...; omegan)
%        p - Parameter struct
%        t - timespan as defined in npend.m
%       Mf - [M] matrix symbolically derived function file name
%       bf - [b] vector symbolically derived function file name
%
% Returns: State vector at new timestep
% ------------------------------------------------------------------------%

function nlink_DAE = func(z, p, t, Mf, bf)

    % Unpack parameters
    g = p.g; L = p.L; d = p.d; m = p.m; I = p.I; n = p.n; rest = p.rest;
    offset2 = p.offset2; wf2 = p.wf2; phi2 = p.phi2;
    
    % Unpack state
    thetas = z(1:n);
    omegas = z(n+1:2*n);
    
    % Optional: Can turn off oscillation after a period, half period, etc.
%     if wf2*t > pi
%         d2_a = zeros(1, 3);
%     else
        dist_ax2 = -wf2^2*offset2.*cos(wf2*t)*sin(phi2);
        dist_ay2 = -wf2^2*offset2.*cos(wf2*t)*cos(phi2);
        d2_a  = [dist_ax2 dist_ay2 0];
%     end
    
    % Symbolic solutions
    M = feval(Mf, I, L, d, g, m, omegas, thetas, rest, d2_a);
    b = feval(bf, I, L, d, g, m, omegas, thetas, rest, d2_a);
    x = M\b;
    alphas = x(2*n+1:2*n+n);
    
    % Return state
    nlink_DAE = [omegas; alphas];
end