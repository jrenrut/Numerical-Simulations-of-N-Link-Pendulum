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
    offset = p.offset; wf = p.wf; phi = p.phi;
    
    % Unpack state
    thetas = z(1:n);
    omegas = z(n+1:2*n);
    
    dist_ax = -wf^2*offset.*cos(wf*t)*sin(phi);
    dist_ay = -wf^2*offset.*cos(wf*t)*cos(phi);
    d_a  = [dist_ax dist_ay 0];
    
    % Optional: Can turn off oscillation after a period, half period, etc.
%     if wf2*t > pi
%         d_a = zeros(1, 3);
%     end
    
    % Symbolic solutions
    M = feval(Mf, I, L, d, g, m, omegas, thetas, rest, d_a);
    b = feval(bf, I, L, d, g, m, omegas, thetas, rest, d_a);
    x = M\b;
    alphas = x(2*n+1:2*n+n);
    
    % Return state
    nlink_DAE = [omegas; alphas];
end