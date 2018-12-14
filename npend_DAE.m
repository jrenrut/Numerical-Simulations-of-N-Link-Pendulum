%-------------------------------------------------------------------------%
% npend_DAE - Jeremy Turner
% 
% Right-hand side function called by ode45 in npend.m to simulate system
% using the symbolically-derived ODEs using the DAE method. Includes
% shaking base.
%
% Input: z - State vector 2nx1 (theta1; ...; thetan; omega1; ...; omegan)
%        p - Parameter struct
%        t - timespan as defined in npend.m
%       Mf - [M] matrix symbolically derived function file name
%       bf - [b] vector symbolically derived function file name
%
% Returns: State vector at new timestep, including the new position of the
% base.
% ------------------------------------------------------------------------%

function npend_DAE = func(z, p, t, Mf, bf)

    % Unpack parameters
    g = p.g; L = p.L; d = p.d; m = p.m; I = p.I; n = p.n;
    offset1 = p.offset1; wf1 = p.wf1; phi1 = p.phi1;
    
    % Unpack state
    thetas = z(1:n);
    omegas = z(n+1:2*n);
    
    % Calculate new position of base
    dist_x = offset1.*sin(wf1*t)*sin(phi1);
    dist_y = offset1.*sin(wf1*t)*cos(phi1);
    d_r = [dist_x dist_y 0];
    
    % Calculate new acceleration of base
    dist_ax = -wf1^2*offset1.*sin(wf1*t)*sin(phi1);
    dist_ay = -wf1^2*offset1.*sin(wf1*t)*cos(phi1);
    d_a = [dist_ax dist_ay 0];
    
    if wf1*t > pi
        d_r = zeros(1, 3);
        d_a = zeros(1, 3);
    end
    
    % Symbolic solutions
    M = feval(Mf, I, L, d, g, m, omegas, thetas, d_a);
    b = feval(bf, I, L, d, g, m, omegas, thetas, d_a);
    x = M\b;
    alphas = x(2*n+1:2*n+n);
    
    % Return state
    npend_DAE = [omegas; alphas; d_r.'];
end