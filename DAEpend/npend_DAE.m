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
    offset = p.offset; wf = p.wf; phi = p.phi;
    
    % Unpack state
    thetas = z(1:n);
    omegas = z(n+1:2*n);
    
    % Optional: Normal modes
%     nm = 3;
%     bj0 = @(x) besselj(0,x);
%     J0 = fzero(bj0,[(nm-1) nm]*pi);
%     wf = J0/2*sqrt(g/sum(L));
    
    % Calculate new position of base
    dist_x = offset.*cos(wf*t)*sin(phi);
    dist_y = offset.*cos(wf*t)*cos(phi);
    d_r = [dist_x dist_y 0];
    
    % Calculate new acceleration of base
    dist_ax = -wf^2*offset.*cos(wf*t)*sin(phi);
    dist_ay = -wf^2*offset.*cos(wf*t)*cos(phi);
    d_a = [dist_ax dist_ay 0];
    
    % Optional: Can turn off oscillation after a period, half period, etc.
%     if wf1*t > 2*pi
%         d_r = zeros(1, 3);
%         d_a = zeros(1, 3);
%     end
    
    % Symbolic solutions
    M = feval(Mf, I, L, d, g, m, omegas, thetas, d_a);
    b = feval(bf, I, L, d, g, m, omegas, thetas, d_a);
    x = M\b;
    alphas = x(2*n+1:2*n+n);
    
    % Return state
    npend_DAE = [omegas; alphas; d_r.'];
end