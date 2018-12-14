%-------------------------------------------------------------------------%
% npend_Lagrange - Jeremy Turner
% 
% Right-hand side function called by ode45 in npend.m to simulate system
% using the symbolically-derived ODEs using the Lagrange method.
%
% Input: z - State vector 2nx1 (theta1; ...; thetan; omega1; ...; omegan)
%        p - Parameter struct
%       Mf - [M] matrix symbolically derived function file name
%       bf - [b] vector symbolically derived function file name
%
% Returns: State vector at new timestep
% ------------------------------------------------------------------------%

function npend_Lagrange = func(z, p, Mf, bf)
   
    % Unpack parameters
    g = p.g; L = p.L; d = p.d; m = p.m; I = p.I; n = p.n;
    
    % Unpack state
    thetas = z(1:n);
    omegas = z(n+1:2*n);
    
    % Symbolic solutions
    M = feval(Mf, I, L, d, g, m, omegas, thetas);
    b = feval(bf, I, L, d, g, m, omegas, thetas);
    alphas = M\b;
    
    % Return new state
    npend_Lagrange = [omegas; alphas];
end