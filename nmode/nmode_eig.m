%-------------------------------------------------------------------------%
% nmode_eig - Jeremy Turner
% 
% Function called in npend.m to find the normal modes of oscillation by
% solving M\K using system parameters.
%
% Input: p - Parameter struct
%       Mf - [M] matrix symbolically derived function file name
%       Kf - [K] matrix symbolically derived function file name
%
% Returns: n x 2n matrix, where first n vectors are the initial theta
% values for the n normal modes, and the last n columns form a diagonal
% matrix whose elements are the squared frequencies of normal mode
% oscillations.
% ------------------------------------------------------------------------%

function nmode_eig = func(p, Mf, Kf)

    % Unpack parameters
    g = p.g; L = p.L; d = p.d; m = p.m; I = p.I; n = p.n;
    
    % Dummy vectors
    omega = zeros(n, 1);
    theta = zeros(n, 1);

    Mf = ['nmodeM_', num2str(n)];
    Kf = ['nmodeK_', num2str(n)];
    M = feval(Mf, I, L, d, g, m, omega, theta);
    K = feval(Kf, I, L, d, g, m, omega, theta);

    % Symbolic solutions
    [v, d] = eig(M\K);
    
    % Return normal modes
    nmode_eig = [v, d];
end