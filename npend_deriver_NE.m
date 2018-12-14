%-------------------------------------------------------------------------%
% npend_deriver_NE - Jeremy Turner
% 
% Symbolic solution for the equations of motion of the n-link pendulum
% using the Newton-Euler method (Linear and Angular momentum balances).
%
% Input: n - number of links
%
% Creates two files, npend_alphas_NE_M_[n].m and
% npend_alphas_Lagrange_NE_[n], which are the [M] and [b] matrices such
% that the solutions [x] = [M]^-1 [b].
%
% Extra: Includes simple friction force.
% ------------------------------------------------------------------------%

function npend_deriver_NE = func(n)

    % Symbolic variables
    syms g;
    m = sym('m', [n 1]);
    L = sym('L', [n 1]);
    d = sym('d', [n 1]);
    I = sym('I', [n 1]);
    theta = sym('theta', [n 1], 'real');
    omega = sym('omega', [n 1], 'real');
    alpha = sym('alpha', [n 1], 'real');
    syms fric

    % Cartesian unit vectors
    y = [0 1 0];
    z = [0 0 1];
    % Position vectors
    rOG = sym(zeros(n, 3));
    rEi = sym(zeros(n+1, 3));
    rEi(1,:) = zeros(1, 3);
    % Acceleration vectors
    aOG = sym(zeros(n, 3));
    aEi = zeros(1, 3);
    % Forces
    Fg = sym(zeros(n, 3));
    F  = sym(zeros(n, 3));
    Ff = sym(zeros(n, 3));
    % Torques
    M  = sym(zeros(n, 1));
    % Changes in angular momentum
    dH  = sym(zeros(n, 1));

    % Loop through all links
    for i=1:n
        
        eri = [sin(theta(i)) -cos(theta(i)) 0]; % Unit vector parallel to link i
        eti = [cos(theta(i)) sin(theta(i))  0]; % Unit vector in direction of link i's rotation
        
        % "One man's 'elbow below' is another man's 'elbow above'."
        % - Isaac Newton, maybe        
        rEG      = d(i)*eri; % Distance from elbow above above to CoM
        rBG      = (L(i) - d(i))*eri; % Distance from CoM to elbow below
        rOG(i,:) = rEi(i,:) + rEG; % Distance from origin to CoM
        rEi(i+1,:) = rEi(i,:) + L(i)*eri; % Distance from elbow below to origin

        aEG      = -rEG*omega(i)^2 + cross(alpha(i)*z, rEG); % Acceleration of elbow above w.r.t. CoM
        aOG(i,:) = aEi + aEG; % Acceleration of CoM w.r.t. origin
        aEi      = aEi - (rEG+rBG)*omega(i)^2 + cross(alpha(i)*z, (rEG+rBG)); % Acceleration of elbow below w.r.t. origin
        
        Fg(i,:)  = -m(i)*g*y; % Force of gravity
        F(i,:)   = m(i)*aOG(i,:); % m a
        Ff(i,:)  = -fric*0.01*(omega(i))*eti; % Force of friction
        
        % Loop through all el above i
        for j=1:i
            rEj   = rEi(j,:); % Distance from elbows above i to origin
            rEGj  = rOG(i,:) - rEj; % Distance from CoM of i to elbow j
            Mi    = cross(rEGj, Fg(i,:)) + cross(rEGj, Ff(i,:)); % Torque about j from forces on i
            M(j)  = M(j) + Mi(3); % Add to total torque about j
            dHi   = cross(rEGj, F(i,:)) + I(i)*alpha(i)*z; % Change in angular momentum about j from motion of i
            dH(j) = dH(j) + dHi(3); % Add to total change in angular momentum about j
        end
        
    end
    
    % Angular momentum balances for each link
    AMB = zeros(n, 1);
    for i=1:n
        AMB(i) = M(i) == dH(i);
    end
    
    % Differential equations for theta in matrix form
    [M, b] = equationsToMatrix(AMB, alpha');
    matlabFunction(M, 'file', ['Npend/npend_alphas_NE_M_', num2str(n)], 'Vars', {[I] [L] [d] g [m] [omega] [theta] fric});
    matlabFunction(b, 'file', ['Npend/npend_alphas_NE_b_', num2str(n)], 'Vars', {[I] [L] [d] g [m] [omega] [theta] fric});
    
end