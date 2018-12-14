%-------------------------------------------------------------------------%
% nmode - Jeremy Turner
% 
% Linearization of symbolic solution for the equations of motion of the 
% n-link pendulum using the Lagrange method.
%
% Input: n - number of links
%
% Creates two files, nmodeM_[n].m and nmodeK_[n].m, which are the [M] and
% [K] matrices such that eig(M\K) returns the normal modes of oscillations
% of the system.
% ------------------------------------------------------------------------%

function nmode = func(n)

    % Symbolic variables
    syms g;
    m = sym('m', [n 1]);
    L = sym('L', [n 1]);
    d = sym('d', [n 1]);
    I = sym('I', [n 1]);
    theta = sym('theta', [n 1], 'real');
    omega = sym('omega', [n 1], 'real');
    alpha = sym('alpha', [n 1], 'real');

    % Cartesian unit vectors
    y = [0 1 0];
    z = [0 0 1];
    % Position vectors
    rOG = sym(zeros(n, 3));
    rEi = sym(zeros(n+1, 3));
    rEi(1,:) = zeros(1, 3);
    % Velocity vectors
    vOG = sym(zeros(n, 3));
    vEi = zeros(1, 3);
    % Kinetic and potential energy
    Ek = 0;
    Ep = 0;

    % Loop through all links
    for i=1:n
        
        eri = [sin(theta(i)) -cos(theta(i)) 0]; % Unit vector parallel to link i
        
        % "One man's 'elbow below' is another man's 'elbow above'."
        % - Joseph-Louis Lagrange, probably
        rEG        = d(i)*eri; % Distance from elbow above to CoM
        rOG(i,:)   = rEi(i,:) + rEG; % Distance from origin to CoM
        rEi(i+1,:) = rEi(i,:) + L(i)*eri; % Distance from origin to elbow below

        vEG      = cross(omega(i)*z, rEG); % Velocity of CoM w.r.t. elbow above
        vOG(i,:) = vEi + vEG; % Velocity of CoM w.r.t. origin
        vEi      = vEi + cross(omega(i)*z, L(i)*eri); % Velocity of elbow below w.r.t. origin

        Ek = Ek + m(i)/2*dot(vOG(i,:), vOG(i,:)) + I(i)/2*omega(i)^2; % 1/2 m v^2 + 1/2 I w^2
        Ep = Ep + m(i)*g*dot(rOG(i,:),y); % m g h
        
    end

    % Lagrangian
    Lg = Ek - Ep;

    % Lagrange's Equation, no non-conservative forces
    dLdomega   = jacobian(Lg, omega'); % dL/dw
    dLdomegadt = jacobian(dLdomega, theta')*omega ... % d/dt(dL/dw)
               + jacobian(dLdomega, omega')*alpha;    % = d/dth(dL/dw)*w + d/dw(dL/dw)*a
    dLdtheta   = jacobian(Lg, theta')'; % dL/dth
    Leq        = dLdomegadt == dLdtheta;
    
    % Linearize about theta_i ~ 0
    for i=1:n
        Leq = subs(Leq, sin(theta(i)), theta(i));
        Leq = subs(Leq, cos(theta(i)), 1);
        Leq = subs(Leq, theta(i)^2, 0);
        for j=1:n
            Leq = subs(Leq, theta(i)*theta(j), 0);
            Leq = subs(Leq, sin(theta(i) - theta(j)), theta(i) - theta(j));
            Leq = subs(Leq, cos(theta(i) - theta(j)), 1);
        end
    end
    
    % Linearized mass and stiffness matrices
    M = equationsToMatrix(Leq, alpha');
    K = equationsToMatrix(Leq, theta');
    matlabFunction(M, 'file', ['nmode/nmodeM_', num2str(n)], 'Vars', {[I] [L] [d] g [m] [omega] [theta]});
    matlabFunction(K, 'file', ['nmode/nmodeK_', num2str(n)], 'Vars', {[I] [L] [d] g [m] [omega] [theta]});
end