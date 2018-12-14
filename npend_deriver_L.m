%-------------------------------------------------------------------------%
% npend_deriver_L - Jeremy Turner
% 
% Symbolic solution for the equations of motion of the n-link pendulum
% using the Lagrange method.
%
% Input: n - number of links
%
% Creates two files, npend_alphas_Lagrange_M_[n].m and
% npend_alphas_Lagrange_b_[n], which are the [M] and [b] matrices such that
% the solutions [x] = [M]^-1 [b].
% ------------------------------------------------------------------------%

function npend_deriver_L = func(n)

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
    
    % Differential equations for theta in matrix form
    [M, b] = equationsToMatrix(Leq, alpha');
    matlabFunction(M, 'file', ['Lpend/npend_alphas_Lagrange_M_', num2str(n)], 'Vars', {[I] [L] [d] g [m] [omega] [theta]});
    matlabFunction(b, 'file', ['Lpend/npend_alphas_Lagrange_b_', num2str(n)], 'Vars', {[I] [L] [d] g [m] [omega] [theta]});  
    
end