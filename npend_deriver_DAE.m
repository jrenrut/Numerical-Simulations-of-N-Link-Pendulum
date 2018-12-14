%-------------------------------------------------------------------------%
% npend_deriver_DAE - Jeremy Turner
% 
% Symbolic solution for the equations of motion of the n-link pendulum
% using the Differential Alebraic Equations method.
%
% Input: n = number of links
%
% Creates two files, npend_alphas_DAE_M_[n].m and
% npend_alphas_Lagrange_DAE_[n], which are the [M] and [b] matrices such
% that the solutions [x] = [M]^-1 [b].
%
% Extra: The base can accelerate sinusoidally in any direction, with any
% freuency, and with any amplitude.
% ------------------------------------------------------------------------%

function npend_deriver_DAE = func(n)

    % Symbolic variables
    syms g;
    m = sym('m', [n 1]);
    L = sym('L', [n 1]);
    d = sym('d', [n 1]);
    I = sym('I', [n 1]);
    ddx = sym('ddx', [n 1]);
    ddy = sym('ddy', [n 1]);
    theta = sym('theta', [n 1], 'real');
    omega = sym('omega', [n 1], 'real');
    alpha = sym('alpha', [n 1], 'real');
    Fx = sym('Fx', [n 1]);
    Fy = sym('Fy', [n 1]);
    d_a = sym('d_a', [1 3], 'real');
    
    % Cartesian unit vector
    z = [0 0 1];
    % Position vectors
    rOG = sym(zeros(n, 3));
    rEi = sym(zeros(n+1, 3));
    % Acceleration vectors
    aOG = sym(zeros(n, 3));
    aEi = zeros(1, 3);
    aET = sym(zeros(n, 3));
    aEB = sym(zeros(n, 3));
    % Equations
    LMBs = [];
    AMBs = [];
    Cs   = [];
    
    % Loop through each link
    for i=1:n
        
        eri = [sin(theta(i)) -cos(theta(i)) 0]; % Unit vector parallel to link i
        
        rEG      = d(i)*eri; % Distance from elbow above to CoM
        rBG      = (L(i) - d(i))*eri; % Distance from CoM to elbow below
        rOG(i,:) = rEi(i,:) + rEG; % Distance from origin to CoM
        rEi(i+1,:) = rEi(i,:) + L(i)*eri; %  Distance from elbow below to origin
        
        aEG      = -rEG*omega(i)^2 + cross(alpha(i)*z, rEG); % Acceleration of elbow above w.r.t. CoM
        aET(i,:) = [ddx(i), ddy(i), 0] - aEG; % Acceleration of elbow above
        aOG(i,:) = aEi + aEG; % Acceleration of CoM w.r.t. origin
        aEi      = aEi - L(i)*eri*omega(i)^2 + cross(alpha(i)*z, L(i)*eri); % Acceleration of elbow below w.r.t. origin
        aEGi     = aEi - aOG(i,:); % Acceleration of elbow below w.r.t. CoM
        aEB(i,:) = [ddx(i), ddy(i), 0] + aEGi; % Acceleration of elbow below
        
        % See general case (i<n) for detail
        if i==1 % Fist link
            if n==1 % Only one link - no elbow below
                LMBx = m(i)*ddx(i) == Fx(i);
                LMBy = m(i)*ddy(i) == -m(i)*g + Fy(i);
                M = cross(-rEG, [Fx(i) Fy(i) 0]);
            else
                LMBx = m(i)*ddx(i) == Fx(i) - Fx(i+1);
                LMBy = m(i)*ddy(i) == -m(i)*g + Fy(i) - Fy(i+1);
                M = cross(-rEG, [Fx(i) Fy(i) 0]) + cross(rBG, -[Fx(i+1) Fy(i+1) 0]);
            end         
            
            dH = I(i)*alpha(i)*z;
            AMB = M(3) == dH(3);
            
            C = aET(i,:) == d_a; % User specified acceleration of base
            Cx = C(1);
            Cy = C(2);
            C = [Cx, Cy];
        elseif i<n % Middle case - most general
            LMBx = m(i)*ddx(i) == Fx(i) - Fx(i+1); % Only forces in x are reactions at elbow above (i) and elbow below (i+1)
            LMBy = m(i)*ddy(i) == -m(i)*g + Fy(i) - Fy(i+1); % Forces in y are gravity and reactions at elbow above (i) and elbow below (i+1)
            
            M = cross(-rEG, [Fx(i) Fy(i) 0]) + cross(rBG, -[Fx(i+1) Fy(i+1) 0]); % Torque about CoM from reaction forces at elbows above and below
            dH = I(i)*alpha(i)*z; % Change in angular momentum at CoM = I a
            AMB = M(3) == dH(3); % All terms are in z direction
            
            C = aET(i,:) == aEB(i-1,:); % Constraint: acceleration of elbow above must be equal to
            Cx = C(1);                  % the acceleration of the elbow below of the previous link
            Cy = C(2);
            C = [Cx, Cy];
        else % Last link - no elbow below
            LMBx = m(i)*ddx(i) == Fx(i);
            LMBy = m(i)*ddy(i) == -m(i)*g + Fy(i);
            
            M = cross(-rEG, [Fx(i) Fy(i) 0]);
            dH = I(i)*alpha(i)*z;
            AMB = M(3) == dH(3);
            
            C = aET(i,:) == aEB(i-1,:);
            Cx = C(1);
            Cy = C(2);
            C = [Cx, Cy];
        end
        
        % Add to linear and angular momentum balance and constraint equations
        LMBs = [LMBs, LMBx, LMBy];
        AMBs = [AMBs, AMB];
        Cs   = [Cs, C];

    end
    
    % All equations and symbolic variables
    eqns = [LMBs, AMBs, Cs];
    vars = [ddx.', ddy.', alpha', Fx.', Fy.'];
    
    % Differential equations for theta in matrix form
    [M, b] = equationsToMatrix(eqns, vars);
    matlabFunction(M, 'file', ['DAEpend/npend_alphas_DAE_M_', num2str(n)], 'Vars', {[I] [L] [d] g [m] [omega] [theta] [d_a]});
    matlabFunction(b, 'file', ['DAEpend/npend_alphas_DAE_b_', num2str(n)], 'Vars', {[I] [L] [d] g [m] [omega] [theta] [d_a]});
    
end