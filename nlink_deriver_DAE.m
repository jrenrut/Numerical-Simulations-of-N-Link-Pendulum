%-------------------------------------------------------------------------%
% nlink_deriver_DAE - Jeremy Turner
% 
% Symbolic solution for the equations of motion of the n-bar linkage using 
% the Differential Alebraic Equations method.
%
% Within the solver, n = n - 1 such that the system is treated like a
% pendulum with one end fixed.
%
% Input: n = number of links
%
% Creates two files, nlink_alphas_DAE_M_[n].m and
% nlink_alphas_Lagrange_DAE_[n], which are the [M] and [b] matrices such
% that the solutions [x] = [M]^-1 [b].
%
% Extra: Restoring force can be applied to the first link. It provides a
% very small driving force in the direction of rotation (opposite of
% friction). The intended purpose is to allow certain systems (e.g. the
% four bar crank-rocker) to make complete rotations. This should not be
% applied over long integration times as it continuously accelerates the
% system.
%
% Extra: The end can accelerate sinusoidally in any direction, with any
% freuency, and with any amplitude. This breaks the constraint of the nth
% bar, so it is really a way to drive the end of an n-link pendulum.
% ------------------------------------------------------------------------%

function nlink_deriver_DAE = func(n)

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
    Fx = sym('Fx', [n+1 1]);
    Fy = sym('Fy', [n+1 1]);
    syms rest
    d_a = sym('d_a', [1 3], 'real');
    
    % Cartesian unit vector
    k = [0 0 1];
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
        eti = [cos(theta(i)) sin(theta(i))  0]; % Unit vector in direction of link i's rotation
                
        rEG      = d(i)*eri; % Distance of elbow above w.r.t. CoM
        rBG      = (L(i) - d(i))*eri; % Distance of elbow below w.r.t CoM
        rOG(i,:) = rEi(i,:) + rEG; % Distance of CoM w.r.t. origin
        rEi(i+1,:) = rEi(i,:) + L(i)*eri; % Distance of elbow below w.r.t. origin

        aEG      = -rEG*omega(i)^2 + cross(alpha(i)*k, rEG); % Acceleration of elbow above w.r.t. CoM
        aET(i,:) = [ddx(i), ddy(i), 0] - aEG; % Acceleration of elbow above
        aOG(i,:) = aEi + aEG; % Acceleration of CoM w.r.t. origin
        aEi      = aEi - L(i)*eri*omega(i)^2 + cross(alpha(i)*k, L(i)*eri); % Acceleration of elbow below w.r.t. origin
        aEGi     = aEi - aOG(i,:); % Acceleration of elbow below w.r.t. CoM
        aEB(i,:) = [ddx(i), ddy(i), 0] + aEGi; % Acceleration of elbow below
        
        % See general case (i<n) for detail
        if i==1 % First link
            Fr = -rest*0.005*(omega(i))*eti; % Restoring force, calibrated to allow the crank-rocker system to make a complete rotation
            
            if n==1 % Only one link - no elbow below
                LMBx = m(i)*ddx(i) == Fx(i);
                LMBy = m(i)*ddy(i) == -m(i)*g + Fy(i);
                M = cross(-rEG, [Fx(i) Fy(i) 0]);
            else
                LMBx = m(i)*ddx(i) == Fx(i) - Fx(i+1);
                LMBy = m(i)*ddy(i) == -m(i)*g + Fy(i) - Fy(i+1);
                M = cross(-rEG, [Fx(i) Fy(i) 0]) + cross(rBG, -[Fx(i+1) Fy(i+1) 0]) + cross(-rEG, Fr);
            end
            
            dH = I(i)*alpha(i)*k;
            AMB = M(3) == dH(3);
            
            C = aET(i,:) == zeros(1, 3); % The base must remain stationary
            Cx = C(1);
            Cy = C(2);
            C = [Cx, Cy];
        elseif i<n % Middle case - most general
            LMBx = m(i)*ddx(i) == Fx(i) - Fx(i+1); % Only forces in x are reactions at elbow above (i) and elbow below (i+1)
            LMBy = m(i)*ddy(i) == -m(i)*g + Fy(i) - Fy(i+1); % Forces in y are gravity and reactions at elbow above (i) and elbow below (i+1)
            
            M = cross(-rEG, [Fx(i) Fy(i) 0]) + cross(rBG, -[Fx(i+1) Fy(i+1) 0]); % Torque about CoM from reaction forces at elbows above and below
            dH = I(i)*alpha(i)*k; % Change in angular momentum at CoM = I a
            AMB = M(3) == dH(3); % All terms are in z direction
            
            C = aET(i,:) == aEB(i-1,:); % Constraint: acceleration of elbow above must be equal to
            Cx = C(1);                  % the acceleration of the elbow below of the previous link
            Cy = C(2);
            C = [Cx, Cy];
        else % Last link - still constrained below, so it has the same equations as the general case
            LMBx = m(i)*ddx(i) == Fx(i) - Fx(i+1);
            LMBy = m(i)*ddy(i) == -m(i)*g + Fy(i) - Fy(i+1);
            
            M = cross(-rEG, [Fx(i) Fy(i) 0]) + cross(rBG, -[Fx(i+1) Fy(i+1) 0]);
            dH = I(i)*alpha(i)*k;
            AMB = M(3) == dH(3);
            
            C1 = aET(i,:) == aEB(i-1,:);
            Cx1 = C1(1);
            Cy1 = C1(2);
            C2 = aEB(i,:) == -d_a; % User specified acceleration of end
            Cx2 = C2(1);
            Cy2 = C2(2);
            C = [Cx1, Cy1, Cx2, Cy2];
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
    matlabFunction(M, 'file', ['DAElink/nlink_alphas_DAE_M_', num2str(n+1)], 'Vars', {[I] [L] [d] g [m] [omega] [theta] rest [d_a]});
    matlabFunction(b, 'file', ['DAElink/nlink_alphas_DAE_b_', num2str(n+1)], 'Vars', {[I] [L] [d] g [m] [omega] [theta] rest [d_a]});
    
end