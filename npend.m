%-------------------------------------------------------------------------%
% npend - Jeremy Turner
% 
% This is the file that drives all the magic. By running npend.m, the user
% is promped on the command line to:
% - select to play an animation
% --- list all available animations
% --- play animations
% --- replay animations
% --- quit
% - create a simulations
% --- enter number of links
% --- choose solution method
% --- specify parameters and initial conditions
% ----- value can be entered, or option selected to assign values
% ------- randomized given a mean and standard deviation
% ------- specified individually
% --- replay animation
% --- save animation
% --- rerun the simulation with new parameters and initial conditions
% ------------------------------------------------------------------------%

clear all
fprintf('===========================================================================\n');

path = pwd;

% ------------------------------------------------------------------------%
% Load animations
% ------------------------------------------------------------------------%

lfq = '';
while isempty(lfq) == 1
    lfq = input('> Load animation? ', 's');
end
if lower(lfq(1)) == 'y'
    if exist('animations', 'dir') ~= 7
        fprintf(' !> There are no saved animations, go make some! <!\n');
    else
        lfn = '';
        while isempty(lfn) == 1
            lfn = input(' > File Name or [L]ist? ', 's');
            if lower(lfn(1)) == 'l'
                files = dir('animations/*.mat');
                flen = 0;
                for i=1:length(files)
                    fname = files(i).name(1:end-4);
                    flen = flen + length(fname) + 4;
                    if flen > 200
                        flen = 0;
                        fprintf('\n');
                    end
                    fprintf(['\t', fname]);
                end
                fprintf('\n');
                lfn = '';
            elseif isfile(['animations/', lfn, '.mat']) == 0
                fprintf(' !> File does not exist, try again.\n');
                lfn = '';
            else
                rp = '';
                while isempty(rp) == 1
                    animate(['animations/', lfn, '.mat'])
                    rp = input('> [R]eplay or [O]ther or [Q]uit? ', 's');
                    if lower(rp(1)) == 'r'
                        rp = '';
                    elseif lower(rp(1)) == 'o'
                        lfn = '';
                    else
                        return
                    end
                end
            end
        end
    end
end          

% ------------------------------------------------------------------------%
% Specify number of links and simulation method
% ------------------------------------------------------------------------%

n = '';
while isempty(n) == 1
    n = input('> Number of Links: ');
end

ptype = '';
while isempty(ptype) == 1
    ptype = lower(input('> [N]ewton-Euler, [L]agrange, [D]AE, or [M]odes? ', 's'));
end
ptype = ptype(1);

prob = '';
if ptype == 'n'
    if exist('Npend', 'dir') ~= 7
        mkdir('Npend');
        tic
        npend_deriver_NE(n)
        toc
    else
        if isfile(['Npend/npend_alphas_NE_M_', num2str(n), '.m']) == 0
            tic
            npend_deriver_NE(n)
            toc
        end
    end
elseif ptype == 'l'
    if exist('Lpend', 'dir') ~= 7
        mkdir('Lpend');
        tic
        npend_deriver_L(n)
        toc
    else
        if isfile(['Lpend/npend_alphas_Lagrange_M_', num2str(n), '.m']) == 0
            tic
            npend_deriver_L(n)
            toc
        end
    end
elseif ptype == 'd'
    while isempty(prob) == 1
        prob = lower(input('> [P]endulum or [L]inkage ([i]nfo on linkage set-up)? ', 's'));
    end
    if prob(1) == 'i'
        fprintf('> A GUI will pop up, you will be prompted to enter the end position in x and y,\n')
        fprintf('> and for the lengths and intial angles of the first n-2 links. You may enter\n')
        fprintf('> them one at a time separated by spaces, or once and for all. The n-1 link will\n')
        fprintf('> be automatically generated to meet the constraints. The nth link is the non-\n')
        fprintf('> base\n')
        prob = '';
        while isempty(prob) == 1
            prob = lower(input('> [P]endulum or [L]inkage? ', 's'));
        end
        prob = prob(1);
    end
    if prob == 'p'
        if exist('DAEpend', 'dir') ~= 7
            mkdir('DAEpend');
            tic
            npend_deriver_DAE(n)
            toc
        else
            if isfile(['DAEpend/npend_alphas_DAE_M_', num2str(n), '.m']) == 0
                tic
                npend_deriver_DAE(n)
                toc
            end
        end
    else
        n = n - 1;
        if exist('DAElink', 'dir') ~= 7
            mkdir('DAElink');
            tic
            nlink_deriver_DAE(n)
            toc
        else
            if isfile(['DAElink/nlink_alphas_DAE_M_', num2str(n), '.m']) == 0
                tic
                nlink_deriver_DAE(n)
                toc
            end
        end
    end
elseif ptype == 'm'
    if exist('nmode', 'dir') ~= 7
        mkdir('nmode');
        tic
        nmode(n)
        toc
    else
        if isfile(['nmode/nmodeM_', num2str(n), '.m']) == 0
            tic
            nmode(n)
            toc
        end
    end
    if exist('Lpend', 'dir') ~= 7
        mkdir('Lpend');
        tic
        npend_deriver_L(n)
        toc
    else
        if isfile(['Lpend/npend_alphas_Lagrange_M_', num2str(n), '.m']) == 0
            tic
            npend_deriver_L(n)
            toc
        end
    end
else
    fprintf('!> Error: You are too stupid to follow the instructions <!\n');
    return
end

% ------------------------------------------------------------------------%
% Enter parameters and initial conditions
% ------------------------------------------------------------------------%

done = 0;
while done == 0
    
% ------------------------------------------------------------------------%
% Normal mode
% ------------------------------------------------------------------------%
    
    if ptype == 'm'
        mode = '';
        while isempty(mode) == 1
            mode = input(['> Normal Mode (1 -> ', num2str(n), '): ']);
        end
    end

% ------------------------------------------------------------------------%
% Friction for Newton-Euler
% ------------------------------------------------------------------------%
    
    if ptype ~= 'd' || prob ~= 'l'
        
        if ptype == 'n'
           Fq = '';
           while isempty(Fq) == 1
               Fq = input('> Friction ([N]o or value)? ', 's');
           end
           if lower(Fq) == 'n'
               p.fric = 0;
           else
               p.fric = str2num(Fq);
           end
        end

% ------------------------------------------------------------------------%
% Link lengths
% ------------------------------------------------------------------------%
        
        Lq = '';
        while isempty(Lq) == 1
            Lq = input('> Link Length or [R]andom or [S]pecified? ', 's');
        end
        [Lval isnum] = str2num(Lq);
        if isnum == 1
            L = Lval.*ones(n, 1);
        elseif lower(Lq) == 'r'
            Lmn = '';
            while isempty(Lmn) == 1
                Lmn = input(' > Link Length Mean: ');
            end
            Lsd = '';
            while isempty(Lsd) == 1
                Lsd = input(' > Link Length Standard Deviation: ');
            end
            L = Lsd.*randn(n, 1) + Lmn;
        elseif lower(Lq) == 's'
            L = ones(n, 1);
            for i=1:n
                Li = '';
                while isempty(Li) == 1
                    Li = input([' > Link ', num2str(i), ' Length: ']);
                end
                L(i) = Li;
            end
        else
            fprintf('!> Error: You are too stupid to follow the instructions <!\n');
            return
        end

% ------------------------------------------------------------------------%
% Initial angles
% ------------------------------------------------------------------------%
        
        if ptype ~= 'm'
            
            theta0 = '';
            while isempty(theta0) == 1
                theta0 = input('> Uniform Initial Angles or [R]andom or [S]pecified? ', 's');
            end
            [tval isnum] = str2num(theta0);
            if isnum == 1
                theta0 = tval*pi/180.*ones(n, 1);
            elseif lower(theta0) == 'r'
                tmn = '';
                while isempty(tmn) == 1
                    tmn = input(' > Initial Angles Mean: ');
                end
                tsd = '';
                while isempty(tsd) == 1
                    tsd = input(' > Initial Angles Standard Deviation: ');
                end
                theta0 = (tsd.*randn(n, 1) + tmn).*pi/180;
            elseif lower(theta0) == 's'
                theta0 = ones(n, 1);
                for i=1:n
                    theta0i = '';
                    while isempty(theta0i) == 1
                        theta0i = input([' > Initial Angle for Link ', num2str(i), ': ']);
                    end
                    theta0(i) = theta0i*pi/180;
                end
            end

% ------------------------------------------------------------------------%
% Initial angular velocities
% ------------------------------------------------------------------------%
        
            omega0 = '';
            while isempty(omega0) == 1
                omega0 = input('> Uniform Initial Angular Velocities or [R]andom or [S]pecified? ', 's');
            end
            [tval isnum] = str2num(omega0);
            if isnum == 1
                omega0 = tval*pi/180.*ones(n, 1);
            elseif lower(omega0) == 'r'
                omn = '';
                while isempty(omn) == 1
                    omn = input(' > Initial Angular Velocities Mean: ');
                end
                osd = '';
                while isempty(osd) == 1
                    osd = input(' > Initial Angular Velocities Standard Deviation: ');
                end
                omega0 = (osd.*randn(n, 1) + omn).*pi/180;
            elseif lower(omega0) == 's'
                omega0 = ones(n, 1);
                for i=1:n
                    omega0i = '';
                    while isempty(omega0i) == 1
                        omega0i = input([' > Initial Anglular Velocity for Link ', num2str(i), ': ']);
                    end
                    omega0(i) = omega0i*pi/180;
                end
            else
                fprintf('!> Error: You are too stupid to follow the instructions <!\n');
                return
            end
                                   
        end
        
        xend = 0;
        yend = 0;
        rq = 0;
        
    else

% ------------------------------------------------------------------------%
% N-bar linkage parameters and initial conditions
% ------------------------------------------------------------------------%
        
        title = [num2str(n), ' Bar Linkage Input (1 -> n-1)'];
        prompt = {'End offset1 x', 'End offset1 y', ['Fisrt ', num2str(n-1), ' Link Lengths (separated by spaces)'], ['Fisrt ', num2str(n-1), ' Initial \theta s (separated by spaces)']};
        definput = {num2str(n-1), num2str(n-1)};
        lin = '';
        tin = '';
        dims = [1 30; 1 30; 1 30; 1 30];
        for i=1:n-1
            lin = [lin, num2str(sqrt(2)), ' '];
            th = 180 - 45/(n-1)*i;
            tin = [tin, [num2str(th), ' ']];
        end
        definput{end + 1} = lin;
        definput{end + 1} = tin;
        opts.Interpreter = 'tex';
        data = inputdlg(prompt, title, dims, definput, opts);
        if isempty(data) == 1
            fprintf('!> Error: Cancelling is for quitters. <!\n');
            return
        end
        xend = str2num(data{1});
        yend = str2num(data{2});
        Lc = strsplit(data{3}, ' ');
        theta0c = strsplit(data{4}, ' ');
        L = zeros(n, 1);
        theta0 = zeros(n, 1);
        xn = -yend;
        yn = xend;
        for i=1:n-1
            if length(Lc) == 1
                Li = str2num(Lc{1});
            else
                Li = str2num(Lc{i});
            end
            if length(theta0c) == 1
                theta0i = str2num(theta0c{1})*pi/180;
            else
                theta0i = str2num(theta0c{i})*pi/180;
            end
            L(i) = Li;
            theta0(i) = theta0i;
            xn = xn - Li*cos(theta0i);
            yn = yn - Li*sin(theta0i);
        end
        L(n) = sqrt(xn^2 + yn^2);
        theta0(n) = atan2(yn, xn);
        omega0 = zeros(n, 1);
        rq = '';
        while isempty(rq) == 1
            rq = input('> Restoring force ([N]o or value)? ', 's');
        end
        if lower(rq) == 'n'
            p.rest = 0;
            rq = 0;
        else
            p.rest = str2num(rq);
        end
        
    end

% ------------------------------------------------------------------------%
% Link masses
% ------------------------------------------------------------------------%

        mq = '';
        while isempty(mq) == 1
            mq = input('> Mass of Links or [P]roportional or [R]andom or [S]pecified? ', 's');
        end
        [mval isnum] = str2num(mq);
        if isnum == 1
            m = mval.*ones(n, 1);
        elseif lower(mq) == 'p'
            m = L.*ones(n, 1);
        elseif lower(mq) == 'r'
            mmn = '';
            while isempty(mmn) == 1
                mmn = input(' > Link Mass Mean: ');
            end
            msd = '';
            while isempty(msd) == 1
                msd = input(' > Link Mass Standard Deviation: ');
            end
            m = msd.*randn(n, 1) + mmn;
        elseif lower(mq) == 's'
            m = ones(n, 1);
            for i=1:n
                mi = '';
                while isempty(mi) == 1
                    mi = input([' > Mass of Link ', num2str(i), ': ']);
                end
                m(i) = mi;
            end
        else
            fprintf('!> Error: You are too stupid to follow the instructions <!\n');
            return
        end

% ------------------------------------------------------------------------%
% Fractional distances to centers of link masses
% ------------------------------------------------------------------------%

        dq = '';
        while isempty(dq) == 1
            dq = input('> Fractional Distance to CoM or [R]andom or [S]pecified? ', 's');
        end
        [dval isnum] = str2num(dq);
        d = ones(n, 1);
        if isnum == 1
            for i=1:n
                d(i) = dval*L(i);
            end
        elseif lower(dq) == 'r'
            dmn = '';
            while isempty(dmn) == 1
                dmn = input(' > Fractional Distance to CoM Mean: ');
            end
            dsd = '';
            while isempty(dsd) == 1
                dsd = input(' > Fractional Distance to COM Standard Deviation: ');
            end
            dr = dsd.*randn(n, 1) + dmn;
            for i=1:n
                d(i) = dr(i)*L(i);
            end
        elseif lower(dq) == 's'
            for i=1:n
                ds = '';
                while isempty(ds) == 1
                    ds = input([' > Fractional Distance to CoM ', num2str(i), ': ']);
                end
                d(i) = ds*L(i);
            end
        else
            fprintf('!> Error: You are too stupid to follow the instructions <!\n');
            return
        end

% ------------------------------------------------------------------------%
% Moment of inertia of links
% ------------------------------------------------------------------------%

        I = ones(n, 1);
        for i=1:n
            I(i) = m(i)/6*(L(i)^2 - 2*L(i)*d(i) + 2*d(i)^2);
        end

% ------------------------------------------------------------------------%
% Shaking base / end for DAE pendulum / linkage
% ------------------------------------------------------------------------%
    
    if ptype == 'd'
        dst = '';
        while isempty(dst) == 1
            if prob == 'p'
                dst = input('> Shaking Base (or Standing Wave [M]odes)? ', 's');
            else
                dst = input('> Shaking End? ', 's');
            end
        end
        if lower(dst(1)) == 'y' || lower(dst(1)) == 'm'
            if lower(dst(1)) == 'm'
                nmode = '';
                while isempty(nmode) == 1
                    nmode = input(' > Standing Wave Mode: ');
                end
                bj0 = @(x) besselj(0, x);
                J0 = fzero(bj0, [(nmode-1) nmode]*pi);
                wf = 0; % Calculated after gravity is defined
            else
                wf = '';
                while isempty(wf) == 1
                    wf = input(' > Driving Frequency: ');
                end
            end
            offset = '';
            while isempty(offset) == 1
                offset = input(' > Disturbance Amplitude: ');
            end
            phi = '';
            while isempty(phi) == 1
                phi = pi/180*input(' > Angle of Disturbance: ');
            end
        else
            offset = 0;
            wf = 0;
            phi = 0;
        end
        p.offset = offset;
        p.wf = wf;
        p.phi = phi;
        d_r = zeros(3, 1);
    else
        dst = 'n';
    end

% ------------------------------------------------------------------------%
% Acceleration due to gravity
% ------------------------------------------------------------------------%
    
    g = '';
    while isempty(g) == 1
        g = input('> Acceleration due to Gravity: ');
    end
    if ptype == 'd' && lower(dst(1)) == 'm'
        p.wf = J0/2*sqrt(g/sum(L));
    end

% ------------------------------------------------------------------------%
% Timespan and tolerance
% ------------------------------------------------------------------------%

    tq = '';
    while isempty(tq) == 1
        tq = input('> Runtime (or [t]ime/tolerance options): ', 's');
    end
    if lower(tq(1)) == 't'
        tmax = '';
        while isempty(tmax) == 1
            tmax = input(' > Runtime: ');
        end
        tstep = '';
        while isempty(tstep) == 1
            tstep = input(' > Number of timesteps: ');
        end
        tol = '';
        while isempty(tol) == 1
            tol = input(' > Absolute/Relative Tolerance: ');
        end
    else
        tmax = str2num(tq);
        tstep = 10*tmax;
        tol = 1e-8;
    end

% ------------------------------------------------------------------------%
% Call ODE solvers
% ------------------------------------------------------------------------%
    
    p.g = g; p.L = L; p.d = d; p.m = m; p.I = I; p.n = n;
    
    tspan = linspace(0, tmax, tstep);
    
    tol = 1e-8;
    opts = odeset('RelTol', tol, 'AbsTol', tol);
    
    if ptype == 'm'
        Mf = ['nmode/nmodeM_', num2str(n)];
        Kf = ['nmode/nmodeK_', num2str(n)];
        z = zeros(2*n, 1);
        cd('nmode');
        eigs = nmode_eig(p, Mf, Kf);
        cd('../');
        theta0 = eigs(:,mode)/10;
        omega0 = zeros(n,1);
        freq = sqrt(sum(eigs(:,n+mode)));
    elseif lower(dst(1)) == 'm'
        freq = p.wf;
        mode = nmode;
    else
        freq = 0;
        mode = 0;
    end
        
    
    z0 = [theta0; omega0];
    
    if ptype == 'n'
        Mf = ['npend_alphas_NE_M_', num2str(n)];
        bf = ['npend_alphas_NE_b_', num2str(n)];
        cd('Npend');
        f = @(t, z) npend_NE(z, p, Mf, bf);
        [t, z] = ode45(f, tspan, z0, opts);
        cd('../');
        [s_y, s_F] = audioread('bioh2.mp3');
    elseif ptype == 'l' || ptype == 'm'
        Mf = ['npend_alphas_Lagrange_M_', num2str(n)];
        bf = ['npend_alphas_Lagrange_b_', num2str(n)];
        cd('Lpend');
        f = @(t, z) npend_Lagrange(z, p, Mf, bf);
        [t, z] = ode45(f, tspan, z0, opts);
        cd('../');
        [s_y, s_F] = audioread('wtlb.mp3');
    elseif prob == 'p'
        Mf = ['npend_alphas_DAE_M_', num2str(n)];
        bf = ['npend_alphas_DAE_b_', num2str(n)];
        z0 = [z0; d_r];
        cd('DAEpend');
        f = @(t, z) npend_DAE(z, p, t, Mf, bf);
        [t, z] = ode45(f, tspan, z0, opts);
        cd('../');
        [s_y, s_F] = audioread('md.mp3');
    else
        Mf = ['nlink_alphas_DAE_M_', num2str(n+1)];
        bf = ['nlink_alphas_DAE_b_', num2str(n+1)];
        p.xend = xend; p.yend = yend;
        cd('DAElink');
        f = @(t, z) nlink_DAE(z, p, t, Mf, bf);
        [t, z] = ode45(f, tspan, z0, opts);
        cd('../');
        [s_y, s_F] = audioread('o.mp3');
    end

% ------------------------------------------------------------------------%
% Calculate link positions from solutions
% ------------------------------------------------------------------------%
        
    xs = ones(length(t), n+1);
    ys = ones(length(t), n+1);
    xs(:,1) = zeros(length(t), 1);
    ys(:,1) = zeros(length(t), 1);
    
    if ptype == 'd' && prob == 'p'
        xs(:,1) = z(:,2*n+1);
        ys(:,1) = z(:,2*n+2);
    end
        
    xGs = ones(length(t), n);
    yGs = ones(length(t), n);
    
    thetas = z(:,1:n);
    omegas = z(:,n+1:2*n);
    
    for i=1:n
       thetan = z(:,i);
       omegan = z(:,n+i);
       xs(:,i+1) = xs(:,i) + L(i)*sin(thetan);
       ys(:,i+1) = ys(:,i) - L(i)*cos(thetan);
    end

% ------------------------------------------------------------------------%
% Calculate total energy throughout simulation
% ------------------------------------------------------------------------%

    i = [1 0 0];
    j = [0 1 0];
    k = [0 0 1];
    
    E = zeros(length(t), 1);
    
    for x=1:length(t)
        
        rOG = zeros(n, 3);
        rEi = zeros(n+1, 3);
        
        vOG = zeros(n, 3);
        vOi = zeros(1, 3);
        
        Ek = 0;
        Ep = 0;
        
        for y=1:n
            er = [sin(thetas(x,y)) -cos(thetas(x,y)) 0];
            
            rEG      = d(y)*er;
            rOG(y,:) = rEi(y,:) + rEG;
            rEi(y+1,:) = rEi(y,:) + L(y)*er;

            vEG      = cross(omegas(x,y)*k, rEG);
            vOG(y,:) = vOi + vEG;
            vOi      = vOi + cross(omegas(x,y)*k, L(y)*er);

            Ek = Ek + m(y)/2*dot(vOG(y,:), vOG(y,:)) + I(y)/2*omegas(x,y)^2;
            Ep = Ep + m(y)*g*dot(rOG(y,:),j);
        end
        
        E(x) = Ek + Ep;
        
    end

% ------------------------------------------------------------------------%
% Animate
% ------------------------------------------------------------------------%

    show = '';
    while isempty(show) == 1

        Lnetx = sum(L) + min(L)/2 - abs(xend);
        Lnety = sum(L) + min(L)/2 - abs(yend);
        xmin = -Lnetx - (xend - abs(xend))/2;
        xmax = Lnetx + (xend + abs(xend))/2;
        ymin = -Lnety - (yend - abs(yend))/2;
        ymax = Lnety + (yend + abs(yend))/2;
        
        m = sum(L) + min(L)/2;

        figure('units', 'normalized', 'outerposition', [0.1 0.1 0.9 0.9])
        set(gca,'Color', 'k')
        hold on
        xlabel('x'); ylabel('y'); title(['n = ', num2str(n)]);
        if ptype == 'd' && prob == 'l' && offset2 == 0
            plot([0 xend], [0 yend],  'color', [0.6350, 0.0780, 0.1840], 'LineWidth', 2);
        end
        sound(s_y, s_F)
        for i=1:length(t)
            pl = plot(xs(i,:), ys(i,:), 'r', 'LineWidth', 2);
            if rq ~= 0
                if i > 1
                    plot([xs(i-1,2) xs(i,2)], [ys(i-1,2) ys(i,2)], ':w')
                end
            elseif lower(dst(1)) ~= 'm' && ptype == 'm'
                if i > 1
                    plot([xs(i-1,n+1) xs(i,n+1)], [ys(i-1,n+1) ys(i,n+1)], ':w')
                end
            end
            axis equal;
            if ptype == 'd' && prob == 'l'
                axis([xmin xmax ymin ymax])
            elseif lower(dst(1)) == 'm'
                axis([-m/2 m/2 -1.1*sum(L) max(L)]);
            else
                axis([-m m -m m])
            end
            pause(0.01)
            shg
            if i < length(t)
                set(pl, 'Visible', 'off')
            end            
        end

        clear sound
        
        figure
        plot(t, E)
        xlabel('t'); ylabel('E'); title('Total Energy');
        
        if ptype == 'm'
            figure
            hold on
            xlabel('t'); ylabel('\theta'); title(['Angles vs. Time: ', num2str(n), ' Links, Normal Mode ', num2str(mode), ', \omega = ', num2str(freq), ' rad/s']);
            legs = {};
            for i=1:n
                plot(t, thetas(:,i))
                legs{i} = ['\theta_', num2str(i)];
            end
            legend(legs)
            
            figure
            hold on
            xlabel('t'); ylabel('\theta'); title({['Linearized - Actual Angles vs. Time: ', num2str(n), ' Links, Normal Mode ', num2str(mode), ', \omega = ', num2str(freq), ' rad/s'], ''});
            legs = {};
            for i=1:n
                plot(t, thetas(:,i)-thetas(1,i)*cos(freq*t))
                legs{i} = ['\Delta \theta_', num2str(i)];
            end
            legend(legs, 'Location', 'northwest')
        elseif lower(dst(1)) == 'm'
            figure('units', 'normalized', 'outerposition', [0.2 0.4 0.8 0.6])
            hold on
            xlabel('t'); ylabel('\theta'); title(['End Angle vs. Time: ', num2str(n), ' Links, Standing Mode ', num2str(mode), ', \omega = ', num2str(freq), ' rad/s']);
            plot(t, thetas(:,n))
            plot(t, cos(freq*t))
            ylim([-1.2 1.4])
            legend('End Frequency', 'Driving Frequency');
        end
        
        sq = '';
        while isempty(sq) == 1
        	sq = input('> Replay? ', 's');
        end
        if lower(sq(1)) ~= 'y'
            show = 'no';
        else
            clf
            close
            clf
            close
        end
    end

% ------------------------------------------------------------------------%
% Save
% ------------------------------------------------------------------------%
    
    svfg = '';
    while isempty(svfg) == 1
        svfg = input('> Save animation? ', 's');
    end
    if lower(svfg(1)) == 'y'
        if exist('animations', 'dir') ~= 7
            mkdir('animations');
        end
        svnm = '';
        while isempty(svnm) == 1
            svnm = input(' > File Name: ', 's');
            if isfile(['animations/', svnm, '.mat']) == 1
                fprintf(' !> Try again, file already exists <!\n');
                svnm = '';
            end
        end
        save([path, '/animations/', svnm, '.mat'], 'n', 't', 'xs', 'ys', 'thetas', 's_y', 's_F', 'ptype', 'prob', 'L', 'xend', 'yend', 'tol', 'E', 'rq', 'dst', 'mode', 'freq')
    end

% ------------------------------------------------------------------------%
% Replay
% ------------------------------------------------------------------------%
    
    done = '';
    while isempty(done) == 1
        done = input('> Run again? ', 's');
    end
    if lower(done(1)) == 'y'
        done = 0;
        clf
        close
        clf
        close
    else
        done = 1;
    end

end