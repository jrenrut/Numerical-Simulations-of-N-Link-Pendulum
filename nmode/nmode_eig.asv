function nmode_eig = func(p, Mf, Kf)

    I = 1/12*ones(n, 1);
    L = ones(n-1, 1);
    d = 0.5*ones(n, 1);
    g = 1;
    m = ones(n, 1);
    omega = zeros(n, 1);
    theta = 0*10*pi/180*ones(n, 1);

    Mf = ['nmodeM_', num2str(n)];
    Kf = ['nmodeK_', num2str(n)];
    M = feval(Mf, I, L, d, g, m, omega, theta);
    K = feval(Kf, I, L, d, g, m, omega, theta);

    [v, d] = eig(M\K);
    
    nmode_eig = [v, d];
end