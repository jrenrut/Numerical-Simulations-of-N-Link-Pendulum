%-------------------------------------------------------------------------%
% comp_results - Jeremy Turner
%
% This file was written solely to compare the time to divergence of a
% certain set of files.
% ------------------------------------------------------------------------%

clear all

n = 10;
m = 5;

for i=1:n
    fs{i,1} = ['animations/D10_', num2str(i), 'o1'];
    for j=2:m
        fs{i,j} = ['animations/D10_', num2str(i), 'o', num2str(j), '1'];
    end
end

load('animations/D10_30s')
t1 = t;
x1 = xs(:,2:end);
y1 = ys(:,2:end);

tm = zeros(n, m);

for i=1:n
    for j=1:m
        
        f = fs{i,j};
        load(f)

        t2 = t;
        x2 = xs(:,2:end);
        y2 = ys(:,2:end);

        xs = zeros(length(t), n);
        ys = zeros(length(t), n);

        Ltot = sqrt(x1(1,end)^2 + y1(1,end)^2);
        L = Ltot/n;

        xe = zeros(length(t), n);
        ye = zeros(length(t), n);
        for k=1:n
            xe(:,k) = x1(:,k) - x2(:,k);
            ye(:,k) = y1(:,k) - y2(:,k);
        end

        tc = 0;
        k = 1;
        while tc < 0.01*L
            diffs = sqrt(xe(k,:).^2 + ye(k,:).^2);
            tc = max(diffs);
            k = k + 1;
        end
        k = k - 1;
        
        fprintf(['Link ', num2str(i),', Offset = 1e-', num2str(j), ': t = ', num2str(t(k)), ' seconds.\n']);
        
        tm(i, j) = t(k);
        
    end
end

figure
hold on
xlabel('Offset in \theta_0'); ylabel('t [s]');
title('Time to Divergence vs. Initial \theta_0 Offset');
for i=1:n
    plot([1e-5 1e-4 1e-3 1e-2 1e-1], fliplr(tm(i,:)))
end
set(gca, 'XScale', 'log');
names = {};
for i=1:n
    names{i} = ['Link ', num2str(i)];
end
legend(names)