%-------------------------------------------------------------------------%
% comp - Jeremy Turner
% 
% Compares either 2 or 3 existing simulations. The simulations must run for
% the same amount of time / have the same number of timesteps. If 3
% simulations are entered, they will be animated together, and their total
% energies will be plotted together on the same figure. If 2 simulations
% are entered, an additional plot will show when the two simulations begin
% to diverge. This is intended for similar simulations. Divergence is
% defined as when any pair the end points of the links between the pendula
% become separated by more than 0.001 of the average link length.
% ------------------------------------------------------------------------%

clear all

n = '';
while isempty(n) == 1
    n = input('> Number of Links: ');
end

f1 = '';
while isempty(f1) == 1
    f1 = input('> First simulation ([L]ist): ', 's');
    if f1 == 'l'
        if exist('animations', 'dir') ~= 7
            fprintf(' !> There are no saved simulations, go make some! <!\n');
            return
        else
            files = dir('animations/*.mat');
            for i=1:length(files)
                fprintf(['\t', files(i).name(1:end-4)]);
                fprintf('\n');
            end
        f1 = '';
        end
    else
        if isfile(['animations/', f1, '.mat']) == 0
            fprintf(' !> File does not exist, try again.\n');
            f1 = '';
        end
    end
end

f2 = '';
while isempty(f2) == 1
    f2 = input('> Second simulation ([L]ist): ', 's');
    if f2 == 'l'
        files = dir('animations/*.mat');
        for i=1:length(files)
            fprintf(['\t', files(i).name(1:end-4)]);
            fprintf('\n');
        end
        f2 = '';
    else
        if isfile(['animations/', f2, '.mat']) == 0
            fprintf(' !> File does not exist, try again.\n');
            f2 = '';
        end
    end
end

f3 = '';
while isempty(f3) == 1
    f3 = input('> Third simulation ([N]one or [L]ist): ', 's');
    if f3 ~= 'n'
        if f3 == 'l'
            files = dir('animations/*.mat');
            for i=1:length(files)
                fprintf(['\t', files(i).name(1:end-4)]);
                fprintf('\n');
            end
            f3 = '';
        else
            if isfile(['animations/', f3, '.mat']) == 0
                fprintf(' !> File does not exist, try again.\n');
                f3 = '';
            end
        end
    end
end

if lower(f3) == 'n'
    fs = {f1; f2};
    
else
    fs = {f1; f2; f3};
end

tm = [];
xm = {};
ym = {};
Em = [];

for i=1:length(fs)
    
    f = fs{i};
    load(['animations/', f, '.mat'])
    
    if i>1 && length(t) ~= length(tm(i-1,:))
        fprintf('!> Simulations have different integration times. <!\n');
        return
    end
    tm(i,:) = t;
    xm{i} = xs(:,2:end);
    ym{i} = ys(:,2:end);
    Em(i,:) = E;
    
end

xs = zeros(length(t), n);
ys = zeros(length(t), n);

figure
set(gca,'Color', 'k')
hold on
xlabel('x'); ylabel('y')
m = sum(L) + min(L)/2;
cs = ['r', 'b', 'w'];
for i=1:length(t)
    ps = [];
    for j=1:length(fs)
        xj = xm{j}(i,:);
        yj = ym{j}(i,:);
        pj = plot([0 xj], [0 yj], 'r', 'LineWidth', 2, 'Color', cs(j));
        ps(j) = pj;
    end
    axis equal; axis([-m m -m m])
    pause(0.001)
    shg
    if i < length(t)
        for j=1:length(ps)
            set(ps(j), 'Visible', 'off')
        end
    end
end

f1 = strrep(f1, '_', '\_');
f2 = strrep(f2, '_', '\_');
if length(fs) == 2
    legend(['\color{white} ', f1], ['\color{white} ', f2])
elseif length(fs) == 3
    f3 = strrep(f3, '_', '\_');
    legend(['\color{white} ', f1], ['\color{white} ', f2], ['\color{white} ', f3])
end

figure
hold on
for i=1:length(fs)

    plot(tm(i,:), Em(i,:)) 
    
end
if length(fs) == 2
    legend(f1, f2)
elseif length(fs) == 3
    legend(f1, f2, f3)
end

if length(fs) == 2
    
    xs1 = xm{1};
    xs2 = xm{2};
    ys1 = ym{1};
    ys2 = ym{2};
    
    Ltot = sqrt(xs1(1,end)^2 + ys1(1,end)^2);
    L = Ltot/n;
    
    for i=1:n
        xe(:,i) = xs1(:,i) - xs2(:,i);
        ye(:,i) = ys1(:,i) - ys2(:,i);
    end
    
    tc = 0;
    i = 1;
    while tc < 0.001*L
        diffs = sqrt(xe(i,:).^2 + ye(i,:).^2);
        tc = max(diffs);
        i = i + 1;
    end
    i = i - 1;

    figure
    hold on
    plot(t, sqrt(xe.^2 + ye.^2))
    plot([t(i) t(i)], ylim, '--k', 'Linewidth', 2)
    nums = {};
    for j=1:n
        nums{j} = num2str(j);
    end
    legend(nums, 'Diverge')
    fprintf(['Diverge: t = ', num2str(t(i)), ' seconds.\n']);

end
