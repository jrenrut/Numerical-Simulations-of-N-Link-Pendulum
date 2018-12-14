%-------------------------------------------------------------------------%
% teval - Jeremy Turner
% 
% This function runs the symbolic deriver methods to analyze their
% performance. Given a number of runs to average over (m) and the numbers
% of links to use (n), teval.m will run the specified method for n links m
% times each. The data will be saved to a file only specifying the method
% used, so that if teval.m is run again with the same method, the data will
% be over-written unless the original file has been renamed. Additionally,
% teval will plot the performance data.
% ------------------------------------------------------------------------%

clear all

ptype = '';
while isempty(ptype) == 1
    ptype = lower(input('> [N]ewton-Euler, [L]agrange, or [D]AE? ', 's'));
end
ptype = ptype(1);

prob = '';
if ptype == 'n'
    f = 'npend_deriver_NE';
    fs = 'Newton-Euler Pendulum';
elseif ptype == 'l'
    f = 'npend_deriver_L';
    fs = 'Lagrange Pendulum';
elseif ptype == 'd'
    while isempty(prob) == 1
        prob = lower(input('> [P]endulum or [L]inkage? ', 's'));
    end
    if prob == 'p'
        f = 'npend_deriver_DAE';
        fs = 'DAE Pendulum';
    else
        f = 'nlink_deriver_DAR';
        fs = 'DAE Linkage';
    end
end

mq = '';
while isempty(mq) == 1
    mq = input('> Iterations: ', 's');
end
m = str2num(mq);

nq = '';
while isempty(nq) == 1
    nq = input('> Numbers of Links (separated by spaces): ', 's');
end

nq = strsplit(nq, ' ');
ns = zeros(1, length(nq));

for i=1:length(nq)
    ns(i) = str2num(nq{i});
end

mt = zeros(m,length(ns));

for i=1:m

    tn = zeros(1, length(ns));
    
    feval(f, 1)

    for j=1:length(ns)
        if fs(1:5) == 'DAE L'
            n = ns(j) + 1;
        else
            n = ns(j);
        end
        pause(0.01)
        tic
        feval(f, n)
        t = toc;
        tn(j) = t;
        fprintf(['Round ', num2str(i), ': n = ', num2str(n), '; t = ', num2str(t), '\n']);
    end

    mt(i,:) = tn;
    
end

tsd = zeros(1, length(ns));
tav = zeros(1, length(ns));

for i=1:length(ns)

    ts = mt(:,i);
    tsd(i) = std(ts);
    tav(i) = mean(ts);

end

save([fs, '_teval.mat'], 'ns', 'tav', 'tsd', 'm', 'fs')

errorbar(ns, tav, tsd); xlabel('n'); ylabel('t [s]');
xlim([0 ns(end)+(ns(end)-ns(end-1))/2]);
title([fs, ' - Computation Time - ', num2str(m), ' Iterations']);