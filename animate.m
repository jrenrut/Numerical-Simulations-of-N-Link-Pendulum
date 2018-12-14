%-------------------------------------------------------------------------%
% animate - Jeremy Turner
% 
% Animates a saved simulation.
%
% Input: fn - filename, e.g. 'animations/DAEp100.mat'
%
% Loads in all the fields of the saved simulation. Animation is the same as
% in npend.m, see there for details.
% ------------------------------------------------------------------------%

function animate = func(fn)
    
    load(fn)    
    
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
    plot([0 xend], [0 yend],  'color', [0.6350, 0.0780, 0.1840], 'LineWidth', 2);
    for i=1:length(t)
        pl = plot(xs(i,:), ys(i,:), 'r', 'LineWidth', 2);
        if rq ~= 0
            if i > 1
                plot([xs(i-1,2) xs(i,2)], [ys(i-1,2) ys(i,2)], ':w')
            end        
        else
            if i > 1
                plot([xs(i-1,n+1) xs(i,n+1)], [ys(i-1,n+1) ys(i,n+1)], ':w')
            end
        end
        axis equal;
        if ptype == 'd' && prob == 'l'
            axis([xmin xmax ymin ymax])
        elseif exist('dst') ~= 0 &&lower(dst(1)) == 'm'
            axis([-m/2 m/2 -1.1*sum(L) max(L)]);
        else
            axis([-m m -m m])
        end
        pause(0.005)
        shg
        if i < length(t)
            set(pl, 'Visible', 'off')
        end            
    end
    
    figure
    plot(t, E)
    xlabel('t'); ylabel('E'); title('Total Energy');
    
    if ptype == 'm'    
        mvar = 'mode';
        md = load(fn, mvar);
        md = md.(mvar);
        figure
        hold on
        xlabel('t'); ylabel('\theta'); title(['Angles vs. Time: ', num2str(n), ' Links, Mode ', num2str(md), ', \omega = ', num2str(freq), ' rad/s']);
        legs = {};
        for i=1:n
            plot(t, thetas(:,i))
            legs{i} = ['\theta_', num2str(i)];
        end
        legend(legs)
            
        figure
        hold on
        xlabel('t'); ylabel('\theta'); title({['Linearized - Actual Angles vs. Time: ', num2str(n), ' Links, Normal Mode ', num2str(md), ', \omega = ', num2str(freq), ' rad/s'], ''});
        legs = {};
        for i=1:n
            plot(t, thetas(:,i)-thetas(1,i)*cos(freq*t))
            legs{i} = ['\Delta \theta_', num2str(i)];
        end
        legend(legs, 'Location', 'northwest')
    elseif exist('dst') ~= 0 && lower(dst(1)) == 'm'
        mvar = 'mode';
        md = load(fn, mvar);
        md = md.(mvar);
        figure('units', 'normalized', 'outerposition', [0.2 0.4 0.8 0.6])
        hold on
        if exist('freq') ~= 0
            xlabel('t'); ylabel('\theta'); title(['Angles vs. Time: ', num2str(n), ' Links, Standing Mode ', num2str(md), ', \omega = ', num2str(freq), ' rad/s']);
        else
            xlabel('t'); ylabel('\theta'); 
        end
        plot(t, thetas(:,n))
        plot(t, cos(freq*t))
        ylim([-1.2 1.4])
        legend('End Frequency', 'Driving Frequency');
    end

end
