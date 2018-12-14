%-------------------------------------------------------------------------%
% antovid - Jeremy Turner
% 
% Saves animation as .avi
%
% Input: fn - filename, e.g. 'animations/DAEp100.mat'
%      : svnm (optional) - filename to save to, fn if unspecified
%
% Loads in all the fields of the saved simulation. Animation is the same as
% in npend.m, see there for details. Each frame saved to .avi file.
% ------------------------------------------------------------------------%

function antovid = func(fn, svnm)

    load(fn)
    
    if ~exist('svnm', 'var')
        svnm = [fn(1:end-3), 'avi'];
    end
    
    v = VideoWriter(svnm);
    open(v);
    
    Lnetx = sum(L) + min(L)/2 - abs(xend);
    Lnety = sum(L) + min(L)/2 - abs(yend);
    xmin = -Lnetx - (xend - abs(xend))/2;
    xmax = Lnetx + (xend + abs(xend))/2;
    ymin = -Lnety - (yend - abs(yend))/2;
    ymax = Lnety + (yend + abs(yend))/2;
       
    m = sum(L) + min(L)/2;

    fig = figure('units', 'normalized', 'outerposition', [0.1 0.1 0.9 0.9]);
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
        pause(0.01);
        frame = getframe(gcf);
        writeVideo(v, frame);
        if i < length(t)
            set(pl, 'Visible', 'off')
        end            
    end
    
    close(v);
    
    close all

end