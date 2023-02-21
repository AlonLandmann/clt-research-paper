% current
overview([NF_DF3 NF_DF5 NF_DF7], 'All Scatters\7 Constraints');

% functions
function [col] = colbycon(con)
    
    index = 0;
    current = 0;
    pvals = [
        con.Mval;
        con.Vval;
        con.Dval;
        con.Lval
    ];
    colors = [
        0.7765 0.3490 0.0667;
        0.7490 0.5608 0.0000;
        0.3294 0.5098 0.2078;
        0.0000 0.0000 0.0000
    ];
    for i = 1:4
        if pvals(i) >= current
            index = i;
            current = pvals(i);
        end
    end
    col = colors(index,:);
    
end
function [con] = conbycol(col)
    
    cons = [
        "Moments"
        "Shear Forces"
        "Deflections"
        "Vibrations"
    ];
    colors = [
        0.7765 0.3490 0.0667;
        0.7490 0.5608 0.0000;
        0.3294 0.5098 0.2078;
        0.0000 0.0000 0.0000
    ];
    for i = 1:4
        if sum(colors(i,:) == col) == 3
            con = cons(i);
        end
    end
    
end
function [mk] = set_marker(n)
    
    switch n
        case 3
            mk = 'o';
        case 5
            mk = 's';
        case 7
            mk = '^';
        otherwise
            error('marker error');
    end
    
end
function [sz] = set_size(w)
    
    sz = (1000 * w)^3 / 3;           
    
end
function [candidate] = set_legend(r)
    
    if r.s > 0
        candidate = conbycol(colbycon(r.opt.con));
        candidate = candidate + " (" + r.opt.n + "-ply)";
    else
        candidate = "";
    end
    
end
function [] = set_axis(L, C)
    
    axis([min(nonzeros(L))-1 max(L, [], 'all')+1 min(nonzeros(C))-30 max(C, [], 'all')+30]);
    grid on
    box on
    xlabel('Span (m)', 'FontSize', 12);
    ylabel('Cost ($/m^2)', 'FontSize', 12);
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.XAxis.TickValues = 0:1:20;
    ax.XAxis.LineWidth = 1;
    ax.YAxis.FontSize = 12;
    ax.YAxis.LineWidth = 1;
    ax.YAxis.TickValues = 0:50:1000;
    ax.YAxis.Label.Units = 'pixels';
    ax.YAxis.Label.Position = [-42 240];
    ax.TickLength = [0 0];
    
end
function [] = set_figure(f)
    
    f.OuterPosition = [0 0 640 640];
    f.InnerPosition = [80 80 560 560];
    
end
function [candidates, handles, L, C] = paint_scatter(R)
    
    Rn = size(R, 1);
    Rlen = size(R, 2);
    L = zeros(Rn, Rlen);
    C = zeros(Rn, Rlen);
    candidates = strings(Rn, Rlen);
    handles = gobjects(Rn, Rlen);
    hold on
    for k = 1:Rn
        for i = 1:Rlen
            if R(k, i).s > 0
                L(k, i) = R(k, i).p.L / 1000;
                C(k, i) = R(k, i).opt.cost;
                candidates(k, i) = set_legend(R(k, i));
                marker = set_marker(R(k, i).opt.n);
                sz = set_size(R(k, i).p.w);
                color = colbycon(R(k, i).opt.con);
                handles(k, i) = scatter(L(k, i), C(k, i), sz, color, marker);
                handles(k, i).LineWidth = 0.65;
            end
        end
    end
    
end
function [] = paint_legend(handles, candidates)
    
    rn = size(handles, 1);
    rlen = size(handles, 2);
    orig = zeros(rn, rlen);
    desc = strings(rn, rlen);
    for k = 1:rn
        for i = 1:rlen
            if not(strcmp(candidates(k, i), "")) && not(ismember(candidates(k, i), desc))
                desc(k, i) = candidates(k, i);
                orig(k, i) = 1;
            end
        end
    end
    d = strings(1, sum(orig, 'all'));
    h = gobjects(1, sum(orig, 'all'));
    b = 0;
    for k = 1:rn
        for i = 1:rlen
            if orig(k, i)
                b = b + 1;
                d(b) = desc(k, i);
                h(b) = handles(k, i);
            end
        end
    end
    % sorting alphabetically
    [ds, change] = sort(d);
    hs = gobjects(1, length(d));
    for g = 1:length(d)
        hs(g) = h(change(g));
    end
    lgd = legend(hs, ds);
    lgd.Location = 'northwest';
    lgd.LineWidth = 0.6;
    lgd.FontSize = 10;
    
end
function [] = overview(R, extension)
    
    switch nargin
        case 1
            f = figure('visible', 'on');
        case 2
            f = figure('visible', 'off');
        otherwise
            error('argument error');
    end
    [candidates, handles, L, C] = paint_scatter(R);
    paint_legend(handles, candidates);
    set_figure(f);
    set_axis(L, C);
    if nargin == 2
        basepath = "C:\Users\rafae\Documents\CLT\Current\Outputs";
        saveas(f, char(basepath + "\" + extension + ".png"));
    end
    
end