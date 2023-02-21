% current
compare([F_S1;Fn_S1;Fs_S1], [4750 7500 8500], 'All Scatters\8 Sym 60');
compare([F_S2;Fn_S2;Fs_S2], [5000 7000 9000], 'All Scatters\9 Sym 120');
compare([F_A1;Fn_A1;Fs_A1], [4250 6500 8500], 'All Scatters\10 Asym 60');
compare([F_A2;Fn_A2;Fs_A2], [4500 6500 8500], 'All Scatters\11 Asym 120');

% functions
function [col] = colorab(ab)
    
    abbreviations = [
        "gs"
        "g1"
        "g2"
        "g3"
        "rs"
        "r1"
        "r2"
        "r3"
        "ys"
        "y1"
        "y2"
        "y3"
    ];
    colors = [
        0.2157 0.3373 0.1373;
        0.3294 0.5098 0.2078;
    	0.6627 0.8157 0.5569;
        0.7765 0.8784 0.7059;
    	0.5137 0.2353 0.0471;
    	0.7765 0.3490 0.0667;
    	0.9569 0.6902 0.5176;
        0.9725 0.7961 0.6784;
    	0.5020 0.3765 0.0000;
    	0.7490 0.5608 0.0000;
    	1.0000 0.8510 0.4000;
        1.0000 0.9020 0.6000
    ];
    if ismember(ab, abbreviations)
        col = [1 1 1];
        for i = 1:length(abbreviations)
            if strcmp(ab, abbreviations(i))
                col = colors(i, :);
            end
        end
    else
        error('color error');
    end
    
end
function [mk, sz] = set_shape(o)
    
    switch sum(o)
        case 2
            mk = 'o';
            sz = 32;
        case 3
            mk = 's';
            sz = 40;
        case 4
            mk = '^';
            sz = 32;
        otherwise
            error('marker error');
    end
    
end
function [ecol, fcol] = set_colors(r)
    
    if isfield(r.p, 'mrange')
        if length(r.p.mrange) == 4
            switch r.p.mrange(1).name
                case "Douglas Fir-Larch"
                    ecol = colorab('g1');
                    fcol = colorab('g2');
                case "Southern Pine (8 in)"
                    ecol = colorab('r1');
                    fcol = colorab('r2');
                case "Northern White Cedar"
                    ecol = colorab('y1');
                    fcol = colorab('y2');
            end
        else
            switch r.p.mrange(1).name
                case "Douglas Fir-Larch"
                    ecol = colorab('g1');
                case "Southern Pine (8 in)"
                    ecol = colorab('r1');
                case "Northern White Cedar"
                    ecol = colorab('y1');
            end
            switch r.p.mrange(5).name
                case "Douglas Fir-Larch"
                    fcol = colorab('g2');
                case "Southern Pine (8 in)"
                    fcol = colorab('r2');
                case "Northern White Cedar"
                    fcol = colorab('y2');
            end
        end
    else
        if strcmp(r.p.set, 'old')
            switch r.opt.m(1).name
                case "Douglas Fir-Larch"
                    ecol = colorab('gs');
                    fcol = colorab('gs');
                case "Southern Pine (8 in)"
                    ecol = colorab('rs');
                    fcol = colorab('rs');
                case "Northern White Cedar"
                    ecol = colorab('ys');
                    fcol = colorab('ys');
            end
        else
            switch r.opt.m(1).name
                case "Douglas Fir-Larch"
                    ecol = colorab('g1');
                    fcol = colorab('g1');
                case "Southern Pine (8 in)"
                    ecol = colorab('r1');
                    fcol = colorab('r1');
                case "Northern White Cedar"
                    ecol = colorab('y1');
                    fcol = colorab('y1');
            end
        end
    end
        
end

function [candidate] = set_legend(r)
    
    if r.s > 0
        if isfield(r.p, 'mrange')
            candidate = "optimized";
        else
            if strcmp(r.p.set, 'old')
                candidate = "old set";
            else
                candidate = "fire set";
            end
        end
        candidate = candidate + " - " + sum(r.opt.o) + " strong layers";
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
function [candidates, handles, L, C] = paint_scatter(R, hl)
    
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
                [marker, sz] = set_shape(R(k, i).opt.o);
                [ecol, fcol] = set_colors(R(k, i));
                handles(k, i) = scatter(L(k, i), C(k, i), sz, marker);
                handles(k, i).LineWidth = 0.9;
                handles(k, i).MarkerEdgeColor = ecol;
                handles(k, i).MarkerFaceColor = fcol;
                if ismember(R(k, i).p.L, hl)
                    sc = scatter(L(k, i), C(k, i), 110, 'b', 'o', 'LineWidth', 0.9);
                    sc.MarkerEdgeColor = 'b';
                    sc.MarkerFaceColor = 'b';
                    sc.MarkerFaceAlpha = 0;
                end
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
function [] = compare(R, hl, extension)
    
    switch nargin
        case 2
            f = figure('visible', 'on');
        case 3
            f = figure('visible', 'off');
        otherwise
            error('argument error');
    end
    [candidates, handles, L, C] = paint_scatter(R, hl);
    paint_legend(handles, candidates);
    set_figure(f);
    set_axis(L, C);
    if nargin == 3
        basepath = "C:\Users\rafae\Documents\CLT\Current\Outputs";
        saveas(f, char(basepath + "\" + extension + ".png"));
    end
    
end