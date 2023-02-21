% current
% pictogram(NF_DF3, [4250 6250 8750], 'All pictograms\1 DF Z');
% pictogram(NFs_DF3, [4250 6250 8750], 'All pictograms\1 DF A');
% pictogram(NF_DF3s, [4250 6250 8750], 'All pictograms\1 DF S');
pictogram(NF_DF3s, [], 'All pictograms\12 DF symmetric');
% 
% pictogram(NF_SP3, [4500 6750 9000], 'All pictograms\2 SP Z');
% pictogram(NFs_SP3, [4500 6750 9000], 'All pictograms\2 SP A');
% 
% pictogram(NF_WC3, [4000 5750 7250], 'All pictograms\3 WC Z');
% pictogram(NFs_WC3, [4000 5750 7250], 'All pictograms\3 WC A');
% pictogram(NF_WC3s, [4000 5750 7250], 'All pictograms\3 WC S');
% 
% % pictogram(NF_DF3, [7000 8750 10500], 'All pictograms\4 DF Z');
% % pictogram(NF_SP3, [7000 8750 10500], 'All pictograms\4 SP Z');
% % pictogram(NF_WC3, [7000 8750 10500], 'All pictograms\4 WC Z');
% 
% pictogram(NF_DF3, [5250 8000 10000], 'All pictograms\5 DF Z');
% pictogram(NF_WC3, [5250 8000 10000], 'All pictograms\5 WC Z');
% pictogram(NF_DFxWC3, [5250 8000 10000], 'All pictograms\5 DFxWC Z');
% 
% pictogram(NF_SP3, [5250 8000 10000], 'All pictograms\6 SP Z');
% pictogram(NF_WC3, [5250 8000 10000], 'All pictograms\6 WC Z');
% pictogram(NF_SPxWC3, [5250 8000 10000], 'All pictograms\6 SPxWC Z');

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
function [ab] = abspecies(spc)
    
    switch spc
        case "Douglas Fir-Larch"
            ab = "DF";
        case "Southern Pine (8 in)"
            ab = "SP";
        case "Northern White Cedar"
            ab = "WC";
        otherwise
            error('species error');
    end
    
end
function [pn] = set_colors(pn)
    
    pn.col = zeros(pn.n, 3);
    for i = 1:pn.n
        switch pn.m(i).name
            case "Douglas Fir-Larch"
                switch pn.m(i).grade
                    case "Sel S"
                        pn.col(i,:) = colorab('gs');
                    case "No. 1"
                        pn.col(i,:) = colorab('g1');
                    case "No. 2"
                        pn.col(i,:) = colorab('g2');
                    case "No. 3"
                        pn.col(i,:) = colorab('g3');
                end
            case "Southern Pine (8 in)"
                switch pn.m(i).grade
                    case "Sel S"
                        pn.col(i,:) = colorab('rs');
                    case "No. 1"
                        pn.col(i,:) = colorab('r1');
                    case "No. 2"
                        pn.col(i,:) = colorab('r2');
                    case "No. 3"
                        pn.col(i,:) = colorab('r3');
                end
            case "Northern White Cedar"
                switch pn.m(i).grade
                    case "Sel S"
                        pn.col(i,:) = colorab('ys');
                    case "No. 1"
                        pn.col(i,:) = colorab('y1');
                    case "No. 2"
                        pn.col(i,:) = colorab('y2');
                    case "No. 3"
                        pn.col(i,:) = colorab('y3');
                end
            otherwise
                error('species error');
        end
    end
    
end
function [handles, maxL] = paint_panels(R, hl)
    
    handles = gobjects(length(R), 7);
    maxL = 4;
    hold on
    for k = 1:length(R)
        if R(k).s > 0
            pn = set_colors(R(k).opt);
            y = 0;
            for i = 1:pn.n
                Xcoor = [125 125 -125 -125] + R(k).p.L;
                Ycoor = [0 pn.h(i) pn.h(i) 0] + y;
                if pn.o(i)
                    handles(k, i) = fill(Xcoor / 1000, Ycoor , pn.col(i,:));
                else
                    handles(k, i) = fill(Xcoor / 1000, Ycoor , pn.col(i,:), 'FaceAlpha', 0.5);
                end
                y = y + pn.h(i);
            end
            if R(k).p.L/1000 > maxL
                maxL = R(k).p.L / 1000;
            end
        end
    end
    for k = 1:length(R)
        if (R(k).s > 0) && ismember(R(k).p.L, hl)
            Xcoor = [125 125 -125 -125] + R(k).p.L + [75 75 -75 -75];
            H = sum(R(k).opt.h);
            Ycoor = [0 H H 0] + [-10 10 10 -10];
            fill(Xcoor / 1000, Ycoor, 'b', 'FaceAlpha', 0, 'EdgeColor', 'b', 'LineWidth', 1.2);
        end
    end
    
end
function [] = paint_legend(R, handles)
    
    orig = zeros(length(R), 7);
    desc = strings(length(R), 7);
    for k = 1:length(R)
        if R(k).s > 0
            for i = 1:R(k).opt.n
                if R(k).opt.o(i)
                    or = "longitudonal";
                else
                    or = "transverse";
                end
                candidate = abspecies(R(k).opt.m(i).name) + ", " + R(k).opt.m(i).grade + ", " + or;
                if not(ismember(candidate, desc))
                    desc(k, i) = candidate;
                    orig(k, i) = 1;
                end
            end
        end
    end
    d = strings(1, sum(orig, 'all'));
    h = gobjects(1, sum(orig, 'all'));
    b = 0;
    for k = 1:length(R)
        if R(k).s > 0
            for i = 1:R(k).opt.n
                if orig(k, i)
                    b = b + 1;
                    d(b) = desc(k, i);
                    h(b) = handles(k, i);
                end
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
function [] = set_axis(f, maxL)
    
    axis([3 (maxL+1) (-100) 550]);
    grid on
    box on
    xlabel('Span (m)', 'FontSize', 12);
    ylabel('Thickness (mm)', 'FontSize', 12);
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.XAxis.TickValues = 0:1:20;
    ax.XAxis.LineWidth = 1;
    ax.YAxis.FontSize = 12;
    ax.YAxis.LineWidth = 1;
    ax.YAxis.TickValues = -500:50:500;
    ax.YAxis.Label.Units = 'pixels';
    ax.YAxis.Label.Position = [-42 f.OuterPosition(4)/2 - 40];
    ax.TickLength = [0 0];
    
end
function [] = set_figure(f)
    
    f.OuterPosition = [0 0 800 500];
    f.InnerPosition = f.OuterPosition - 80;
    
end
function [] = show(R, hl)
    
    f = figure('visible', 'on');
    [handles, maxL] = paint_panels(R, hl);
    paint_legend(R, handles);
    set_axis(f, maxL);
    set_figure(f);
    
end
function [] = save(R, hl, extension)
    
    basepath = "C:\Users\rafae\Documents\Career\CLT research\Current\Outputs";
    f = figure('visible', 'off');
    [handles, maxL] = paint_panels(R, hl);
    paint_legend(R, handles);
    set_axis(f, maxL);
    set_figure(f);
    saveas(f, char(basepath + "\" + extension + ".png"));
    
end
function [] = pictogram(R, hl, extension)
    
    switch nargin
        case 2
            show(R, hl);
        case 3
            save(R, hl, extension);
        otherwise
            error('argument error');
    end
    
end