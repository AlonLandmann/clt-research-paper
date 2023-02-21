% current
pictogram(F_S1, [4750 7500 8500], 1, 'All pictograms\8 Sym 60 Z charred');
pictogram(F_S2, [5000 7000 9000], 1, 'All pictograms\9 Sym 120 Z charred');
pictogram(F_A1, [4250 6500 8500], 1, 'All pictograms\10 Asym 60 Z charred');
pictogram(F_A2, [4500 6500 8500], 1, 'All pictograms\11 Asym 120 Z charred');

pictogram(Fn_S1, [4750 7500 8500], 1, 'All pictograms\8 Sym 60 B charred');
pictogram(Fn_S2, [5000 7000 9000], 1, 'All pictograms\9 Sym 120 B charred');
pictogram(Fn_A1, [4250 6500 8500], 1, 'All pictograms\10 Asym 60 B charred');
pictogram(Fn_A2, [4500 6500 8500], 1, 'All pictograms\11 Asym 120 B charred');

pictogram(Fs_S1, [4750 7500 8500], 1, 'All pictograms\8 Sym 60 A charred');
pictogram(Fs_S2, [5000 7000 9000], 1, 'All pictograms\9 Sym 120 A charred');
pictogram(Fs_A1, [4250 6500 8500], 1, 'All pictograms\10 Asym 60 A charred');
pictogram(Fs_A2, [4500 6500 8500], 1, 'All pictograms\11 Asym 120 A charred');

pictogram(F_S1, [4750 7500 8500], 0, 'All pictograms\8 Sym 60 Z');
pictogram(F_S2, [5000 7000 9000], 0, 'All pictograms\9 Sym 120 Z');
pictogram(F_A1, [4250 6500 8500], 0, 'All pictograms\10 Asym 60 Z');
pictogram(F_A2, [4250 6500 8500], 0, 'All pictograms\11 Asym 120 Z');

pictogram(Fn_S1, [4750 7500 8500], 0, 'All pictograms\8 Sym 60 B');
pictogram(Fn_S2, [5000 7000 9000], 0, 'All pictograms\9 Sym 120 B');
pictogram(Fn_A1, [4250 6500 8500], 0, 'All pictograms\10 Asym 60 B');
pictogram(Fn_A2, [4500 6500 8500], 0, 'All pictograms\11 Asym 120 B');

pictogram(Fs_S1, [4750 7500 8500], 0, 'All pictograms\8 Sym 60 A');
pictogram(Fs_S2, [5000 7000 9000], 0, 'All pictograms\9 Sym 120 A');
pictogram(Fs_A1, [4250 6500 8500], 0, 'All pictograms\10 Asym 60 A');
pictogram(Fs_A2, [4500 6500 8500], 0, 'All pictograms\11 Asym 120 A');

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
function [handles, maxL] = paint_panels(R, hl, charred)
    
    handles = gobjects(length(R), 9);
    maxL = 4;
    hold on
    for k = 1:length(R)
        if R(k).s > 0
            pn = set_colors(R(k).opt);
            y = 0;
            for i = 1:pn.n
                Xcoor = ([125 125 -125 -125] + R(k).p.L);
                Ycoor = [0 pn.h(i) pn.h(i) 0] + y;
                if pn.o(i)
                    handles(k, i) = fill(Xcoor / 1000, Ycoor / 25.4 , pn.col(i,:));
                else
                    handles(k, i) = fill(Xcoor / 1000, Ycoor / 25.4 , pn.col(i,:), 'FaceAlpha', 0.5);
                end
                y = y + pn.h(i);
            end
            if charred
                if pn.AF.xcb > 0
                    Ycoor = [0 pn.AF.xcb pn.AF.xcb 0];
                    fill(Xcoor / 1000, Ycoor / 25.4 , 'k', 'FaceAlpha', 0.8);
                end
                if pn.AF.xct > 0
                    Ycoor = [0 (-pn.AF.xct) (-pn.AF.xct) 0] + y;
                    fill(Xcoor / 1000, Ycoor / 25.4 , 'k', 'FaceAlpha', 0.8);
                end
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
            fill(Xcoor / 1000, Ycoor / 25.4, 'b', 'FaceAlpha', 0, 'EdgeColor', 'b', 'LineWidth', 1.2);
        end
    end
    
end
function [] = paint_legend(R, handles)
    
    orig = zeros(length(R), 9);
    desc = strings(length(R), 9);
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
    lgd = legend(h, d);
    lgd.Location = 'northwest';
    lgd.LineWidth = 0.6;
    lgd.FontSize = 10;
    
end
function [] = set_axis(f, maxL)
    
    axis([3 (maxL+1) (-4) 31]);
    grid on
    box on
    xlabel('Span (m)', 'FontSize', 12);
    ylabel('Depth (in)', 'FontSize', 12);
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.XAxis.TickValues = 0:1:20;
    ax.XAxis.LineWidth = 1;
    ax.YAxis.FontSize = 12;
    ax.YAxis.LineWidth = 1;
    ax.YAxis.TickValues = -50:2:50;
    ax.YAxis.Label.Units = 'pixels';
    ax.YAxis.Label.Position = [-42 f.OuterPosition(4)/2 - 40];
    ax.TickLength = [0 0];
    
end
function [] = set_figure(f)
    
    f.OuterPosition = [0 0 800 500];
    f.InnerPosition = f.OuterPosition - 80;
    
end
function [] = show(R, hl, charred)
    
    f = figure('visible', 'on');
    [handles, maxL] = paint_panels(R, hl, charred);
    paint_legend(R, handles);
    set_axis(f, maxL);
    set_figure(f);
    
end
function [] = save(R, hl, charred, extension)
    
    basepath = "C:\Users\Alon\Documents\Saved Documents\Career\Mass Timber\CLT optimization\Current\Outputs";
    f = figure('visible', 'off');
    [handles, maxL] = paint_panels(R, hl, charred);
    paint_legend(R, handles);
    set_axis(f, maxL);
    set_figure(f);
    saveas(f, char(basepath + "\" + extension + ".png"));
    
end
function [] = pictogram(R, hl, charred, extension)
    
    switch nargin
        case 3
            show(R, hl, charred);
        case 4
            save(R, hl, charred, extension);
        otherwise
            error('argument error');
    end
    
end