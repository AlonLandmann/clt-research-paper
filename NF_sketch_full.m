% current
% full(NFs_DF3, [4250 6250 8750], 0.003, 'All Full Sketches', '1 DF3', 'A');
% full(NF_DF3, [4250 6250 8750], 0.003, 'All Full Sketches', '1 DF3', 'Z');
% full(NFs_SP3, [4500 6750 9000], 0.003, 'All Full Sketches', '2 SP3', 'A');
% full(NF_SP3, [4500 6750 9000], 0.003, 'All Full Sketches', '2 SP3', 'Z');
% full(NFs_WC3, [4000 5750 7250], 0.003, 'All Full Sketches', '3 WC3', 'A');
% full(NF_WC3, [4000 5750 7250], 0.003, 'All Full Sketches', '3 WC3', 'Z');
% full(NF_DF3, [7000 8750 10500], 0.003, 'All Full Sketches', '4 DF3', 'Z');
% full(NF_SP3, [7000 8750 10500], 0.003, 'All Full Sketches', '4 SP3', 'Z');
% full(NF_WC3, [7000 8750 10500], 0.003, 'All Full Sketches', '4 WC3', 'Z');
% full(NF_DF3, [5250 8000 10000], 0.003, 'All Full Sketches', '5 DF3', 'Z');
% full(NF_WC3, [5250 8000 10000], 0.003, 'All Full Sketches', '5 WC3', 'Z');
% full(NF_DFxWC3, [5250 8000 10000], 0.003, 'All Full Sketches', '5 DFxWC3', 'Z');
% full(NF_SP3, [5250 8000 10000], 0.003, 'All Full Sketches', '6 SP3', 'Z');
% full(NF_WC3, [5250 8000 10000], 0.003, 'All Full Sketches', '6 WC3', 'Z');
% full(NF_SPxWC3, [5250 8000 10000], 0.003, 'All Full Sketches', '6 SPxWC3', 'Z');
full(NF_DF3s, [4250 6250 8750], 0.003, 'All Full Sketches', '1 DF3', 'X');

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
function [nm] = set_name(p, prefix, suffix)
    
    if not(strcmp(prefix, ""))
        nm = prefix + " ";
    else
        nm = "";
    end
    nm = nm + p.L + " mm ";
    nm = nm + 1000 * p.w + " kPa";
    if not(strcmp(suffix, ""))
        nm = nm + " " + suffix;
    end
    
end
function [handles] = paint_panel(pn)
    
    handles = gobjects(1, pn.n);
    pn = set_colors(pn);
    y = 0;
    hold on
    for i = 1:pn.n
        handles(i) = fill([20 20 0 0], [y (y + pn.h(i)) (y + pn.h(i)) y] / 25.4 , pn.col(i,:));
        if pn.o(i)
            fill([-0.7 -0.7 -0.5 -0.5], [y (y + pn.h(i)) (y + pn.h(i)) y] / 25.4 , 'k');
        else
            fill([-0.7 -0.7 -0.5 -0.5], [y (y + pn.h(i)) (y + pn.h(i)) y] / 25.4 , [0.96 0.96 0.96])
        end
        y = y + pn.h(i);
    end
    axis([-2 22 -2 30]);
    axis on
    grid on
    box on
    
end
function [] = paint_info(pn)
    
    cost = "Cost: " + string(round(pn.cost)) + " $/m^2";
    depth = "Depth: " + round(sum(pn.h(1:pn.n)) / 25.4) + " in, " + round(sum(pn.h(1:pn.n)), 1) + " mm";
    Minfo = "Moment: " + round(pn.con.M / 1000) + " kNm/m | " + round(pn.con.Mmax / 1000) + " kNm/m";
    Vinfo = "Shear: " + round(pn.con.V) + " kN/m | " + round(pn.con.Vmax) + " kN/m";
    Dinfo = "Deflections: " + round(pn.con.D) + " mm | " + round(pn.con.Dmax) + " mm";
    Linfo = "Vibrations: " + round(pn.con.L) + " mm | " + round(pn.con.Lmax) + " mm";
    dim = [0.138 0.713 0.2 0.2];
    str = {"PANEL INFORMATION", cost, depth, Minfo, Vinfo, Dinfo, Linfo};
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'BackgroundColor', 'w');
    
end
function [] = paint_legend(pn, handles)
    
    orig = zeros(1, pn.n);
    desc = strings(1, pn.n);
    for i = 1:pn.n
        if not(ismember(pn.m(i).name + ", " + pn.m(i).grade, desc))
            desc(i) = pn.m(i).name + ", " + pn.m(i).grade;
            orig(i) = 1;
        end
    end
    d = strings(1, sum(orig));
    h = gobjects(1, sum(orig));
    for i = 1:pn.n
        if orig(i)
            d(sum(orig(1:i))) = desc(i);
            h(sum(orig(1:i))) = handles(i);
        end
    end
    legend(h, d, 'Fontsize', 12);
    
end
function [] = show(R, L, w)
    
    for i = 1:length(R)
        if (R(i).s > 0) && ismember(R(i).p.w, w(1)) && ismember(R(i).p.L, L(1))
            figure('visible', 'on', 'units', 'normalized', 'outerposition', [0 0 1 1]);
            handles = paint_panel(R(i).opt);
            paint_info(R(i).opt);
            paint_legend(R(i).opt, handles);
        end
    end
    
end
function [] = save(R, L, w, folders, prefix, suffix)
    
    basepath = "C:\Users\rafae\Documents\CLT\Current\Outputs";
    for i = 1:length(R)
        if (R(i).s > 0) && ismember(R(i).p.w, w) && ismember(R(i).p.L, L)
            f = figure('visible', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1]);
            handles = paint_panel(R(i).opt);
            paint_info(R(i).opt);
            paint_legend(R(i).opt, handles);
            name = set_name(R(i).p, prefix, suffix);
            saveas(f, char(basepath + "\" + folders + "\" + name + ".png"));
        end
    end
    
end
function [] = full(R, L, w, folders, prefix, suffix)
    
    switch nargin
        case 3
            show(R, L, w);
        case 4
            save(R, L, w, folders, '', '');
        case 6
            save(R, L, w, folders, prefix, suffix);
        otherwise
            error('argument error');
    end
    
end