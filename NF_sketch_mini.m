% current
mini(NFs_DF3, [4250 6250 8750], 0.003, 'Mini Sketches', '1 DF3', 'A');
mini(NF_DF3, [4250 6250 8750], 0.003, 'Mini Sketches', '1 DF3', 'Z');
mini(NFs_SP3, [4500 6750 9000], 0.003, 'Mini Sketches', '2 SP3', 'A');
mini(NF_SP3, [4500 6750 9000], 0.003, 'Mini Sketches', '2 SP3', 'Z');
mini(NFs_WC3, [4000 5750 7250], 0.003, 'Mini Sketches', '3 WC3', 'A');
mini(NF_WC3, [4000 5750 7250], 0.003, 'Mini Sketches', '3 WC3', 'Z');
mini(NF_DF3, [7000 8750 10500], 0.003, 'Mini Sketches', '4 DF3', 'Z');
mini(NF_SP3, [7000 8750 10500], 0.003, 'Mini Sketches', '4 SP3', 'Z');
mini(NF_WC3, [7000 8750 10500], 0.003, 'Mini Sketches', '4 WC3', 'Z');
mini(NF_DF3, [5250 8000 10000], 0.003, 'Mini Sketches', '5 DF3', 'Z');
mini(NF_WC3, [5250 8000 10000], 0.003, 'Mini Sketches', '5 WC3', 'Z');
mini(NF_DFxWC3, [5250 8000 10000], 0.003, 'Mini Sketches', '5 DFxWC3', 'Z');
mini(NF_SP3, [5250 8000 10000], 0.003, 'Mini Sketches', '6 SP3', 'Z');
mini(NF_WC3, [5250 8000 10000], 0.003, 'Mini Sketches', '6 WC3', 'Z');
mini(NF_SPxWC3, [5250 8000 10000], 0.003, 'Mini Sketches', '6 SPxWC3', 'Z');

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
function [] = paint_panel(pn)
    
    pn = set_colors(pn);
    y = 0;
    hold on
    for i = 1:pn.n
        fill([4 4 0 0], [y (y + pn.h(i)) (y + pn.h(i)) y] / 25.4 , pn.col(i,:));
        if pn.o(i)
            fill([-0.3 -0.3 -0.2 -0.2], [y (y + pn.h(i)) (y + pn.h(i)) y] / 25.4 , 'k');
        else
            fill([-0.3 -0.3 -0.2 -0.2], [y (y + pn.h(i)) (y + pn.h(i)) y] / 25.4 , [0.96 0.96 0.96])
        end
        y = y + pn.h(i);
    end
    axis([-1 5 -1 22]);
    axis off
    
end
function [] = set_figure(f)
    
    f.OuterPosition = [0 0 150 200];
    
end
function [] = show(R, L, w)
    
    for i = 1:length(R)
        if (R(i).s > 0) && ismember(R(i).p.w, w(1)) && ismember(R(i).p.L, L(1))
            f = figure('visible', 'on');
            paint_panel(R(i).opt);
            set_figure(f);
        end
    end
    
end
function [] = save(R, L, w, folders, prefix, suffix)
    
    basepath = "C:\Users\rafae\Documents\CLT\Current\Outputs";
    for i = 1:length(R)
        if (R(i).s > 0) && ismember(R(i).p.w, w) && ismember(R(i).p.L, L)
            f = figure('visible', 'off');
            paint_panel(R(i).opt);
            set_figure(f);
            name = set_name(R(i).p, prefix, suffix);
            saveas(f, char(basepath + "\" + folders + "\" + name + ".png"));
        end
    end
    
end
function [] = mini(R, L, w, folders, prefix, suffix)
    
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