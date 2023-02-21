% current
mini(Fs_S1, [4750 7500 8500], 0, 60, 'All Mini Sketches', '8 Sym 60', 'A');
mini(Fn_S1, [4750 7500 8500], 0, 60, 'All Mini Sketches', '8 Sym 60', 'B');
mini(F_S1, [4750 7500 8500], 0, 60, 'All Mini Sketches', '8 Sym 60', 'Z');
mini(Fs_S2, [5000 7000 9000], 0, 120, 'All Mini Sketches', '9 Sym 120', 'A');
mini(Fn_S2, [5000 7000 9000], 0, 120, 'All Mini Sketches', '9 Sym 120', 'B');
mini(F_S2, [5000 7000 9000], 0, 120, 'All Mini Sketches', '9 Sym 120', 'Z');
mini(Fs_A1, [4250 6500 8500], -1, 60, 'All Mini Sketches', '10 Asym 60', 'A');
mini(Fn_A1, [4250 6500 8500], -1, 60, 'All Mini Sketches', '10 Asym 60', 'B');
mini(F_A1, [4250 6500 8500], -1, 60, 'All Mini Sketches', '10 Asym 60', 'Z');
mini(Fs_A2, [4500 6500 8500], -1, 120, 'All Mini Sketches', '11 Asym 120', 'A');
mini(Fn_A2, [4500 6500 8500], -1, 120, 'All Mini Sketches', '11 Asym 120', 'B');
mini(F_A2, [4500 6500 8500], -1, 120, 'All Mini Sketches', '11 Asym 120', 'Z');

% function
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
    nm = nm + 1000 * p.w + " kPa ";
    switch p.Fire
        case -1
            Finfo = "below";
        case 1
            Finfo = "above";
        case 0
            Finfo = "sym";
        otherwise
            error('Fire error');
    end
    nm = nm + Finfo + " ";
    nm = nm + p.t + " min";
    if not(strcmp(suffix, ""))
        nm = nm + " " + suffix;
    end
    
end
function [] = paint_panel(pn)
    
    pn = set_colors(pn);
    pn.AF = set_colors(pn.AF);
    y = 0;
    hold on
    for i = 1:pn.n
        fill([2 2 0 0], [y (y + pn.h(i)) (y + pn.h(i)) y] / 25.4 , pn.col(i,:));
        if pn.o(i)
            fill([-0.3 -0.3 -0.2 -0.2], [y (y + pn.h(i)) (y + pn.h(i)) y] / 25.4 , 'k');
        else
            fill([-0.3 -0.3 -0.2 -0.2], [y (y + pn.h(i)) (y + pn.h(i)) y] / 25.4 , [0.96 0.96 0.96])
        end
        y = y + pn.h(i);
    end
    y = 0;
    if pn.AF.xcb > 0
        fill([4 4 2 2], [y (y + pn.AF.xcb) (y + pn.AF.xcb) y] / 25.4 , 'k');
        y = pn.AF.xcb;
    end
    for i = 1:pn.AF.n
        fill([4 4 2 2], [y (y + pn.AF.h(i)) (y + pn.AF.h(i)) y] / 25.4 , pn.AF.col(i,:));
        y = y + pn.AF.h(i);
    end
    if pn.AF.xct > 0
        fill([4 4 2 2], [y (y + pn.AF.xct) (y + pn.AF.xct) y] / 25.4 , 'k');
    end
    axis([-1 5 -1 28]);
    axis off
    
end
function [] = set_figure(f)
    
    f.OuterPosition = [0 0 150 200];
    
end
function [] = show(R, L, Fire, t)
    
    for i = 1:length(R)
        if (R(i).s > 0) && ismember(R(i).p.L, L(1)) && ismember(R(i).p.Fire, Fire(1)) && ismember(R(i).p.t, t(1))
            f = figure('visible', 'on');
            paint_panel(R(i).opt);
            set_figure(f);
        end
    end
    
end
function [] = save(R, L, Fire, t, folders, prefix, suffix)
    
    basepath = "C:\Users\Alon\Documents\Saved Documents\Career\Mass Timber\CLT optimization\Current\Outputs";
    for i = 1:length(R)
        if (R(i).s > 0) && ismember(R(i).p.L, L) && ismember(R(i).p.Fire, Fire) && ismember(R(i).p.t, t)
            f = figure('visible', 'off');
            paint_panel(R(i).opt);
            set_figure(f);
            name = set_name(R(i).p, prefix, suffix);
            saveas(f, char(basepath + "\" + folders + "\" + name + ".png"));
        end
    end
    
end
function [] = mini(R, L, Fire, t, folders, prefix, suffix)
    
    switch nargin
        case 4
            show(R, L, Fire, t);
        case 5
            save(R, L, Fire, t, folders, '', '');
        case 7
            save(R, L, Fire, t, folders, prefix, suffix);
        otherwise
            error('argument error');
    end
    
end