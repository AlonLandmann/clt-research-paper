% current
Fn_S1 = findopts();

% input
function [m] = new_material(name, grade, E, fb, fv, sg)
    
    m.name = name;
    m.grade = grade;
    m.E = E * 0.006895;
    m.fb = fb * 0.006895;
    m.fv = fv * 0.006895;
    m.sg = sg;
    
    m.cost = 250 + E * 0.000125 + fb * 0.0125;
    if grade == "Sel S"
        m.cost = m.cost * 1.5;
    end
    
end
function [pn] = set_panels_old(p)
    
    switch p.species
        case "DF"
            g(1) = new_material("Douglas Fir-Larch", "Sel S", 1900000, 1500, 180, 0.50);
            g(2) = new_material("Douglas Fir-Larch", "No. 1", 1700000, 1000, 180, 0.50);
            g(3) = new_material("Douglas Fir-Larch", "No. 2", 1600000,  900, 180, 0.50);
            g(4) = new_material("Douglas Fir-Larch", "No. 3", 1400000,  525, 180, 0.50);
        case "SP"
            g(1) = new_material("Southern Pine (8 in)", "Sel S", 1800000, 1950, 175, 0.55);
            g(2) = new_material("Southern Pine (8 in)", "No. 1", 1600000, 1250, 175, 0.55);
            g(3) = new_material("Southern Pine (8 in)", "No. 2", 1400000,  925, 175, 0.55);
            g(4) = new_material("Southern Pine (8 in)", "No. 3", 1300000,  525, 175, 0.55);
        case "WC"
            g(1) = new_material("Northern White Cedar", "Sel S",  800000,  775, 120, 0.31);
            g(2) = new_material("Northern White Cedar", "No. 1",  700000,  575, 120, 0.31);
            g(3) = new_material("Northern White Cedar", "No. 2",  700000,  550, 120, 0.31);
            g(4) = new_material("Northern White Cedar", "No. 3",  600000,  325, 120, 0.31);
        otherwise
            error('select a species');
    end
    
    pn(84).n = [];
    pn(84).h = [];
    pn(84).m = [];
    pn(84).o = [];
    
    i = 0;
    for h = [38.1 63.5]
        for strong = 1:4
            for weak = strong:4
                i = i + 1;
                pn(i).n = 3;
                pn(i).h = [h h h];
                pn(i).m = [g(strong) g(weak) g(strong)];
                pn(i).o = [1 0 1];
            end
        end
    end
    for h = [38.1 63.5]
        for strong = 1:4
            for weak = strong:4
                i = i + 1;
                pn(i).n = 5;
                pn(i).h = [h h h h h];
                pn(i).m = [g(strong) g(weak) g(strong) g(weak) g(strong)];
                pn(i).o = [1 0 1 0 1];
                if not(strong == weak)
                    i = i + 1;
                    pn(i).n = 5;
                    pn(i).h = [h h h h h];
                    pn(i).m = [g(strong) g(weak) g(weak) g(weak) g(strong)];
                    pn(i).o = [1 0 1 0 1];
                end
            end
        end
    end
    for h = [38.1 63.5]
        for strong = 1:4
            for weak = strong:4
                i = i + 1;
                pn(i).n = 7;
                pn(i).h = [h h h h h h h];
                pn(i).m = [g(strong) g(weak) g(strong) g(weak) g(strong) g(weak) g(strong)];
                pn(i).o = [1 0 1 0 1 0 1];
                if not(strong == weak)
                    i = i + 1;
                    pn(i).n = 7;
                    pn(i).h = [h h h h h h h];
                    pn(i).m = [g(strong) g(weak) g(weak) g(weak) g(weak) g(weak) g(strong)];
                    pn(i).o = [1 0 1 0 1 0 1];
                end
            end
        end
    end
    
end
function [pn] = set_panels_new(p)
    
    switch p.species
        case "DF"
            g(1) = new_material("Douglas Fir-Larch", "Sel S", 1900000, 1500, 180, 0.50);
            g(2) = new_material("Douglas Fir-Larch", "No. 1", 1700000, 1000, 180, 0.50);
            g(3) = new_material("Douglas Fir-Larch", "No. 2", 1600000,  900, 180, 0.50);
            g(4) = new_material("Douglas Fir-Larch", "No. 3", 1400000,  525, 180, 0.50);
        case "SP"
            g(1) = new_material("Southern Pine (8 in)", "Sel S", 1800000, 1950, 175, 0.55);
            g(2) = new_material("Southern Pine (8 in)", "No. 1", 1600000, 1250, 175, 0.55);
            g(3) = new_material("Southern Pine (8 in)", "No. 2", 1400000,  925, 175, 0.55);
            g(4) = new_material("Southern Pine (8 in)", "No. 3", 1300000,  525, 175, 0.55);
        case "WC"
            g(1) = new_material("Northern White Cedar", "Sel S",  800000,  775, 120, 0.31);
            g(2) = new_material("Northern White Cedar", "No. 1",  700000,  575, 120, 0.31);
            g(3) = new_material("Northern White Cedar", "No. 2",  700000,  550, 120, 0.31);
            g(4) = new_material("Northern White Cedar", "No. 3",  600000,  325, 120, 0.31);
        otherwise
            error('select a species');
    end
    
    pn(336).n = [];
    pn(336).h = [];
    pn(336).m = [];
    pn(336).o = [];
    
    i = 0;
    for h = [38.1 63.5]
        for strong = 1:4
            for weak = strong:4
                i = i + 1;
                pn(i).n = 3;
                pn(i).h = [h h h];
                pn(i).m = [g(strong) g(weak) g(strong)];
                pn(i).o = [1 0 1];
                i = i + 1;
                pn(i).n = 4;
                pn(i).h = [h h h h];
                pn(i).m = [g(4) g(strong) g(weak) g(strong)];
                pn(i).o = [0 1 0 1];
                i = i + 1;
                pn(i).n = 4;
                pn(i).h = [h h h h];
                pn(i).m = [g(strong) g(weak) g(strong) g(4)];
                pn(i).o = [1 0 1 0];
                i = i + 1;
                pn(i).n = 5;
                pn(i).h = [h h h h h];
                pn(i).m = [g(4) g(strong) g(weak) g(strong) g(4)];
                pn(i).o = [0 1 0 1 0];
            end
        end
    end
    for h = [38.1 63.5]
        for strong = 1:4
            for weak = strong:4
                i = i + 1;
                pn(i).n = 5;
                pn(i).h = [h h h h h];
                pn(i).m = [g(strong) g(weak) g(strong) g(weak) g(strong)];
                pn(i).o = [1 0 1 0 1];
                i = i + 1;
                pn(i).n = 6;
                pn(i).h = [h h h h h h];
                pn(i).m = [g(4) g(strong) g(weak) g(strong) g(weak) g(strong)];
                pn(i).o = [0 1 0 1 0 1];
                i = i + 1;
                pn(i).n = 6;
                pn(i).h = [h h h h h h];
                pn(i).m = [g(strong) g(weak) g(strong) g(weak) g(strong) g(4)];
                pn(i).o = [1 0 1 0 1 0];
                i = i + 1;
                pn(i).n = 7;
                pn(i).h = [h h h h h h h];
                pn(i).m = [g(4) g(strong) g(weak) g(strong) g(weak) g(strong) g(4)];
                pn(i).o = [0 1 0 1 0 1 0];
                if not(strong == weak)
                    i = i + 1;
                    pn(i).n = 5;
                    pn(i).h = [h h h h h];
                    pn(i).m = [g(strong) g(weak) g(weak) g(weak) g(strong)];
                    pn(i).o = [1 0 1 0 1];
                    i = i + 1;
                    pn(i).n = 6;
                    pn(i).h = [h h h h h h];
                    pn(i).m = [g(4) g(strong) g(weak) g(weak) g(weak) g(strong)];
                    pn(i).o = [0 1 0 1 0 1];
                    i = i + 1;
                    pn(i).n = 6;
                    pn(i).h = [h h h h h h];
                    pn(i).m = [g(strong) g(weak) g(weak) g(weak) g(strong) g(4)];
                    pn(i).o = [1 0 1 0 1 0];
                    i = i + 1;
                    pn(i).n = 7;
                    pn(i).h = [h h h h h h h];
                    pn(i).m = [g(4) g(strong) g(weak) g(weak) g(weak) g(strong) g(4)];
                    pn(i).o = [0 1 0 1 0 1 0];
                end
            end
        end
    end
    for h = [38.1 63.5]
        for strong = 1:4
            for weak = strong:4
                i = i + 1;
                pn(i).n = 7;
                pn(i).h = [h h h h h h h];
                pn(i).m = [g(strong) g(weak) g(strong) g(weak) g(strong) g(weak) g(strong)];
                pn(i).o = [1 0 1 0 1 0 1];
                i = i + 1;
                pn(i).n = 8;
                pn(i).h = [h h h h h h h h];
                pn(i).m = [g(4) g(strong) g(weak) g(strong) g(weak) g(strong) g(weak) g(strong)];
                pn(i).o = [0 1 0 1 0 1 0 1];
                i = i + 1;
                pn(i).n = 8;
                pn(i).h = [h h h h h h h h];
                pn(i).m = [g(strong) g(weak) g(strong) g(weak) g(strong) g(weak) g(strong) g(4)];
                pn(i).o = [1 0 1 0 1 0 1 0];
                i = i + 1;
                pn(i).n = 9;
                pn(i).h = [h h h h h h h h h];
                pn(i).m = [g(4) g(strong) g(weak) g(strong) g(weak) g(strong) g(weak) g(strong) g(4)];
                pn(i).o = [0 1 0 1 0 1 0 1 0];
                if not(strong == weak)
                    i = i + 1;
                    pn(i).n = 7;
                    pn(i).h = [h h h h h h h];
                    pn(i).m = [g(strong) g(weak) g(weak) g(weak) g(weak) g(weak) g(strong)];
                    pn(i).o = [1 0 1 0 1 0 1];
                    i = i + 1;
                    pn(i).n = 8;
                    pn(i).h = [h h h h h h h h];
                    pn(i).m = [g(4) g(strong) g(weak) g(weak) g(weak) g(weak) g(weak) g(strong)];
                    pn(i).o = [0 1 0 1 0 1 0 1];
                    i = i + 1;
                    pn(i).n = 8;
                    pn(i).h = [h h h h h h h h];
                    pn(i).m = [g(strong) g(weak) g(weak) g(weak) g(weak) g(weak) g(strong) g(4)];
                    pn(i).o = [1 0 1 0 1 0 1 0];
                    i = i + 1;
                    pn(i).n = 9;
                    pn(i).h = [h h h h h h h h h];
                    pn(i).m = [g(4) g(strong) g(weak) g(weak) g(weak) g(weak) g(weak) g(strong) g(4)];
                    pn(i).o = [0 1 0 1 0 1 0 1 0];
                end
            end
        end
    end
    
end
function [p] = set_parameters()
    
    p.species = "SP";
    p.set = "new";
    
    p.adcost = 3;
    p.lycost = 5;
    
    p.w = [];
    p.L = [];
    p.Ltodmin = [];
    p.Hmax = [];
    
    p.Fire = [];
    p.t = [];
    
end
function [s] = set_scenarios()
    
    s.w = 0.003;
    s.L = 4000:250:13000;
    
    s.Ltodmin = 360;
    s.Hmax = Inf;
    
    s.Fire = 0;
    s.t = 60;
    
end

% set up
function [n, h, m, o, allburnt, xcb] = after_fire_bottom(p, n, h, m, o)
    
    if p.t < 20
        xt = 7 * p.t / 20;
    else
        xt = 7;
    end
    
    if 0.65 * p.t < h(1)
        xcb = 0.65 * p.t + xt;
    else
        xcb = 0.80 * p.t + xt;
    end
    
    if xcb >= sum(h)
        nchar = n-1;
        nremain = 1;
        hremain = 0;
        allburnt = 1;
    else
        y = 0;
        for i = 1:n
            y = y + h(i);
            if xcb < y
                nchar = i - 1;
                nremain = n - nchar;
                hremain = y - xcb;
                break;
            end
        end
        allburnt = 0;
    end
    
    norig = n;
    
    n = nremain;
    h = h(nchar+1:norig);
    h(1) = hremain;
    m = m(nchar+1:norig);
    o = o(nchar+1:norig);
    
end
function [n, h, m, o, allburnt, xct] = after_fire_top(p, n, h, m, o)
    
    if p.t < 20
        xt = 7 * p.t / 20;
    else
        xt = 7;
    end
    
    if 0.65 * p.t < h(n)
        xct = 0.65 * p.t + xt;
    else
        xct = 0.80 * p.t + xt;
    end
    
    if xct >= sum(h)
        nremain = 1;
        hremain = 0;
        allburnt = 1;
    else
        y = 0;
        for i = n:-1:1
            y = y + h(i);
            if xct < y
                nremain = i;
                hremain = y - xct;
                break;
            end
        end
        allburnt = 0;
    end
    
    n = nremain;
    h = h(1:n);
    h(n) = hremain;
    m = m(1:n);
    o = o(1:n);
    
end
function [n, h, m, o, allburnt, xcb, xct] = after_fire_both(p, n, h, m, o)
    
    [n, h, m, o, allburnt, xcb] = after_fire_bottom(p, n, h, m, o);
    
    if not(allburnt)
        [n, h, m, o, allburnt, xct] = after_fire_top(p, n, h, m, o);
    else
        xct = 0;
    end
    
end
function [c] = cost(p, n, h, m)
    
    c = 0;
    for i = 1:n
        c = c + m(i).cost * (h(i) / 1000);
    end
    c = c + p.adcost * (n - 1);
    c = c + p.lycost * n;
    
end
function [wt] = load(p, n, h, m)
    
    ws = 0;
    for i = 1:n
        ws = ws + 0.00001 * m(i).sg * h(i);
    end
    wt = p.w + ws;
    
end
function [con] = constraints(p, wt, n, h, m, o)
    
    E = zeros(1, n);
    G = zeros(1, n);
    for i = 1:n
        E(i) = E(i) + m(i).E;
        G(i) = E(i) / 16;
        if o(i) == 0
            E(i) = E(i) / 30;
            G(i) = G(i) / 10;
        end
    end
    
    yb = 0;
    s1 = 0;
    s2 = 0;
    for i = 1:n
        s1 = s1 + E(i) * h(i) * (yb + h(i)/2);
        s2 = s2 + E(i) * h(i);
        yb = yb + h(i);
    end
    na = s1 / s2;
    
    y = zeros(1, n);
    z = zeros(1, n);
    I = zeros(1, n);
    yb = 0;
    for i = 1:n
        y(i) = yb + h(i) / 2;
        z(i) = y(i) - na;
        I(i) = h(i)^3 / 12;
        yb = yb + h(i);
    end
    
    BA = zeros(1, n);
    BB = zeros(1, n);
    for i = 1:n
        BA(i) = E(i) * I(i);
        BB(i) = E(i) * h(i) * z(i)^2;
    end
    EIeff = sum(BA) + sum(BB);
    
    a = 0;
    d = 0;
    for i = 1:n
        if (i == 1) || (i == n)
            a = a + h(i) / 2;
            d = d + h(i) / (2 * G(i));
        else
            a = a + h(i);
            d = d + h(i) / G(i);
        end
    end
    GAeff = a^2 / d;
    
    mcoeff = 5/384;
    scoeff = 1/8;
    gamma = (mcoeff * p.L^2 / sum(BA)) / ((mcoeff * p.L^2 / sum(BB)) + (scoeff / GAeff));
    percA = 1 / (1 + gamma);
    percB = gamma / (1 + gamma);
    
    cA = ((BA / sum(BA)) .* (h / 2) ./ I) * percA;
    if sum(BB) > 0
        cB = (E .* z) / sum(BB) * percB;
    else
        cB = 0;
    end
    cf = zeros(1, n);
    fmax = zeros(1, n);
    for i = 1:n
        r1 = abs(cB + cA);
        r2 = abs(cB - cA);
        cf(i) = max(r1(i), r2(i));
        fmax(i) = m(i).fb;
    end
    
    con.M = (1/8) * wt * p.L^2;
    con.Mmax = min(fmax ./ cf);
    con.Mval = con.M / con.Mmax;
    
    mcoeff = 999801 / 2400000000;
    scoeff = 99 / 20000;
    gamma = (mcoeff * p.L^2 / sum(BA)) / ((mcoeff * p.L^2 / sum(BB)) + (scoeff / GAeff));
    percA = 1 / (1 + gamma);
    percB = gamma / (1 + gamma);
    
    cA = (((BA / sum(BA)) * 1.5) ./ h) * percA;
    if n > 1
        cC = zeros(1, n-1);
        for i = 1:n-1
            sub = 0;
            for j = i+1:n
                sub = sub + E(j) * h(j) * z(j);
            end
            cC(i) = (percB / sum(BB)) * sub;
        end
        cB = zeros(1, n);
        cB(1) = cC(1);
        for i = 2:n-1
            cB(i) = max(cC(i-1), cC(i));
        end
        cB(n) = cC(n-1);
        ct = cA + cB;
    else
        ct = cA;
    end
    tmax = zeros(1, n);
    for i = 1:n
        tmax(i) = m(i).fv;
    end
    
    con.V = 0.49 * wt * p.L;
    con.Vmax = min(tmax ./ ct);
    con.Vval = con.V / con.Vmax;
    
    con.D = (5/384) * (wt * p.L^4 / EIeff) + (1/8) * (wt * p.L^2 / GAeff);
    con.Dmax = p.L / p.Ltodmin;
    con.Dval = con.D / con.Dmax;
    
    ld = 0;
    for i = 1:n
        ld = ld + h(i) * m(i).sg;
    end
    
    con.L = p.L;
    con.Lmax = 1000 * 0.11 * (EIeff / 10^3)^0.29 / ld^0.12;
    con.Lval = con.L / con.Lmax;
    
    con.H = sum(h(1:n));
    con.Hmax = p.Hmax;
    con.Hval = con.H / con.Hmax;
    
end
function [ch] = check(panel)
    
    ineqcon(1) = panel.con.M - panel.con.Mmax;
    ineqcon(2) = panel.con.V - panel.con.Vmax;
    ineqcon(3) = panel.con.D - panel.con.Dmax;
    ineqcon(4) = panel.con.L - panel.con.Lmax;
    ineqcon(5) = panel.con.H - panel.con.Hmax;
    
    if not(panel.AF.allburnt)
        ineqcon(6) = panel.AF.con.M - panel.AF.con.Mmax;
        ineqcon(7) = panel.AF.con.V - panel.AF.con.Vmax;
    else
        ineqcon(6) = 1;
        ineqcon(7) = 1;
    end
    
    if max(ineqcon) <= 0
        ch = 1;
    else
        ch = 0;
    end
    
end

% execution
function [r] = findopt()
    
    global p
    
    switch p.set
        case "old"
            pn = set_panels_old(p);
        case "new"
            pn = set_panels_new(p);
        otherwise
            error('set error');
    end
    costs = Inf(1, length(pn));
    
    for i = 1:length(pn)
        pn(i).cost = cost(p, pn(i).n, pn(i).h, pn(i).m);
        pn(i).wt = load(p, pn(i).n, pn(i).h, pn(i).m);
        pn(i).con = constraints(p, pn(i).wt, pn(i).n, pn(i).h, pn(i).m, pn(i).o);
        switch p.Fire
            case -1
                [pn(i).AF.n, pn(i).AF.h, pn(i).AF.m, pn(i).AF.o, pn(i).AF.allburnt, pn(i).AF.xcb] = after_fire_bottom(p, pn(i).n, pn(i).h, pn(i).m, pn(i).o);
                pn(i).AF.xct = 0;
            case 1
                [pn(i).AF.n, pn(i).AF.h, pn(i).AF.m, pn(i).AF.o, pn(i).AF.allburnt, pn(i).AF.xct] = after_fire_top(p, pn(i).n, pn(i).h, pn(i).m, pn(i).o);
                pn(i).AF.xcb = 0;
            case 0
                [pn(i).AF.n, pn(i).AF.h, pn(i).AF.m, pn(i).AF.o, pn(i).AF.allburnt, pn(i).AF.xcb, pn(i).AF.xct] = after_fire_both(p, pn(i).n, pn(i).h, pn(i).m, pn(i).o);
            otherwise
                pn(i).AF.n = 0; 
                pn(i).AF.h = 0; 
                pn(i).AF.m = 0; 
                pn(i).AF.o = 0; 
                pn(i).AF.allburnt = 1; 
                pn(i).AF.xcb = 0;
                pn(i).AF.xct = 0;
        end
        if not(pn(i).AF.allburnt)
            pn(i).AF.con = constraints(p, pn(i).wt, pn(i).AF.n, pn(i).AF.h, pn(i).AF.m, pn(i).AF.o);
        else
            pn(i).AF.con = 0;
        end
        if check(pn(i))
            costs(i) = pn(i).cost;
        end
    end
    
    index = 0;
    current = Inf;
    for i = 1:length(costs)
        if costs(i) < current
            index = i;
            current = costs(i);
        end
    end
    
    r.p = p;
    
    if index > 0
        r.s = 1;
        r.k = 1;
        r.q = 1;
        r.c = 1;
        r.opt = pn(index);
    else
        r.s = 0;
        r.k = 1;
        r.q = 0;
        r.c = 0;
        r.opt.n = 0;
        r.opt.h = 0;
        r.opt.m = 0;
        r.opt.o = 0;
        r.opt.cost = Inf;
        r.opt.wt = 0;
        r.opt.con = 0;
        r.opt.AF = 0;
    end
    
end
function [R] = findopts()
    
    s = set_scenarios();
    
    handles = ["w", "L", "Hmax", "Ltodmin", "Fire", "t"];
    d = length(handles);
    
    lengths = ones(1, d);
    for i = 1:d
        lengths(i) = length(s.(handles(i)));
    end
    
    values = zeros(d, max(lengths));
    for i = 1:d
        for j = 1:lengths(i)
            values(i, j) = s.(handles(i))(j);
        end
    end
    
    R = multirun(handles, lengths, values);
    
end
function [R] = multirun(handles, lengths, values)
    
    global p count total
    p = set_parameters();
    count = 0;
    total = prod(lengths);
    
    dim = length(handles);
    i = ones(1, dim);
    
    R(total).p = [];
    R(total).s = [];
    R(total).k = [];
    R(total).q = [];
    R(total).c = [];
    R(total).opt = [];
    
    R = dive(R, dim, handles, lengths, values, i, 1);

end
function [R, i] = dive(R, dim, handles, lengths, values, i, d)
    
    global p count total
    
    if d < dim
        for k = 1:lengths(d)
            i(d) = k;
            p.(handles(d)) = values(d, i(d));
            [R, i] = dive(R, dim, handles, lengths, values, i, d + 1);
        end
    else
        for k = 1:lengths(d)
            i(d) = k;
            p.(handles(d)) = values(d, i(d));
            count = count + 1;
            R(count) = findopt();
            disp(count + " / " + total);
        end
    end
    
end