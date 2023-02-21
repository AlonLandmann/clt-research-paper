% current
F_DF3_S60 = findopts();

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
function [mrange, unav] = set_available()
    
    mrange = [
        new_material("Douglas Fir-Larch", "Sel S", 1900000, 1500, 180, 0.50)
        new_material("Douglas Fir-Larch", "No. 1", 1700000, 1000, 180, 0.50)
        new_material("Douglas Fir-Larch", "No. 2", 1600000,  900, 180, 0.50)
        new_material("Douglas Fir-Larch", "No. 3", 1400000,  525, 180, 0.50)
    ];
    unav = [
        new_material("Southern Pine (8 in)", "Sel S", 1800000, 1950, 175, 0.55)
        new_material("Southern Pine (8 in)", "No. 1", 1600000, 1250, 175, 0.55)
        new_material("Southern Pine (8 in)", "No. 2", 1400000,  925, 175, 0.55)
        new_material("Southern Pine (8 in)", "No. 3", 1300000,  525, 175, 0.55)
        new_material("Northern White Cedar", "Sel S",  800000,  775, 120, 0.31)
        new_material("Northern White Cedar", "No. 1",  700000,  575, 120, 0.31)
        new_material("Northern White Cedar", "No. 2",  700000,  550, 120, 0.31)
        new_material("Northern White Cedar", "No. 3",  600000,  325, 120, 0.31)
    ];
    
end
function [p] = set_parameters()
    
    p.tmax = 300;
    p.kmax = 100;
    p.smax = [10 20 30];
    p.cmin = [0.8 0.5 0.2];
    
    p.nrange = [3 4 5 6 7 8 9];
    p.hrange = [38.1 63.5];
    p.mrange = set_available();
    
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
function [n, h, m, o] = convert(p, x)
    
    nmax = p.nrange(end);
    
    n = p.nrange(x(2*nmax+1));
    h(1:n) = p.hrange(x(1:n));
    m(1:n) = p.mrange(x(nmax+1:nmax+n));
    
    if x(2*nmax+2) == 1
        o(1:n) = 1;
        for i = 1:n
            if mod(i, 2) == 0
                o(i) = 0;
            end
        end
    else
        o(1:n) = 0;
        for i = 1:n
            if mod(i, 2) == 0
                o(i) = 1;
            end
        end
    end
    
end
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
function [c] = cost(x)
    
    global p
    [n, h, m] = convert(p, x);
    
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
function [pn] = panelinfo(p, x)
    
    [pn.n, pn.h, pn.m, pn.o] = convert(p, x);
    pn.cost = cost(x);
    pn.wt = load(p, pn.n, pn.h, pn.m);
    pn.con = constraints(p, pn.wt, pn.n, pn.h, pn.m, pn.o);
    
    switch p.Fire
        case -1
            [pn.AF.n, pn.AF.h, pn.AF.m, pn.AF.o, pn.AF.allburnt, pn.AF.xcb] = after_fire_bottom(p, pn.n, pn.h, pn.m, pn.o);
            pn.AF.xct = 0;
        case 1
            [pn.AF.n, pn.AF.h, pn.AF.m, pn.AF.o, pn.AF.allburnt, pn.AF.xct] = after_fire_top(p, pn.n, pn.h, pn.m, pn.o);
            pn.AF.xcb = 0;
        case 0
            [pn.AF.n, pn.AF.h, pn.AF.m, pn.AF.o, pn.AF.allburnt, pn.AF.xcb, pn.AF.xct] = after_fire_both(p, pn.n, pn.h, pn.m, pn.o);
        otherwise
            pn.AF.n = 0; 
            pn.AF.h = 0; 
            pn.AF.m = 0; 
            pn.AF.o = 0; 
            pn.AF.allburnt = 1;
            pn.AF.xcb = 0;
            pn.AF.xct = 0;
    end
    if not(pn.AF.allburnt)
        pn.AF.con = constraints(p, pn.wt, pn.AF.n, pn.AF.h, pn.AF.m, pn.AF.o);
    else
        pn.AF.con = 0;
    end
    
end
function [ineqcon, eqcon] = nonlcon(x)
    
    global p
    [n, h, m, o] = convert(p, x);
    wt = load(p, n, h, m);
    
    [con] = constraints(p, wt, n, h, m, o);
    
    ineqcon(1) = con.M - con.Mmax;
    ineqcon(2) = con.V - con.Vmax;
    ineqcon(3) = con.D - con.Dmax;
    ineqcon(4) = con.L - con.Lmax;
    ineqcon(5) = con.H - con.Hmax;
    
    if n == 9 && o(1) == 1
        ineqcon(6) = 1;
    else
        ineqcon(6) = -1;
    end
    
    switch p.Fire
        case -1
            [n, h, m, o, allburnt] = after_fire_bottom(p, n, h, m, o);
        case 1
            [n, h, m, o, allburnt] = after_fire_top(p, n, h, m, o);
        case 0
            [n, h, m, o, allburnt] = after_fire_both(p, n, h, m, o);
        otherwise
            allburnt = 0;
    end
    
    if not(allburnt)
        [con] = constraints(p, wt, n, h, m, o);
        ineqcon(7) = con.M - con.Mmax;
        ineqcon(8) = con.V - con.Vmax;
    else
        ineqcon(7) = 1;
        ineqcon(8) = 1;
    end
    
    eqcon = [];
    
end

% execution
function [r] = findopt()
    
    global p
    
    nmax = p.nrange(end);
    nvar = 2 * nmax + 2;
    
    lb = ones(nvar, 1);
    ub = ones(nvar, 1);
    ub(1:nmax) = length(p.hrange);
    ub(nmax+1:2*nmax) = length(p.mrange);
    ub(2*nmax+1) = length(p.nrange);
    lb(2*nmax+2) = 0;
    
    intCon = 1:nvar;
    
    options = optimoptions('ga');
    options.Display = 'off';
    
    tic
    s = 0;
    k = 0;
    br = 0;
    
    optnum = 0;
    optcost = Inf;
    optx = zeros(1, nvar);
    costs = NaN(1, p.kmax);
    
    confidence = 0;
    confidence1 = 0;
    confidence2 = 0;
    
    while k < p.kmax && toc < p.tmax
        [x, fval, exitflag] = ga(@cost, nvar, [], [], [], [], lb, ub, @nonlcon, intCon, options);
        if exitflag == 1
            s = s + 1;
            costs(s) = fval;
            if fval < optcost
                optcost = fval;
                optx = x;
                optnum = 1;
            elseif fval == optcost
                optnum = optnum + 1;
            end
            confidence1 = optnum / s;
            confidence2 = exp((-2) * ((mean(costs, 'omitnan') - optcost) / optcost));
            confidence = sqrt(confidence1 * confidence2);
            for i = 1:length(p.smax)
                if s >= p.smax(i) && confidence >= p.cmin(i)
                    br = 1;
                end
            end
        end
        k = k + 1;
        if br == 1
            break;
        end
    end
    
    r.p = p;
    r.s = s;
    r.k = k;
    r.q = s / k;
    r.c = [confidence confidence1 confidence2];
    if s > 0
        r.opt = panelinfo(p, optx);
    else
        r.opt = 0;
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