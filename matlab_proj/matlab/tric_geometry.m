function geo = tric_geometry(u, pc)
%TRIC_GEOMETRY Алгебраическая преамбула модели. Владелец: A. ГОТОВО.
%
%   u  — шина входов (BusTricInputs): Q, dn, C, Ro, Rw, eps_*, x50, d50
%   pc — шина параметров (BusTricParams): Rd, Wsc, Lpond, alpha, Lsep,
%        mu_s, u_conv, rho_*, ...
%   geo — шина BusTricGeo
%
% Перенос @property из питон-класса Params. Эталон: ref/geometry_*.csv
%#codegen

G = 9.81;
n = tric_dims();

% --- вращение и геометрия шнека ---
omega = sqrt(u.C * G / pc.Rd);
beta  = atan2(pc.Wsc, 2*pi*pc.Rd);

% --- осевая длина ячейки ---
Lax = pc.Lpond / n;

% --- вклад конуса в осветление, в эквивалентных метрах цилиндра ---
% захват на осевой срез ~ R(x)^2, отсюда int R^2 dx / Rd^2
t      = tan(pc.alpha * pi/180);
L_cone = (pc.Rd^3 - u.Rw^3) / (3*t) / pc.Rd^2;

if pc.Lsep > 0
    L_eff = pc.Lsep;
else
    L_eff = pc.Lpond + L_cone;
end
f_clar = L_eff / pc.Lpond;

% --- осевая скорость выноса кека шнеком ---
kappa = atan(pc.mu_s) + beta;
eff   = 1 / (1 + tan(beta)*tan(kappa));
u_ax  = pc.Wsc * u.dn/60 * eff * pc.u_conv;

% --- есть ли нефтяная фаза ---
has_oil = (u.eps_o + u.eps_wd) > 0;

% --- граница нефть-вода из баланса давлений двух сливов ---
% ВНИМАНИЕ: r_i не зависит от кека и не является состоянием.
if has_oil
    r2  = (pc.rho_w*u.Rw^2 - pc.rho_o*u.Ro^2) / (pc.rho_w - pc.rho_o);
    r_i = sqrt(max(r2, 0));
    r_i = min(max(r_i, u.Ro), pc.Rd);
else
    r_i = u.Rw;
end

% --- расходы по кольцам (капли разбавлены) ---
Qo = u.Q * (u.eps_o + u.eps_wd);
Qw = u.Q * (u.eps_wf + u.eps_s);

% --- сборка выходной шины ---
geo.omega   = omega;
geo.beta    = beta;
geo.Lax     = Lax;
geo.L_cone  = L_cone;
geo.f_clar  = f_clar;
geo.u_ax    = u_ax;
geo.r_i     = r_i;
geo.Qo      = Qo;
geo.Qw      = Qw;
geo.has_oil = has_oil;
end
