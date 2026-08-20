function geo = tric_geometry(u, pc)
%TRIC_GEOMETRY Алгебраическая преамбула модели. Владелец: A.
%
%   u  — шина входов (BusTricInputs)
%   pc — шина параметров (BusTricParams)
%   geo — шина BusTricGeo
%
% Перенос @property из питон-класса Params (tricanter3_oil.py).
% Эталон: ref_dyn/geometry_*.csv (при frac_s_in_oil = 0).
%
% ОБНОВЛЕНО под трёхпопуляционную модель: твёрдое подачи делится между
% несущими фазами долей frac_s_in_oil, поэтому изменились Qo, Qw и все
% входные концентрации. При frac_s_in_oil = 0 всё вырождается в прежние
% формулы до последнего бита.
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

% --- есть ли нефтяное КОЛЬЦО, а не только нефть ---
% Отдельный флаг: нефть в подаче может быть, а кольца уже нет, если
% r_i прижался к нефтяному сливу. Тогда в нефтяном кольце ничего
% не оседает, и обе «нефтяные» популяции идут насквозь.
oil = has_oil && (r_i - u.Ro > 1e-6);

% --- деление твёрдого подачи между несущими фазами ---
% Эмпирический вход (смачиваемость, история эмульсии), из состава
% не выводится. Без нефти диспергировать не во что.
if has_oil
    f_so = min(max(u.frac_s_in_oil, 0), 1);
else
    f_so = 0;
end
eps_s_o = u.eps_s * f_so;           % твёрдое, пришедшее внутри нефти
eps_s_w = u.eps_s * (1 - f_so);     % твёрдое, пришедшее внутри воды

% --- расходы по кольцам ---
% Твёрдое, диспергированное в нефти, едет с нефтяным расходом, поэтому
% входит в Qo, а не в Qw. Переток через границу расходы не меняет:
% дисперсии разбавлены.
Qo = u.Q * (u.eps_o + u.eps_wd + eps_s_o);
Qw = u.Q * (u.eps_wf + eps_s_w);

% --- входные доли дисперсных фаз в своих несущих ---
d_w = u.eps_wf + eps_s_w;
if d_w > 0
    phi_s0 = eps_s_w / d_w;         % твёрдое в воде
else
    phi_s0 = 0;
end

d_o = u.eps_o + u.eps_wd + eps_s_o;
if d_o > 0
    phi_so0 = eps_s_o / d_o;        % твёрдое в нефти
    phi_w0  = u.eps_wd / d_o;       % капли воды в нефти
else
    phi_so0 = 0;
    phi_w0  = 0;
end

% --- сборка выходной шины. Порядок обязан совпасть с BusTricGeo ---
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
geo.oil     = oil;
geo.f_so    = f_so;
geo.phi_s0  = phi_s0;
geo.phi_so0 = phi_so0;
geo.phi_w0  = phi_w0;
end
