function y = tric_outputs(m, Rtr, dH, cap_s, phi_s, phi_so, phi_w, ...
                          phi_s_out, phi_so_out, phi_w_out, tr, geo, u, pc)
%TRIC_OUTPUTS Наблюдаемые величины и диагностика. Владелец: D.
%
%   cap_s              — n×1, захват на кек по ячейкам (НЕ E_s_i!)
%   phi_s, phi_so, phi_w — n×J, концентрации ПОСЛЕ обновления
%   y                  — шина BusTricOut: cake, streams, quality, flags
%
% ПОКА_ЗАГЛУШКА: заполняет то, что уже посчитано, остальное нулями.
% Сигнатура и новые поля обновлены под трёхпопуляционную модель
% (добавилось твёрдое в нефти); сама реализация за D.
% Эталон: ref_dyn/steady_out_*.csv, trajectory_*.csv.
%
% Формулы:
%   U          = mean(Rd^2 - Rtr.^2) / max(Rd^2 - r_i^2, 1e-12)
%   E_s        = 1 - phi_s_out / phi_s0      осветление ВОДЫ
%   E_so       = 1 - phi_so_out / phi_so0    осветление НЕФТИ
%   E_w        = 1 - phi_w_out / phi_w0
%   E_s_tot    = 1 - (Qw*phi_s_out + Qo*phi_so_out) / (Q*eps_s)
%   bsw        = phi_w_out + phi_so_out      вода и мехпримеси в нефти
%
%   phi_s_prof = sum(phi_s, 2)     n×1, профиль вдоль машины
%   phi_so_prof= sum(phi_so, 2)    n×1
%   phi_w_prof = sum(phi_w, 2)     n×1
%   q_cake     = tr(1) / rho_s     твёрдое в выгрузке, м3/с
%   tau_c      = sum(tau_s) — полное время пребывания жидкости
%   mdot_cake  = tr(1)      — уже кг/с, пересчёт не нужен
%   mdot_oil   = Qo * плотность смеси в нефтяном кольце
%   mdot_water = Qw * плотность смеси в водяном кольце
%
% E_s при заметном frac_s_in_oil может стать ОТРИЦАТЕЛЬНЫМ: нефть
% досыпает твёрдого в воду, и водяной поток осветляется хуже, чем
% его собственный вход. Честная общая величина — E_s_tot.
%#codegen

n = tric_dims();

% Входные доли теперь считает tric_geometry: знаменатели зависят
% от деления твёрдого между фазами, руками их не повторять.
phi_s0  = geo.phi_s0;
phi_so0 = geo.phi_so0;
phi_w0  = geo.phi_w0;

% --- cake ---
y.cake.Rtr        = Rtr;
y.cake.dH         = dH;
y.cake.m_cake     = m;
y.cake.m_cake_tot = sum(m);
y.cake.tau_s      = zeros(n,1);
y.cake.q_cake     = tr(1) / pc.rho_s;      % твёрдое в выгрузке, м3/с

% --- streams ---
y.streams.mdot_oil     = 0.0;
y.streams.mdot_water   = 0.0;
y.streams.mdot_cake    = tr(1);
y.streams.x_w_in_oil   = 0.0;
y.streams.x_o_in_water = 0.0;
y.streams.x_s_in_cake  = 0.0;

% --- quality ---
% ПОРЯДОК ПРИСВАИВАНИЯ = порядок полей структуры, и он обязан совпасть
% с порядком элементов в BusQuality. Новые поля идут В КОНЦЕ, после
% tau_c, а не по смыслу — см. WIRING_OIL.md.
% При нулевом входе эффективность по соглашению нуль, а не NaN.
y.quality.U          = 0.0;
y.quality.E_s        = 0.0;
if phi_s0 > 0
    y.quality.E_s = 1 - phi_s_out / phi_s0;
end
y.quality.E_w        = 0.0;
if phi_w0 > 0
    y.quality.E_w = 1 - phi_w_out / phi_w0;
end
y.quality.phi_s_out  = phi_s_out;
y.quality.phi_w_out  = phi_w_out;
y.quality.phi_s_prof = sum(phi_s, 2);
y.quality.phi_w_prof = sum(phi_w, 2);
% E_s_i — доля осевшего в ячейке от вошедшего в неё. Это НЕ cap_s:
% cap_s падает вдоль машины в тысячи раз, E_s_i — в десяток.
y.quality.E_s_i      = zeros(n,1);
y.quality.tau_c      = 0.0;

% --- новые поля: третья популяция ---
y.quality.E_so       = 0.0;
if phi_so0 > 0
    y.quality.E_so = 1 - phi_so_out / phi_so0;
end
y.quality.E_s_tot    = 0.0;
if u.eps_s > 0
    esc = geo.Qw * phi_s_out + geo.Qo * phi_so_out;    % унос твёрдого, м3/с
    y.quality.E_s_tot = 1 - esc / (u.Q * u.eps_s);
end
y.quality.phi_so_out  = phi_so_out;
y.quality.bsw         = phi_w_out + phi_so_out;
y.quality.phi_so_prof = sum(phi_so, 2);

% --- flags ---
flags = tric_limits(geo, u, pc);
y.flags.bal_resid        = 0.0;
y.flags.M_scroll         = 0.0;
y.flags.f_oil_collapse   = flags.oil_collapse;
y.flags.f_water_collapse = flags.water_collapse;
y.flags.f_overload       = false;
y.flags.f_overfill       = false;
y.flags.f_M1_limit = false;
y.flags.f_M2_limit = false;
% Проверки состава подачи — из tric_limits, тоже в конце шины.
y.flags.f_feed_sum   = flags.feed_sum;
y.flags.f_frac_range = flags.frac_range;

% --- drives ---
y.drives.M1 = 0;
y.drives.P1 = 0;
y.drives.M2 = 0;
y.drives.P2 = 0;
end
