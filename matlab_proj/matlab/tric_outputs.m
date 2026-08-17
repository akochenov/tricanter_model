function y = tric_outputs(m, Rtr, dH, E_s_i, phi_s, phi_w, phi_s_out, phi_w_out, tr, geo, u, pc)
%TRIC_OUTPUTS Наблюдаемые величины и диагностика. Владелец: D.
%
%   phi_s, phi_w — n×J, концентрации ПОСЛЕ обновления
%   y            — шина BusTricOut: cake, streams, quality, flags
%
% ПОКА_ЗАГЛУШКА: заполняет то, что уже посчитано, остальное нулями.
% Эталон: ref_dyn/steady_out_*.csv, trajectory_*.csv.
%
% Формулы:
%   U          = mean(Rd^2 - Rtr.^2) / max(Rd^2 - r_i^2, 1e-12)
%   E_s        = 1 - phi_s_out / phi_s0
%   E_w        = 1 - phi_w_out / phi_w0
%   phi_s_prof = sum(phi_s, 2)    n×1, профиль вдоль машины
%   phi_w_prof = sum(phi_w, 2)    n×1
%   tau_c      = sum(tau_s) — полное время пребывания жидкости
%   mdot_cake  = tr(1)      — уже кг/с, пересчёт не нужен
%   mdot_oil   = Qo * плотность смеси в нефтяном кольце (через phi_w_out)
%   mdot_water = Qw * плотность смеси в водяном кольце
%#codegen

n = tric_dims();

phi_s0 = u.eps_s / (u.eps_wf + u.eps_s);
if geo.has_oil
    phi_w0 = u.eps_wd / (u.eps_o + u.eps_wd);
else
    phi_w0 = 0.0;
end

% --- cake ---
y.cake.Rtr        = Rtr;
y.cake.dH         = dH;
y.cake.m_cake     = m;
y.cake.m_cake_tot = sum(m);
y.cake.tau_s      = zeros(n,1);

% --- streams ---
y.streams.mdot_oil     = 0.0;
y.streams.mdot_water   = 0.0;
y.streams.mdot_cake    = tr(1);
y.streams.x_w_in_oil   = 0.0;
y.streams.x_o_in_water = 0.0;
y.streams.x_s_in_cake  = 0.0;

% --- quality ---
y.quality.U          = 0.0;
y.quality.E_s        = 0.0;
y.quality.E_w        = 0.0;
if phi_s0 > 0
    y.quality.E_s = 1 - phi_s_out / phi_s0;
end
if phi_w0 > 0
    y.quality.E_w = 1 - phi_w_out / phi_w0;
end
y.quality.phi_s_out  = phi_s_out;
y.quality.phi_w_out  = phi_w_out;
y.quality.phi_s_prof = sum(phi_s, 2);
y.quality.phi_w_prof = sum(phi_w, 2);
y.quality.E_s_i      = E_s_i;
y.quality.tau_c      = 0.0;

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

% --- drives ---
y.drives.M1 = 0;
y.drives.P1 = 0;
y.drives.M2 = 0;
y.drives.P2 = 0;
end
