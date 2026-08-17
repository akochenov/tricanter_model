function [phi_w_next, phi_w_out] = tric_cascade_drops(phi_w, geo, u, pc, dt)
%TRIC_CASCADE_DROPS Осаждение капель воды в нефтяном кольце. Владелец: B.
%
%   phi_w      — n×J, концентрации с предыдущего шага (СОСТОЯНИЕ)
%   phi_w_next — n×J, концентрации на следующий шаг
%   phi_w_out  — остаточная обводнённость нефти = sum(phi_w_next(end,:))
%
% ЗАГЛУШКА: разделения нет.
% Эталон: ref_dyn/step1_phiw_next_*.csv, steady_out_*.csv (phi_w_out).
%
% Кольцо Ro .. r_i, несущая — нефть:
%   tau_o  = pi*(r_i^2 - Ro^2)*Lax / Qo          <- ОДИНАКОВО для всех ячеек
%   k_w    = (rho_w - rho_o)*omega^2 / (18*eta_o)  <- eta_o из u, не из pc
%   pref_w = r_i^2 / max(r_i^2 - Ro^2, 1e-12)
%   phi_w0 = eps_wd / (eps_o + eps_wd)
%   опорная доля — phi_max (0.64), а не phi_ref
%
% Три ветки, как в питоне:
%   has_oil && (r_i - Ro > 1e-6)  -> полный расчёт
%   has_oil                       -> кольца нет, phi_w_out = phi_w0
%   иначе                         -> всё по нулям (двухфазный режим)
%#codegen

[n, J]     = tric_dims();
phi_w_next = phi_w;

if geo.has_oil
    phi_w_out = u.eps_wd / (u.eps_o + u.eps_wd);
else
    phi_w_out = 0.0;
end
end
