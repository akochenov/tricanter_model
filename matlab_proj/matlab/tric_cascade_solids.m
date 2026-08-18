function [cap_s, phi_s_next, phi_s_out] = tric_cascade_solids(phi_s, Rtr, geo, u, pc, dt)
%TRIC_CASCADE_SOLIDS Осаждение твёрдого в водяном кольце. Владелец: B.
%
%   phi_s      — n×J, концентрации с предыдущего шага (СОСТОЯНИЕ)
%   Rtr        — n×1, радиус поверхности кека, м
%   dt         — шаг по времени, с
%   E_s_i      — n×1, захват по ячейкам
%   phi_s_next — n×J, концентрации на следующий шаг
%   phi_s_out  — скаляр, доля твёрдого в фугате = sum(phi_s_next(end,:))
%
% ПОКА_ЗАГЛУШКА: захвата нет.
% Эталон: ref_dyn/step1_vec_*.csv (E_s_i), step1_sca_*.csv (k_s, pref_s).
%
% Кольцо r_i .. Rd, несущая — вода:
%   tau_s  = pi*(Rtr.^2 - r_i^2)*Lax / Qw        <- зависит от ячейки И от кека
%   k_s    = (rho_s - rho_w)*omega^2 / (18*eta_w)
%   pref_s = Rd^2 / max(Rd^2 - r_i^2, 1e-12)     <- Rd, НЕ Rtr
%   phi_s0 = eps_s / (eps_wf + eps_s)
% Затем tric_grade(xs, k_s, tau_s*f_clar, pref_s, phi_s, phi_ref)
% и tric_settle_dyn(phi_s0, ws, T_s, phi_s, tau_s, dt).
%
% Классы xs, ws — из tric_psd(u.x50, pc.bs, pc.rrsb).
%#codegen

[n, J]     = tric_dims();
cap_s      = zeros(n, 1);
phi_s_next = phi_s;
phi_s_out  = u.eps_s / (u.eps_wf + u.eps_s);
end
