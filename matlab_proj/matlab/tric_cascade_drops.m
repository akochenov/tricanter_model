function [phi_w_next, phi_w_out] = tric_cascade_drops(phi_w, geo, u, pc, dt)
%TRIC_CASCADE_DROPS Осаждение капель воды в нефтяном кольце. Владелец: B.
%
%   phi_w      — n×J, концентрации с предыдущего шага (СОСТОЯНИЕ)
%   geo        — шина BusTricGeo
%   u, pc      — шины входов и параметров
%   dt         — шаг по времени, с
%   phi_w_next — n×J, концентрации на следующий шаг
%   phi_w_out  — остаточная обводнённость нефти
%
% Кольцо Ro .. r_i, несущая — НЕФТЬ (eta_o из ВХОДОВ, не из параметров:
% сценарий 13.5 меняет её при охлаждении сырья).
%
% Захвата в возврате нет: капли садятся на границу нефть-вода, в кек
% не идут и в баланс массы твёрдого не входят. Этим они отличаются
% от твёрдого в нефти, которое границу пересекает — см.
% tric_cascade_solids_oil.
%
% ВНИМАНИЕ: phi_w0 теперь считается в tric_geometry и знаменатель у него
% eps_o + eps_wd + eps_s_o, а не eps_o + eps_wd. При frac_s_in_oil = 0
% это одно и то же.
%#codegen

if ~geo.oil
    % Кольца нет: либо нефти вообще нет (двухфазный режим), либо r_i
    % прижался к сливу. Разделяться негде — состояние не обновляется,
    % на выход идёт входная концентрация. В двухфазном режиме
    % geo.phi_w0 = 0, поэтому получаются честные нули.
    phi_w_next = phi_w;
    phi_w_out  = geo.phi_w0;
    return
end

[n, J] = tric_dims();

% Капли всегда сигмоидные: pc.rrsb относится только к твёрдому.
[xw, ww] = tric_psd(u.d50, pc.bw, false);

tau_o = (pi * (geo.r_i^2 - u.Ro^2) * geo.Lax / geo.Qo) * ones(n, 1);

% Несущая — нефть. Перепутать с eta_w значит ошибиться в тридцать раз.
k_w = (pc.rho_w - pc.rho_o) * geo.omega^2 / (18 * u.eta_o);

pref_o = geo.r_i^2 / max(geo.r_i^2 - u.Ro^2, 1e-12);

% Опорная концентрация — phi_max (плотная упаковка капель), не phi_ref.
T_w = tric_grade(xw, k_w, tau_o * geo.f_clar, pref_o, phi_w, pc.phi_max);

[~, ~, phi_w_next] = tric_settle_dyn(geo.phi_w0, ww, zeros(n, J), ...
                                     T_w, phi_w, tau_o, dt);

phi_w_out = sum(phi_w_next(end, :));
end
