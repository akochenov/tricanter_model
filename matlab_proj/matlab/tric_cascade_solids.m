function [cap_s, phi_s_next, phi_s_out] = tric_cascade_solids(phi_s, Rtr, geo, u, pc, dt)
%TRIC_CASCADE_SOLIDS Осаждение твёрдого в водяном кольце. Владелец: B.
%
%   phi_s      — n×J, концентрации с предыдущего шага (СОСТОЯНИЕ)
%   Rtr        — n×1, радиус поверхности кека, м
%   geo        — шина BusTricGeo
%   u, pc      — шины входов и параметров
%   dt         — шаг по времени, с
%   cap_s      — n×1, захват по ячейкам в долях расхода
%   phi_s_next — n×J, концентрации на следующий шаг
%   phi_s_out  — скаляр, доля твёрдого в фугате
%
% Обёртка: собирает коэффициенты, зовёт tric_grade и tric_settle_dyn.
% Питон-оригинал — строки 208-222 simulate_dyn.
% Эталон: ref_dyn/step1_sca_*.csv (k_s, pref_s, phi_s0).
%
% Кольцо r_i .. Rtr, несущая фаза — ВОДА (eta_w из параметров).
% Для капель несущая нефть, вязкости отличаются в тридцать раз.
%#codegen

% Классы твёрдого. rrsb — параметр, у капель он всегда false.
[xs, ws] = tric_psd(u.x50, pc.bs, pc.rrsb);

% Время пребывания в ячейке: зависит и от ячейки, и от намытого кека.
% Это единственный канал обратной связи «кек -> осаждение».
tau_s = pi * (Rtr.^2 - geo.r_i^2) * geo.Lax / geo.Qw;

% Стоксова константа без x^2. Несущая — вода.
k_s = (pc.rho_s - pc.rho_w) * geo.omega^2 / (18 * pc.eta_w);

% КРИТИЧНО: префактор по Rd, а не по Rtr. В печатных формулах Gleiss
% в обоих местах стоит Rtr, и это ошибка публикации: с Rtr модель
% теряет самоограничение и U уходит в единицу.
pref_s = pc.Rd^2 / max(pc.Rd^2 - geo.r_i^2, 1e-12);

% Знаменатель — свободная вода плюс твёрдое (eps_wf, не eps_o).
phi_s0 = u.eps_s / (u.eps_wf + u.eps_s);

% В grade tau идёт ДОМНОЖЕННЫМ на f_clar (конус тоже осветляет),
% в settle_dyn — голым: экспонента затухания идёт по реальному
% времени пребывания. В CSV записан голый.
T_s = tric_grade(xs, k_s, tau_s * geo.f_clar, pref_s, phi_s, pc.phi_ref);

[cap_s, phi_s_next] = tric_settle_dyn(phi_s0, ws, T_s, phi_s, tau_s, dt);

% Фугат — то, что вышло из последней ячейки.
phi_s_out = sum(phi_s_next(end, :));
end
