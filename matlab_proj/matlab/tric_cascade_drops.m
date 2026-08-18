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
% Питон-оригинал — строки 225-241 simulate_dyn.
% Эталон: ref_dyn/step1_phiw_next_*.csv, steady_out_*.csv (phi_w_out).
%
% Кольцо Ro .. r_i, несущая фаза — НЕФТЬ (eta_o из ВХОДОВ, не из
% параметров: сценарий 13.5 меняет её при охлаждении сырья).
%
% Захвата в возврате нет: капли садятся на границу нефть-вода,
% в кек не идут и в баланс массы не входят.
%#codegen

[n, J] = tric_dims();

% Три ветки вместо питоновских двух с тернарником. Для codegen каждая
% обязана присвоить phi_w_next одного и того же размера n×J.
if geo.has_oil && (geo.r_i - u.Ro > 1e-6)

    % Классы капель. rrsb здесь всегда false, форма распределения
    % сигмоидная — параметр pc.rrsb относится только к твёрдому.
    [xw, ww] = tric_psd(u.d50, pc.bw, false);

    % Кольцо нефти от кека не зависит, поэтому tau одинаково во всех
    % ячейках. Раздуть до n×1 всё равно обязательно: grade и settle_dyn
    % ждут вектор по ячейкам.
    tau_o = (pi * (geo.r_i^2 - u.Ro^2) * geo.Lax / geo.Qo) * ones(n, 1);

    % Несущая — нефть. eta_o, не eta_w: разница в тридцать раз.
    k_w = (pc.rho_w - pc.rho_o) * geo.omega^2 / (18 * u.eta_o);

    % Префактор через r_i, а не через Rd: капли оседают на границу
    % нефть-вода, внешняя стенка кольца для них — это r_i.
    pref_w = geo.r_i^2 / max(geo.r_i^2 - u.Ro^2, 1e-12);

    phi_w0 = u.eps_wd / (u.eps_o + u.eps_wd);

    % Опорная концентрация — phi_max (плотная упаковка капель),
    % а не phi_ref. phi_ref относится к твёрдому в воде.
    T_w = tric_grade(xw, k_w, tau_o * geo.f_clar, pref_w, phi_w, pc.phi_max);

    [~, phi_w_next] = tric_settle_dyn(phi_w0, ww, T_w, phi_w, tau_o, dt);

    phi_w_out = sum(phi_w_next(end, :));

elseif geo.has_oil

    % Нефть есть, а кольца нет: r_i прижался к нефтяному сливу.
    % Разделяться негде, состояние не обновляется вообще,
    % на выход идёт входная концентрация.
    phi_w_next = phi_w;
    phi_w_out  = u.eps_wd / (u.eps_o + u.eps_wd);

else

    % Двухфазный режим: нефти нет, капель нет.
    phi_w_next = zeros(n, J);
    phi_w_out  = 0;

end
end
