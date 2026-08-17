function [dm, phi_s_next, phi_w_next, y] = tric_rhs(m, phi_s, phi_w, u, pc, dt)
%TRIC_RHS Правая часть модели трикантера, динамический вариант. Владелец: A.
%
% Тело единственного MATLAB Function block. Три состояния:
%   m      — n×1,  масса кека            -> Discrete-Time Integrator
%   phi_s  — n×J,  твёрдое по ячейкам    -> Unit Delay
%   phi_w  — n×J,  капли по ячейкам      -> Unit Delay
%
% Масса кека интегрируется (dm/dt), концентрации обновляются точной
% экспонентой и потому возвращаются как готовое следующее значение,
% а не как производная. Отсюда разные блоки памяти.
%
%   dt — шаг по времени, с. Задаётся как Parameter в Edit Data блока,
%        значение берётся из объекта Ts в словаре. ОБЯЗАН совпадать
%        с Fixed-step size решателя, иначе экспонента посчитается не
%        для того интервала и результат тихо разойдётся с эталоном.
%
% ЭТОТ ФАЙЛ ПРАВИТ ТОЛЬКО A.
%#codegen

geo = tric_geometry(u, pc);

[Rtr, dH] = tric_cake_geom(m, geo, pc);

[cap_s, phi_s_next, phi_s_out] = tric_cascade_solids(phi_s, Rtr, geo, u, pc, dt);
[phi_w_next, phi_w_out]        = tric_cascade_drops(phi_w, geo, u, pc, dt);

[dm, tr] = tric_cake_balance(Rtr, cap_s, geo, pc);

y = tric_outputs(m, Rtr, dH, cap_s, phi_s_next, phi_w_next, ...
                 phi_s_out, phi_w_out, tr, geo, u, pc);
end
