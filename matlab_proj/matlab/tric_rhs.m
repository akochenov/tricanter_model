function [dm, phi_s_next, phi_so_next, phi_w_next, y] = tric_rhs(m, phi_s, phi_so, phi_w, u, pc, dt)
%TRIC_RHS Правая часть модели трикантера, динамический вариант. Владелец: A.
%
% Тело единственного MATLAB Function block. ЧЕТЫРЕ состояния:
%   m       — n×1,  масса кека                 -> Discrete-Time Integrator
%   phi_s   — n×J,  твёрдое в воде             -> Unit Delay
%   phi_so  — n×J,  твёрдое в нефти            -> Unit Delay   <- НОВОЕ
%   phi_w   — n×J,  капли воды в нефти         -> Unit Delay
%
% Масса кека интегрируется (dm/dt), концентрации обновляются точной
% экспонентой и потому возвращаются как готовое следующее значение.
%
%   dt — шаг по времени, с. ОБЯЗАН совпадать с Fixed-step size решателя,
%        иначе экспонента посчитается не для того интервала.
%
% ПОРЯДОК ВЫЗОВОВ ЗНАЧИМ. Сначала нефтяное кольцо целиком, потом водяное:
% твёрдое, захваченное на границу r_i, пересекает её и достаётся водяному
% кольцу ТОЙ ЖЕ ячейки. Обратной связи «вода -> нефть» нет, поэтому кольца
% и разделяются на два независимых вызова.
%
% ЭТОТ ФАЙЛ ПРАВИТ ТОЛЬКО A.
%#codegen

geo = tric_geometry(u, pc);

[Rtr, dH] = tric_cake_geom(m, geo, pc);

% --- нефтяное кольцо: твёрдое и капли, обе популяции на границу r_i ---
[moved, phi_so_next, phi_so_out] = tric_cascade_solids_oil(phi_so, geo, u, pc, dt);
[phi_w_next, phi_w_out]          = tric_cascade_drops(phi_w, geo, u, pc, dt);

% --- водяное кольцо: своя подача плюс переток из нефти ---
[cap_s, phi_s_next, phi_s_out] = tric_cascade_solids(phi_s, moved, Rtr, geo, u, pc, dt);

[dm, tr] = tric_cake_balance(Rtr, cap_s, geo, pc);

y = tric_outputs(m, Rtr, dH, cap_s, phi_s_next, phi_so_next, phi_w_next, ...
                 phi_s_out, phi_so_out, phi_w_out, tr, geo, u, pc);
end
