function [Rtr, dH] = tric_cake_geom(m, geo, pc)
%TRIC_CAKE_GEOM Геометрия кека из его массы. Владелец: C.
%
%   m   — n×1, масса кека по ячейкам, кг
%   Rtr — n×1, радиус поверхности кека, м
%   dH  — n×1, толщина слоя осадка, м
%
% ЗАГЛУШКА: пустая машина, Rtr = Rd, dH = 0. Заглушка ОБЯЗАНА возвращать
% именно Rd, а не нули: при Rtr = 0 время пребывания уйдёт в минус и B
% с D будут отлаживать не свой код.
%
% Эталон: ref/cake_empty_*.csv, ref/cake_steady_*.csv (колонки Rtr, dH).
%
%   Rtr = sqrt(max(Rd^2 - m/(rho_s*phi_sed*pi*Lax), r_i^2))
%   dH  = Rd - Rtr
%#codegen

n = tric_dims();

Rtr = pc.Rd * ones(n,1);
dH  = zeros(n,1);
end
