function [dm, tr] = tric_cake_balance(Rtr, cap_s, geo, pc)
%TRIC_CAKE_BALANCE Транспорт кека шнеком и баланс масс. Владелец: C.
%
%   Rtr   — n×1, радиус поверхности кека, м
%   cap_s — n×1, захват твёрдого по ячейкам, доли расхода
%   dm    — n×1, dm/dt, кг/с
%   tr    — n×1, поток кека из ячейки i в сторону выгрузки, кг/с
%
% ЗАГЛУШКА: нулевой баланс. Эталон: ref/cake_*.csv (колонки tr, dm).
%
%   tr = phi_sed*rho_s*pi*(Rd^2 - Rtr^2)*u_ax
%   dm = rho_s*Qw*cap_s + [tr(2:end); 0] - tr
% Кек выносится в сторону i = 1, граничное условие tr(n+1) = 0.
%#codegen

n = tric_dims();

dm = zeros(n,1);
tr = zeros(n,1);
end
