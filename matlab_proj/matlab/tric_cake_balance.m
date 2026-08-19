function [dm, tr] = tric_cake_balance(Rtr, cap_s, geo, pc)
%TRIC_CAKE_BALANCE Транспорт кека шнеком и баланс масс. Владелец: C.
%
%   Rtr   — n×1, радиус поверхности кека, м
%   E_s_i — n×1, захват твёрдого по ячейкам, доли расхода
%   dm    — n×1, dm/dt, кг/с
%   tr    — n×1, поток кека из ячейки i в сторону выгрузки, кг/с
%
% GOOD: нулевой баланс. Эталон: ref/cake_*.csv (колонки tr, dm).
%
%   tr = phi_sed*rho_s*pi*(Rd^2 - Rtr^2)*u_ax
%   dm = rho_s*Qw*cap_s + [tr(2:end); 0] - tr
% Кек выносится в сторону i = 1, граничное условие tr(n+1) = 0.
%#codegen

n = tric_dims();

tr = pc.phi_sed*pc.rho_s*pi*(pc.Rd^2 - Rtr.^2)*geo.u_ax;
dm = pc.rho_s*geo.Qw*cap_s + [tr(2:end); 0] - tr;
end
