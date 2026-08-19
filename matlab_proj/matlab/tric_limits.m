function flags = tric_limits(geo, u, pc)
%TRIC_LIMITS Проверка существования режима. Владелец: D.
%   ПОКА_ЗАГЛУШКА: надо потестировать то, как че работает
%   flags.oil_collapse   — нефтяное кольцо схлопнулось: вода в нефтяной слив
%   flags.water_collapse — водяное кольцо схлопнулось: нефть в водяной слив
%   flags.near_limit     — режим на грани (кольцо тоньше 1 мм)
%   flags.feed_sum       — доли подачи не дают единицу
%   flags.frac_range     — frac_s_in_oil вне [0,1], зажимается в geometry
%
% Питон-оригинал — check_limits(). Баланс чувствителен: сдвиг Ro на 2 мм
% двигает r_i на 7-10 мм.
%#codegen

flags.oil_collapse   = false;
flags.water_collapse = false;
flags.near_limit     = false;
flags.feed_sum       = false;
flags.frac_range     = false;

% Доли подачи обязаны давать единицу. Твёрдое делится между несущими
% фазами отдельным входом frac_s_in_oil, поэтому состав перебалансировать
% не нужно и сумма остаётся инвариантом.
if abs((u.eps_o + u.eps_wd + u.eps_wf + u.eps_s) - 1.0) > 1e-9
    flags.feed_sum = true;
end

if u.frac_s_in_oil < 0 || u.frac_s_in_oil > 1
    flags.frac_range = true;
end

if ~geo.has_oil
    % Нефти нет — диспергировать твёрдое не во что, frac_s_in_oil
    % игнорируется, всё твёрдое считается пришедшим с водой.
    return
end

if geo.r_i <= u.Ro + 1e-6
    flags.oil_collapse = true;
elseif geo.r_i - u.Ro < 1e-3
    flags.near_limit = true;
end

if geo.r_i >= pc.Rd - 1e-6
    flags.water_collapse = true;
end
end
