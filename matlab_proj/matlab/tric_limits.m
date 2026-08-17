function flags = tric_limits(geo, u, pc)
%TRIC_LIMITS Проверка существования режима. Владелец: D.
%   ПОКА_ЗАГЛУШКА: надо потестировать то, как че работает
%   flags.oil_collapse   — нефтяное кольцо схлопнулось: вода в нефтяной слив
%   flags.water_collapse — водяное кольцо схлопнулось: нефть в водяной слив
%   flags.near_limit     — режим на грани (кольцо тоньше 1 мм)
%
% Питон-оригинал — check_limits(). Баланс чувствителен: сдвиг Ro на 2 мм
% двигает r_i на 7-10 мм.
%#codegen

flags.oil_collapse   = false;
flags.water_collapse = false;
flags.near_limit     = false;

if ~geo.has_oil
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
