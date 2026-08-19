function [cap_s, phi_s_next, phi_s_out] = tric_cascade_solids(phi_s, moved, Rtr, geo, u, pc, dt)
%TRIC_CASCADE_SOLIDS Осаждение твёрдого в водяном кольце. Владелец: B.
%
%   phi_s      — n×J, концентрации с предыдущего шага (СОСТОЯНИЕ)
%   moved      — n×J, поклассовый переток из нефтяного кольца
%                (из tric_cascade_solids_oil), в долях НЕФТЯНОГО расхода
%   Rtr        — n×1, радиус поверхности кека, м
%   geo, u, pc — шины
%   dt         — шаг по времени, с
%   cap_s      — n×1, захват на кек по ячейкам, в долях водяного расхода
%   phi_s_next — n×J, концентрации на следующий шаг
%   phi_s_out  — доля твёрдого в фугате
%
% Кольцо r_i .. Rtr, несущая — ВОДА (eta_w из параметров).
%
% В это кольцо впадает твёрдое, перешедшее из нефти: обе жидкости текут
% соосно к сливам, поэтому захваченное на границу достаётся водяному
% кольцу ТОЙ ЖЕ ячейки и дальше идёт обычным каскадом в кек. Слой на
% границе раздела (rag layer) не моделируется.
%
% ЛОВУШКИ. Префактор pref_s по Rd, а не по Rtr: с Rtr модель теряет
% самоограничение и U уходит в единицу. В tric_grade tau идёт домноженным
% на f_clar, в tric_settle_dyn — голым.
%#codegen

[n, J] = tric_dims();

% Классы общие с нефтяной популяцией — переток поклассовый.
[xs, ws, ~] = tric_psd_solids(u, pc);

% Время пребывания зависит и от ячейки, и от намытого кека.
% Это единственный канал обратной связи «кек -> осаждение».
tau_s = pi * (Rtr.^2 - geo.r_i^2) * geo.Lax / geo.Qw;

k_s    = (pc.rho_s - pc.rho_w) * geo.omega^2 / (18 * pc.eta_w);
pref_s = pc.Rd^2 / max(pc.Rd^2 - geo.r_i^2, 1e-12);

% --- переток из нефти, пересчитанный в доли водяного расхода ---
if geo.oil
    q = geo.Qo / geo.Qw;
else
    q = 0;
end
inj = moved * q;

% Водяное кольцо схлопнулось: граница раздела совпала с поверхностью
% кека, и перетекшему твёрдому осаждаться уже негде — оно попадает
% на кек сразу, минуя каскад.
collapsed = (Rtr - geo.r_i) < 1e-6;
cap_extra = sum(inj, 2) .* collapsed;
inj       = inj .* (1 - collapsed);          % обнулить строки схлопнувшихся

T_s = tric_grade(xs, k_s, tau_s * geo.f_clar, pref_s, phi_s, pc.phi_ref);

if pc.inject_at_interface
    % Плуг: перетекшее в каскад НЕ подмешивается, оно идёт отдельной
    % непереме́шанной фракцией от границы r_i. Оценка снизу для захвата.
    [cap, ~, phi_s_next] = tric_settle_dyn(geo.phi_s0, ws, zeros(n, J), ...
                                           T_s, phi_s, tau_s, dt);

    plug_q  = zeros(1, J);
    lr0     = log(geo.r_i);
    plug_lr = lr0 * ones(1, J);
    plug_out = 0;
    x2      = (xs.^2).';                      % 1×J

    for c = 1:n
        if c == 1
            cin = (geo.phi_s0 * ws).';
        else
            cin = phi_s_next(c-1, :);
        end
        % Стеснение для плуга считается по своему сечению: содержимое
        % ячейки плюс сам плуг.
        R = max(0, 1 - (sum(cin) + plug_out) / pc.phi_ref) ^ 4.65;
        travel = k_s * x2 * R * tau_s(c) * geo.f_clar;
        [cp, plug_q, plug_lr] = tric_plug_step(plug_q, plug_lr, inj(c,:), ...
                                               travel, Rtr(c), lr0);
        cap(c)   = cap(c) + cp;
        plug_out = sum(plug_q);
    end

    cap_s     = cap + cap_extra;
    phi_s_out = sum(phi_s_next(end, :)) + plug_out;
else
    % Основной вариант: перетекшее домешивается равномерно по сечению
    % кольца — то же допущение, на котором стоит вся ячеечная модель.
    [cap, ~, phi_s_next] = tric_settle_dyn(geo.phi_s0, ws, inj, ...
                                           T_s, phi_s, tau_s, dt);
    cap_s     = cap + cap_extra;
    phi_s_out = sum(phi_s_next(end, :));
end
end
