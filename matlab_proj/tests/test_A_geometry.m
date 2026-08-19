function test_A_geometry(ddPath)
%TEST_A_GEOMETRY Тест участника A. Запуск: test_A_geometry
%
% Проверяет tric_geometry против эталона из питона и то, что
% tric_rhs проходит насквозь через заглушки без ошибок.
%
% Данные читаются из словаря, а не из base workspace.

if nargin < 1 || isempty(ddPath), ddPath = 'tricanter_data.sldd'; end
tol = 1e-10;

[pc, u, Ts] = read_dd(ddPath);

%% 1. Геометрия против эталона (трёхфазный режим)
r   = ref_load('geometry_3ph');
geo = tric_geometry(u, pc);

fprintf('Геометрия 3ph:\n');
chk_named('omega',  geo.omega,  r.omega,  tol);
chk_named('beta',   geo.beta,   r.beta,   tol);
chk_named('Lax',    geo.Lax,    r.Lax,    tol);
chk_named('L_cone', geo.L_cone, r.L_cone, tol);
chk_named('f_clar', geo.f_clar, r.f_clar, tol);
chk_named('u_ax',   geo.u_ax,   r.u_ax,   tol);
chk_named('r_i',    geo.r_i,    r.r_i,    tol);
chk_named('Qo',     geo.Qo,     r.Qo,     tol);
chk_named('Qw',     geo.Qw,     r.Qw,     tol);

%% 2. Сведение к двухфазному
% При eps_o = eps_wd = 0 нефтяного кольца нет, и внутренний радиус
% обязан совпасть с радиусом слива воды до машинной точности.
% Расхождение здесь означает ошибку в балансе давлений.
u2 = u;  u2.eps_o = 0;  u2.eps_wd = 0;  u2.eps_wf = 0.98;
geo2 = tric_geometry(u2, pc);
fprintf('Сведение к 2ph:\n');
chk_named('r_i = Rw', geo2.r_i, u2.Rw, 1e-14);

%% 3. Прогон tric_rhs через заглушки
[n, J] = tric_dims();
m   = zeros(n, 1);
phi = zeros(n, J);
% Четыре состояния: масса кека и три дисперсные популяции.
[dm, phi_s_next, phi_so_next, phi_w_next, y] = tric_rhs(m, phi, phi, phi, u, pc, Ts);

fprintf('Сквозной прогон tric_rhs:\n');
assert(isequal(size(dm), [n 1]),         'dm: ожидалось %dx1', n);
assert(isequal(size(phi_s_next), [n J]), 'phi_s_next: ожидалось %dx%d', n, J);
assert(isequal(size(phi_so_next), [n J]), 'phi_so_next: ожидалось %dx%d', n, J);
assert(isequal(size(phi_w_next), [n J]), 'phi_w_next: ожидалось %dx%d', n, J);
assert(all(isfinite(dm)),                'dm содержит NaN или Inf');
assert(isequal(size(y.cake.Rtr), [n 1]), 'y.cake.Rtr: ожидалось %dx1', n);
assert(isequal(size(y.quality.E_s_i), [n 1]), 'y.quality.E_s_i: ожидалось %dx1', n);
fprintf('  ок    размеры верны, NaN нет\n');

fprintf('\nПроверки A пройдены.\n');
end


% =====================================================================
function [pc, u, Ts] = read_dd(ddPath)
%READ_DD Достаёт значения из словаря данных.
dd  = Simulink.data.dictionary.open(ddPath);
sec = dd.getSection('Design Data');
pc  = unwrap(sec.getEntry('TricParams').getValue());
u   = unwrap(sec.getEntry('TricInputs0').getValue());
Ts  = unwrap(sec.getEntry('Ts').getValue());
dd.close();
end

function v = unwrap(v)
%UNWRAP Разворачивает Simulink.Parameter до самого значения.
if isa(v, 'Simulink.Parameter'), v = v.Value; end
end

function chk_named(name, got, want, tol)
err = abs(got - want) / max(abs(want), eps);
if err > tol
    error('  ПРОВАЛ %s: получено %.17g, эталон %.17g (отн. %.3g)', ...
          name, got, want, err);
end
fprintf('  ок    %-12s %.10g\n', name, got);
end
