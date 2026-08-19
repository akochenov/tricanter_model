function test_B_oil()
%TEST_B_OIL Тест третьей популяции: твёрдое, пришедшее внутри нефти.
%   Запуск: test_B_oil
%
% Эталонов на эту физику НЕТ — она появилась позже, чем выгружен ref_dyn/.
% Поэтому проверяются свойства, а не числа. Главное из них — вырождение:
% при frac_s_in_oil = 0 модель обязана совпасть с прежней.

fprintf('\n=== твёрдое в нефти ===\n');

[u, pc, dt] = ref_case('3ph');       % ref_case зануляет frac_s_in_oil
[n, J] = tric_dims();
v = ref_load('step1_vec_3ph');
Z = zeros(n, J);

%% --- вырождение при frac_s_in_oil = 0 ---
geo0 = tric_geometry(u, pc);
[moved0, so_next0, so_out0] = tric_cascade_solids_oil(Z, geo0, u, pc, dt);

assert(geo0.f_so == 0,        'f_so обязан быть нулём');
assert(geo0.phi_so0 == 0,     'phi_so0 обязан быть нулём');
assert(all(moved0(:) == 0),   'при frac_s_in_oil = 0 перетока быть не может');
assert(so_out0 == 0,          'твёрдого в нефти нет');
assert(isequal(so_next0, Z),  'состояние не должно меняться');
fprintf('  вырождение при frac = 0: ок\n');

% Расходы и доли обязаны совпасть со старыми формулами
g = ref_load('geometry_3ph');
fprintf('  Qo  vs эталон: ');  chk(geo0.Qo, g.Qo);
fprintf('  Qw  vs эталон: ');  chk(geo0.Qw, g.Qw);
s = ref_load('step1_sca_3ph');
fprintf('  phi_s0 vs эт.: ');  chk(geo0.phi_s0, s.phi_s0);
fprintf('  phi_w0 vs эт.: ');  chk(geo0.phi_w0, s.phi_w0);

%% --- рабочий режим: frac_s_in_oil > 0 ---
u2 = u;  u2.frac_s_in_oil = 0.25;
geo2 = tric_geometry(u2, pc);

assert(abs(geo2.f_so - 0.25) < 1e-15, 'f_so не подхватился');
assert(geo2.phi_so0 > 0,              'твёрдое в нефти обязано появиться');
% Сумма долей подачи не меняется: твёрдое делится, а не добавляется
assert(abs((u2.eps_o + u2.eps_wd + u2.eps_wf + u2.eps_s) - 1) < 1e-12, ...
       'состав подачи разбалансирован');

[moved2, so_next2, so_out2] = tric_cascade_solids_oil(Z, geo2, u2, pc, dt);
assert(any(moved2(:) > 0), 'переток обязан быть ненулевым');
assert(so_out2 > 0,        'что-то обязано уйти с нефтью');

% Захваченное плюс оставшееся = вошедшее (первая ячейка, машина пустая)
[xs, ~, wso] = tric_psd_solids(u2, pc);
in1  = geo2.phi_so0 * wso.';
out1 = moved2(1,:) + so_next2(1,:);
fprintf('  баланс первой ячейки: ');  chk(out1, in1, 1e-14);

%% --- нефть чистится на порядок хуже воды ---
% (rho_s-rho_o)/eta_o против (rho_s-rho_w)/eta_w
k_s  = (pc.rho_s - pc.rho_w) * geo2.omega^2 / (18 * pc.eta_w);
k_so = (pc.rho_s - pc.rho_o) * geo2.omega^2 / (18 * u2.eta_o);
fprintf('  k_s / k_so = %.1f (ожидается ~28)\n', k_s / k_so);
assert(k_s / k_so > 20, 'осаждение в нефти обязано быть много медленнее');

%% --- переток попадает в водяное кольцо ---
[cap_a, ~, ~] = tric_cascade_solids(Z, Z,      v.Rtr, geo2, u2, pc, dt);
[cap_b, ~, ~] = tric_cascade_solids(Z, moved2, v.Rtr, geo2, u2, pc, dt);
assert(all(cap_b >= cap_a - 1e-18), 'переток не может уменьшить захват');
assert(sum(cap_b) > sum(cap_a),     'переток обязан увеличить захват');
fprintf('  переток увеличивает захват: ок\n');

%% --- схлопнувшееся водяное кольцо ---
% Кек дорос до границы раздела: перетекшему осаждаться негде,
% оно попадает на кек сразу, минуя каскад.
Rtr_c = geo2.r_i * ones(n, 1);
[cap_c, ~, ~] = tric_cascade_solids(Z, moved2, Rtr_c, geo2, u2, pc, dt);
q = geo2.Qo / geo2.Qw;
fprintf('  схлопнулось, ячейка 1: ');  chk(cap_c(1), sum(moved2(1,:)) * q, 1e-15);

%% --- размеры ---
assert(isequal(size(moved2),   [n J]), 'moved: ожидалось %dx%d', n, J);
assert(isequal(size(so_next2), [n J]), 'phi_so_next: ожидалось %dx%d', n, J);
assert(isscalar(so_out2),              'phi_so_out обязан быть скаляром');

fprintf('\nТвёрдое в нефти: проверки пройдены.\n');
end
