function test_B_grade()
%TEST_B_GRADE Тест tric_grade. Запуск: test_B_grade
%
% Эта функция от редакции каскада не зависит: T считается по концентрациям
% состояния, а не по прогонке. Поэтому сходится с ref_dyn/step1_Ts_*.csv
% до машинной точности в обоих режимах.

fprintf('\n=== tric_grade ===\n');

for tag = {'3ph', '2ph'}
    t = tag{1};
    fprintf('\n[%s]\n', t);

    [~, pc] = ref_case(t);
    [n, J]  = tric_dims();
    g  = ref_load(['geometry_' t]);
    v  = ref_load(['step1_vec_' t]);
    s  = ref_load(['step1_sca_' t]);
    ps = ref_load(['psd_solids_' t]);

    % --- step1: машина пустая, стеснение неактивно ---
    % ВНИМАНИЕ: в CSV tau_s записан ГОЛЫЙ, а grade принимает домноженный
    % на f_clar. Забыть — потерять полдня.
    T = tric_grade(ps.x, s.k_s, v.tau_s * g.f_clar, s.pref_s, ...
                   zeros(n, J), pc.phi_ref);
    fprintf('  step1 T_s vs эталон: ');  chk(T, ref_mat(['step1_Ts_' t]));

    % --- steady: концентрации ненулевые, стеснение оживает ---
    % Прямого эталона на T нет, поэтому проверяем свойством: захват
    % обязан упасть всюду, где концентрация ненулевая.
    vs = ref_load(['steady_vec_' t]);
    P  = ref_mat(['steady_phis_' t]);
    T2 = tric_grade(ps.x, s.k_s, vs.tau_s * g.f_clar, s.pref_s, P, pc.phi_ref);
    T0 = tric_grade(ps.x, s.k_s, vs.tau_s * g.f_clar, s.pref_s, ...
                    zeros(n, J), pc.phi_ref);
    assert(all(T2(:) <= T0(:) + 1e-15), 'стеснение не снижает захват');
    drop = 1 - min(T2(T0 > 1e-12) ./ T0(T0 > 1e-12));
    fprintf('  steady стеснение живо: максимальное падение T %.3f%%\n', drop*100);
    assert(drop > 1e-6, ...
        ['множитель стеснения не написан: на steady T не изменилось. ' ...
         'На step1 такой код проходит, потому что R == 1.']);

    % --- ось суммирования ---
    % Сумма идёт ПО КЛАССАМ (в питоне sum(1), в MATLAB sum(...,2)).
    % Ставим стеснение только в одну ячейку: если ось перепутана,
    % просядут все ячейки сразу либо не просядет ни одна.
    Q = zeros(n, J);
    Q(3, :) = pc.phi_ref / (2*J);
    T3 = tric_grade(ps.x, s.k_s, vs.tau_s * g.f_clar, s.pref_s, Q, pc.phi_ref);
    hit = abs(T3 - T0) > 1e-15;
    % Насыщенные классы исключаем: там min(1,...) срезает и падение
    % захвата не проявляется, сколько ни добавляй стеснения.
    act = T0(3,:) > 1e-12 & T0(3,:) < 1 - 1e-12;
    assert(any(act), 'нет ненасыщенных классов, проверка оси невозможна');
    assert(all(hit(3, act)),              'ячейка 3 не отреагировала: ось суммы');
    assert(~any(any(hit([1:2 4:n], :))),  'отреагировали чужие ячейки: ось суммы');
    fprintf('  ось суммирования: ок (просела только ячейка 3)\n');

    % --- границы ---
    assert(all(T(:) >= 0) && all(T(:) <= 1), 'T вне [0,1]');
    assert(all(diff(T(1, :)) >= -1e-15),     'T не растёт с размером класса');
    fprintf('  границы и монотонность: ок\n');
end

fprintf('\ntric_grade: проверки пройдены.\n');
end
