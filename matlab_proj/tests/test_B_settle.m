function test_B_settle()
%TEST_B_SETTLE Тест tric_settle_dyn. Запуск: test_B_settle
%
% Реализована прогоночная редакция каскада: вход ячейки берётся с текущего
% шага, up = phi_next(c-1,:) + inj(c,:). Эталоны в ref_dyn/ перевыгружены
% под неё, сверка идёт напрямую по CSV на обоих уровнях.
%
% Каскад теперь общий на все три дисперсные популяции, отсюда аргумент
% inj (переток из нефти) и второй выход cap_mat (захват по классам).
% Здесь inj нулевой: сверка идёт при frac_s_in_oil = 0.
%
% Если кто-то вернёт старую векторную редакцию (up = phi_prev(c-1,:)),
% уровень step1 развалится: ячейки 2..25 разойдутся на ~5e-04.
% Уровень steady к редакции нечувствителен и такую подмену НЕ ловит.
% См. REDACTION_B.md.

fprintf('\n=== tric_settle_dyn ===\n');

for tag = {'3ph', '2ph'}
    t = tag{1};
    fprintf('\n[%s]\n', t);

    [~, pc, dt] = ref_case(t);
    [n, J] = tric_dims();
    v  = ref_load(['step1_vec_' t]);
    s  = ref_load(['step1_sca_' t]);
    ps = ref_load(['psd_solids_' t]);

    % MATLAB не разрешает индексировать результат вызова, поэтому эталонные
    % матрицы читаются в переменные заранее.
    Ts1_mat   = ref_mat(['step1_Ts_' t]);
    phis1_ref = ref_mat(['step1_phis_next_' t]);

    %% --- step1: машина пустая ---
    % Ловит опечатки. Множитель стеснения здесь равен единице, поэтому
    % его отсутствие этим уровнем НЕ ловится — на то есть steady.
    % inj — приход из нефтяного кольца. При frac_s_in_oil = 0 его нет.
    [cap, ~, phi_next] = tric_settle_dyn(s.phi_s0, ps.w, zeros(n, J), ...
                                         Ts1_mat, zeros(n, J), v.tau_s, dt);

    fprintf('  step1 cap vs эталон:         ');  chk(cap, v.cap_s);
    fprintf('  step1 phi_next vs эталон:    ');  chk(phi_next, phis1_ref);

    % Прогонка обязана дотянуться до последней ячейки: в векторной
    % редакции при пустой машине cap(2:end) строго нулевой.
    assert(any(cap(2:end) > 0), ...
        ['cap нулевой начиная со второй ячейки — похоже, написана ' ...
         'векторная редакция вместо прогоночной.']);

    %% --- steady: концентрации ненулевые, оживает стеснение ---
    g   = ref_load(['geometry_' t]);
    vs  = ref_load(['steady_vec_' t]);
    ssv = ref_load(['steady_step_vec_' t]);
    P   = ref_mat(['steady_phis_' t]);

    % В grade tau идёт домноженным на f_clar, в settle_dyn — голым.
    % В CSV записан голый.
    T = tric_grade(ps.x, s.k_s, vs.tau_s * g.f_clar, s.pref_s, P, pc.phi_ref);
    [cap2, ~, phi2] = tric_settle_dyn(s.phi_s0, ps.w, zeros(n, J), ...
                                      T, P, vs.tau_s, dt);

    fprintf('  steady cap vs эталон:        ');  chk(cap2, ssv.cap_s);
    fprintf('  steady phi_next vs эталон:   ');
    chk(phi2, ref_mat(['steady_step_phis_next_' t]));

    %% --- вырождение в квазистатику ---
    % При dt >> tau экспонента гаснет и остаётся phi_next = star.
    [~, ~, phiQ] = tric_settle_dyn(s.phi_s0, ps.w, zeros(n, J), ...
                                   Ts1_mat, zeros(n, J), v.tau_s, 1e9);
    up   = s.phi_s0 * ps.w.';
    star = up .* (1 - Ts1_mat(1,:));
    fprintf('  dt -> inf, ячейка 1 = star:  ');  chk(phiQ(1,:), star);
end

fprintf('\ntric_settle_dyn: проверки пройдены.\n');
end
