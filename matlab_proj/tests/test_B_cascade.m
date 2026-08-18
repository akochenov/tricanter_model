function test_B_cascade()
%TEST_B_CASCADE Тест tric_cascade_solids и tric_cascade_drops.
%   Запуск: test_B_cascade
%
% Единственная чужая зависимость — tric_geometry (готова). Rtr берётся
% из CSV, а не из tric_cake_geom: модуль C проверяется отдельно.
%
% Эталоны перевыгружены под прогоночную редакцию каскада, сверка прямая.
% См. REDACTION_B.md.

fprintf('\n=== tric_cascade_solids / tric_cascade_drops ===\n');

for tag = {'3ph', '2ph'}
    t = tag{1};
    fprintf('\n[%s]\n', t);

    [u, pc, dt] = ref_case(t);
    [n, J] = tric_dims();
    geo = tric_geometry(u, pc);

    v   = ref_load(['step1_vec_' t]);
    so  = ref_load(['steady_out_' t]);
    vs  = ref_load(['steady_vec_' t]);
    ssv = ref_load(['steady_step_vec_' t]);

    phis1_ref = ref_mat(['step1_phis_next_' t]);
    phiw1_ref = ref_mat(['step1_phiw_next_' t]);
    Psteady   = ref_mat(['steady_phis_' t]);
    Wsteady   = ref_mat(['steady_phiw_' t]);

    %% ---------------- твёрдое ----------------
    fprintf('  -- solids --\n');

    % step1 держит все ловушки твёрдого: k_s через eta_w,
    % pref_s через Rd (не Rtr), phi_s0 через eps_wf, tau с f_clar в grade.
    [cap_s, phi_s_next, ~] = tric_cascade_solids(zeros(n, J), v.Rtr, geo, u, pc, dt);
    fprintf('    step1 cap_s vs эталон:        ');  chk(cap_s, v.cap_s);
    fprintf('    step1 phi_s_next vs эталон:   ');  chk(phi_s_next, phis1_ref);

    % steady. Здесь оживает множитель стеснения.
    [cap_s2, phi_s2, phi_s_out] = tric_cascade_solids(Psteady, vs.Rtr, geo, u, pc, dt);
    fprintf('    steady cap_s vs эталон:       ');  chk(cap_s2, ssv.cap_s);
    fprintf('    steady phi_s_next vs эталон:  ');
    chk(phi_s2, ref_mat(['steady_step_phis_next_' t]));
    fprintf('    steady phi_s_out vs эталон:   ');  chk(phi_s_out, so.phi_s_out);

    % Установившийся режим обязан быть неподвижной точкой каскада.
    fprintf('    steady — неподвижная точка:   ');  chk(phi_s2, Psteady, 1e-9);

    %% ---------------- капли ----------------
    fprintf('  -- drops --\n');

    [phi_w_next, phi_w_out] = tric_cascade_drops(zeros(n, J), geo, u, pc, dt);
    % Ловушки капель: eta_o из ВХОДОВ (не eta_w), pref_w через r_i
    % (не через Rd), опорная концентрация phi_max (не phi_ref).
    fprintf('    step1 phi_w_next vs эталон:   ');  chk(phi_w_next, phiw1_ref);

    [phi_w2, phi_w_out2] = tric_cascade_drops(Wsteady, geo, u, pc, dt);
    fprintf('    steady — неподвижная точка:   ');  chk(phi_w2, Wsteady, 1e-9);
    fprintf('    steady phi_w_out vs эталон:   ');  chk(phi_w_out2, so.phi_w_out);

    %% ---------------- ветвление по режиму ----------------
    if geo.has_oil
        % Кольца нет. Поднимаем нефтяной слив до границы раздела: баланс
        % давлений даёт отрицательный r^2, sqrt срезается нулём, и клип
        % снизу сажает r_i ровно на Ro. Состояние не обновляется вообще,
        % на выход идёт входная концентрация.
        u3 = u;  u3.Ro = geo.r_i - 1e-9;
        geo3 = tric_geometry(u3, pc);
        [pw3, po3] = tric_cascade_drops(Wsteady, geo3, u3, pc, dt);
        assert(isequal(pw3, Wsteady), 'без кольца состояние обязано остаться прежним');
        assert(abs(po3 - u3.eps_wd/(u3.eps_o + u3.eps_wd)) < 1e-15, ...
               'без кольца на выход идёт входная концентрация');
        fprintf('    ветка «кольца нет»: ок\n');
    else
        assert(all(phi_w_next(:) == 0) && phi_w_out == 0, ...
               'двухфазный режим: капель быть не должно');
        fprintf('    ветка «нефти нет»: ок (всё по нулям)\n');
    end

    %% ---------------- размеры ----------------
    assert(isequal(size(cap_s),      [n 1]), 'cap_s: ожидалось %dx1', n);
    assert(isequal(size(phi_s_next), [n J]), 'phi_s_next: ожидалось %dx%d', n, J);
    assert(isequal(size(phi_w_next), [n J]), 'phi_w_next: ожидалось %dx%d', n, J);
    assert(isscalar(phi_s_out) && isscalar(phi_w_out), 'выходы обязаны быть скалярами');
end

fprintf('\nКаскады: проверки пройдены.\n');
end
