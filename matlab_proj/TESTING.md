# Как тестировать свой кусок

Каждая функция проверяется **в вакууме**: входы читаются числами из CSV,
выход сравнивается с CSV. Готовность чужих файлов не имеет значения.

Эталоны выгружены из питона (`python/dump_reference_dyn.py`) с точностью
`%.17g`. Критерий сходимости — `1e-15`. Не «примерно», а до последнего бита.
Расхождение `1e-6` — это ошибка в формуле, а не округление.

## Два уровня, оба обязательны

**step1** — машина пустая: `m = 0`, `phi_s = 0`, `phi_w = 0`.
Ловит опечатки. Но множитель стеснения `(1 - sum(phi)/phi_ref)^4.65`
при нулевых концентрациях равен единице — то есть **функцию можно
написать без него, и step1 всё равно пройдёт**.

**steady** — установившийся режим, концентрации ненулевые.
Ловит недописанную физику.

Пройти надо оба. Если сошёлся только step1 — работа не сделана.

## Файлы эталонов

Суффикс `_3ph` — трёхфазный режим (основной), `_2ph` — двухфазный
(для контрольных значений §11). Тестировать на обоих.

| файл | что внутри |
|---|---|
| `geometry_3ph.csv` | omega, beta, Lax, L_cone, f_clar, u_ax, r_i, Qo, Qw, has_oil, dt |
| `psd_solids_3ph.csv` | колонки x, w — 40 строк, твёрдое |
| `psd_drops_3ph.csv` | колонки x, w — 40 строк, капли |
| `step1_vec_3ph.csv` | Rtr, dH, tau_s, tau_o, cap_s, tr, dm — по 25 строк |
| `step1_sca_3ph.csv` | k_s, k_w, pref_s, pref_w, phi_s0, phi_w0 |
| `step1_Ts_3ph.csv` | матрица 25×40, результат grade |
| `step1_phis_next_3ph.csv` | матрица 25×40 |
| `step1_phiw_next_3ph.csv` | матрица 25×40 |
| `steady_phis_3ph.csv` | матрица 25×40, установившиеся концентрации |
| `steady_phiw_3ph.csv` | матрица 25×40 |
| `steady_vec_3ph.csv` | m, Rtr, dH, tau_s, cap_s, tr, dm — установившиеся |
| `steady_step_phis_next_3ph.csv` | шаг из установившегося состояния |
| `steady_step_vec_3ph.csv` | cap_s, dm — шаг из установившегося |
| `steady_out_3ph.csv` | U, E_s, E_w, phi_s_out, phi_w_out, m_tot |
| `trajectory_3ph.csv` | t, U, E_s, E_w, phi_w — каждые 50 шагов |

## Загрузчики

```matlab
r = ref_load('step1_vec_3ph');   % колонки -> поля структуры
M = ref_mat('step1_Ts_3ph');     % матрица целиком
```

---

# A — каркас

`tric_geometry(u, pc)` → `geo`

```matlab
g   = ref_load('geometry_3ph');
geo = tric_geometry(TricInputs0, TricParams);
chk(geo.omega, g.omega); chk(geo.r_i, g.r_i);
chk(geo.Qw, g.Qw); chk(geo.Qo, g.Qo);
chk(geo.f_clar, g.f_clar); chk(geo.u_ax, g.u_ax);
```

Отдельная проверка сведения: при `eps_o = eps_wd = 0` должно быть
`r_i == Rw` до машинной точности. Если нет — ошибка в балансе давлений.

`tric_rhs` тестируется последним, когда готовы все остальные:
один вызов из нулевого состояния должен дать `dm` из `step1_vec`,
`phi_s_next` из `step1_phis_next`.

---

# B — дисперсная фаза

## tric_psd(d50, b, rrsb) → [x, w]

```matlab
p = ref_load('psd_solids_3ph');
[x, w] = tric_psd(TricInputs0.x50, TricParams.bs, TricParams.rrsb);
chk(x, p.x); chk(w, p.w);

p = ref_load('psd_drops_3ph');
[x, w] = tric_psd(TricInputs0.d50, TricParams.bw, false);
chk(x, p.x); chk(w, p.w);
```

Зависимостей нет вообще. Начинать с этого.

## tric_grade(x, k0, tau, pref, phi_cells, phi_ref) → T_mat

```matlab
g = ref_load('geometry_3ph');
v = ref_load('step1_vec_3ph');
s = ref_load('step1_sca_3ph');
p = ref_load('psd_solids_3ph');

% step1
T = tric_grade(p.x, s.k_s, v.tau_s * g.f_clar, s.pref_s, ...
               zeros(25,40), 0.13);
chk(T, ref_mat('step1_Ts_3ph'));

% steady — здесь оживает стеснение
vs = ref_load('steady_vec_3ph');
T2 = tric_grade(p.x, s.k_s, vs.tau_s * g.f_clar, s.pref_s, ...
                ref_mat('steady_phis_3ph'), 0.13);
```

**В CSV записан голый `tau_s`, а `grade` вызывается с уже
домноженным на `f_clar`.** Забыть — потерять полдня.

`phi_ref = 0.13` для твёрдого, `phi_max` для капель.

## tric_settle_dyn(phi0, w, T_mat, phi_prev, tau, dt) → [cap, phi_next]

```matlab
[cap, phi_next] = tric_settle_dyn(s.phi_s0, p.w, ref_mat('step1_Ts_3ph'), ...
                                  zeros(25,40), v.tau_s, 1.0);
chk(cap, v.cap_s);
chk(phi_next, ref_mat('step1_phis_next_3ph'));
```

Здесь `tau` без `f_clar` — экспонента затухания идёт по реальному
времени пребывания.

## tric_cascade_solids / tric_cascade_drops

Обёртки: собирают `tau`, `k0`, `pref`, зовут `grade` и `settle_dyn`.

```matlab
geo = tric_geometry(TricInputs0, TricParams);
[cap_s, phi_s_next, ~] = tric_cascade_solids(zeros(25,40), v.Rtr, ...
                                             geo, TricInputs0, TricParams, 1.0);
chk(cap_s, v.cap_s);
chk(phi_s_next, ref_mat('step1_phis_next_3ph'));
```

Единственное место с чужой зависимостью — `tric_geometry`, она готова.
`Rtr` берётся из CSV, а не из `tric_cake_geom`.

**Критично:** в префакторе `pref_s` стоит `Rd`, не `Rtr`.
При подстановке `Rtr` модель теряет самоограничение и `U` уходит в 1.

Для капель — `steady_phiw`, `step1_phiw_next`, и несущая вязкость
`eta_o`, не `eta_w`. Перепутать — ошибка в 30 раз.

---

# C — кек

## tric_cake_geom(m, geo, pc) → [Rtr, dH]

```matlab
v = ref_load('step1_vec_3ph');
[Rtr, dH] = tric_cake_geom(zeros(25,1), geo, TricParams);
chk(Rtr, v.Rtr); chk(dH, v.dH);

vs = ref_load('steady_vec_3ph');
[Rtr, dH] = tric_cake_geom(vs.m, geo, TricParams);
chk(Rtr, vs.Rtr); chk(dH, vs.dH);
```

Второй тест важнее: при `m = 0` формула вырождается в `Rtr = Rd`,
и ошибку в объёме кольца не видно.

Нижняя отсечка по `r_i` внутри `sqrt(max(..., r_i^2))` — обязательна.

## tric_cake_balance(Rtr, cap_s, geo, pc) → [dm, tr]

```matlab
[dm, tr] = tric_cake_balance(v.Rtr, v.cap_s, geo, TricParams);
chk(dm, v.dm); chk(tr, v.tr);

vs = ref_load('steady_vec_3ph');
[dm, tr] = tric_cake_balance(vs.Rtr, vs.cap_s, geo, TricParams);
chk(tr, vs.tr);
chk(dm, vs.dm);         % ~0, установившийся режим
```

В установившемся режиме `dm` близко к нулю — приход равен уносу.
Это хорошая проверка знака: если перепутать направление сдвига
в `[tr(2:end); 0] - tr`, установившегося режима не будет вообще.

---

# D — выходы, пределы, верификация

## tric_outputs — 12 аргументов

Собирается последним. Проверка по `steady_out_3ph.csv`:

```matlab
o = ref_load('steady_out_3ph');
% U = 0.0057, E_s = 0.9978, phi_w_out = 0.0938, m_tot = 0.00188
```

## tric_limits(geo, u, pc) → flags

Проверяется не эталоном, а провокацией: подать заведомо плохие входы
и убедиться, что флаг встал. Границы — §12 постановки.

## Верификация (не тест, отдельный слой)

Сверка с питоном доказывает, что перенос верный.
Она **не** доказывает, что модель верна физически.

По §11, двухфазный режим:

| Q, л/ч | модель | эксперимент |
|---|---|---|
| 30 | 0.708 | 0.72 |
| 45 | 0.583 | — |
| 60 | 0.491 | — |

Это второй слой проверки, поверх сверки с питоном. На защите
спрашивают оба отдельно.

## trajectory_3ph.csv

Динамика: `U` и `E_s` во времени. Сравнивать не поточечно, а
по форме и времени выхода на полку. Расхождение в третьем знаке
на переходном участке допустимо — накопление ошибки за 1500 шагов.

---

# Хелпер chk

Положить в `tests/chk.m`:

```matlab
function chk(got, want, tol)
if nargin < 3, tol = 1e-12; end
e = max(abs(got(:) - want(:)));
st = 'ПРОВАЛ'; if e < tol, st = 'ок'; end
fprintf('  %-6s  max|d| = %.3e\n', st, e);
end
```

# Правило

**Параметры не трогать до недели 3.** Эталоны привязаны к
дефолтному `Params()`. Поменяли `u_conv` или `phi_sed` — все CSV
недействительны, нужен перезапуск `dump_reference_dyn.py`
и коммит новых CSV целиком.
