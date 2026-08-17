"""
Выгрузка эталонов из питон-модели для тестов MATLAB.

Запуск:  python dump_reference.py
Результат: ref/*.csv — эталонные значения, против которых каждый
участник проверяет свою функцию.

Файл должен лежать рядом с tricanter3.py (или поправьте импорт).
"""
import csv
import os

import numpy as np

import tricanter3 as T


REF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ref")
os.makedirs(REF, exist_ok=True)


def save(name, **cols):
    """Пишет CSV с именованными колонками одинаковой длины."""
    keys = list(cols)
    ncol = len(keys)
    n = max(np.atleast_1d(cols[k]).size for k in keys)
    rows = []
    for i in range(n):
        row = []
        for k in keys:
            v = np.atleast_1d(cols[k])
            row.append(f"{v[i]:.17g}" if i < v.size else "")
        rows.append(row)
    path = os.path.join(REF, name + ".csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(keys)
        w.writerows(rows)
    print(f"  {name}.csv  ({n} x {ncol})")


def scalars(name, **vals):
    path = os.path.join(REF, name + ".csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["name", "value"])
        for k, v in vals.items():
            w.writerow([k, f"{float(v):.17g}"])
    print(f"  {name}.csv  ({len(vals)} значений)")


def dump(p, tag):
    print(f"\n[{tag}]")

    # --- 1. геометрия (tric_geometry) ---
    scalars(f"geometry_{tag}",
            omega=p.omega, beta=p.beta, Lax=p.Lax, L_cone=p.L_cone,
            f_clar=p.f_clar, u_ax=p.u_ax, r_i=p.r_i, Qo=p.Qo, Qw=p.Qw,
            has_oil=float(p.has_oil))

    # --- 2. распределения размеров (tric_psd) ---
    xs, ws = T.psd(p.x50, p.bs, p.J, p.rrsb)
    xw, ww = T.psd(p.d50, p.bw, p.J)
    save(f"psd_solids_{tag}", x=xs, w=ws)
    save(f"psd_drops_{tag}", x=xw, w=ww)

    # --- 3. геометрия кека при пустой машине и в установившемся режиме ---
    for label, m in (("empty", np.zeros(p.n)), ("steady", steady_mass(p))):
        Rtr = np.sqrt(np.maximum(
            p.Rd ** 2 - m / (p.rho_s * p.phi_sed * np.pi * p.Lax), p.r_i ** 2))
        dH = p.Rd - Rtr
        tau_s = np.pi * (Rtr ** 2 - p.r_i ** 2) * p.Lax / p.Qw
        tr = p.phi_sed * p.rho_s * np.pi * (p.Rd ** 2 - Rtr ** 2) * p.u_ax
        dm, U, E_s, E_w, phi_w, _ = T.cascade(p, m, xs, ws, xw, ww)
        save(f"cake_{label}_{tag}",
             m=m, Rtr=Rtr, dH=dH, tau_s=tau_s, tr=tr, dm=dm)

        # захват по ячейкам — вход для теста каскада
        k_s = (p.rho_s - p.rho_w) * p.omega ** 2 / (18 * p.eta_w)
        pref_s = p.Rd ** 2 / max(p.Rd ** 2 - p.r_i ** 2, 1e-12)
        phi_s0 = p.eps_s / (p.eps_w + p.eps_s)
        cap_s, phi_s_out = T.settle(phi_s0, ws, xs, k_s,
                                    tau_s * p.f_clar, pref_s, p.phi_ref)
        save(f"cascade_solids_{label}_{tag}", cap=cap_s)
        scalars(f"outputs_{label}_{tag}",
                phi_s_out=phi_s_out, phi_s0=phi_s0,
                U=U, E_s=E_s, E_w=E_w, phi_w_out=phi_w,
                k_s=k_s, pref_s=pref_s)

    # --- 4. траектория (тест схемы целиком) ---
    r = T.simulate(p, t_end=1500.0, dt=1.0)
    idx = np.arange(0, r["t"].size, 50)
    save(f"trajectory_{tag}",
         t=r["t"][idx], U=r["U"][idx], E_s=r["E_s"][idx],
         E_w=r["E_w"][idx], phi_w=r["phi_w"][idx])
    save(f"profile_final_{tag}", dH=r["dH"][-1])


def steady_mass(p, t_end=1500.0, dt=1.0):
    """Масса кека в установившемся режиме — вход для тестов."""
    xs, ws = T.psd(p.x50, p.bs, p.J, p.rrsb)
    xw, ww = T.psd(p.d50, p.bw, p.J)
    m = np.zeros(p.n)
    for _ in range(int(t_end / dt)):
        dm, *_ = T.cascade(p, m, xs, ws, xw, ww)
        m = np.maximum(m + dt * dm, 0.0)
    return m


if __name__ == "__main__":
    print("Выгрузка эталонов в", os.path.abspath(REF))
    dump(T.Params(), "3ph")
    dump(T.two_phase(), "2ph")

    # контрольная таблица §11: E_s при трёх расходах, двухфазный режим
    rows = []
    for Q_lph in (30.0, 45.0, 60.0):
        p = T.two_phase(Q=Q_lph / 3.6e6)
        r = T.simulate(p, 1500.0)
        rows.append((Q_lph, r["E_s"][-1], r["dH"][-1, 0]))
    with open(os.path.join(REF, "gleiss_table.csv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Q_lph", "E_s", "dH_inlet_m"])
        for row in rows:
            w.writerow([f"{v:.17g}" for v in row])
    print("\n  gleiss_table.csv — контрольная таблица §11")
    print("\nГотово. Закоммитьте папку ref/ в репозиторий.")
