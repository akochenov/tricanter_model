"""
Эталоны для динамического варианта (simulate_dyn, §7.2 документа).

Концентрации phi_s и phi_w здесь — состояния размером n x J = 25 x 40.
Обновляются точной экспонентой, а не интегрированием.

Запуск:  python dump_reference_dyn.py
"""
import csv
import os

import numpy as np

import tricanter3 as T

REF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ref_dyn")
os.makedirs(REF, exist_ok=True)


def save_mat(name, arr):
    """Матрица n x J -> CSV без заголовка."""
    path = os.path.join(REF, name + ".csv")
    np.savetxt(path, np.atleast_2d(arr), delimiter=",", fmt="%.17g")
    print(f"  {name}.csv  {np.atleast_2d(arr).shape}")


def save_cols(name, **cols):
    keys = list(cols)
    n = max(np.atleast_1d(cols[k]).size for k in keys)
    path = os.path.join(REF, name + ".csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for i in range(n):
            w.writerow([f"{np.atleast_1d(cols[k])[i]:.17g}" for k in keys])
    print(f"  {name}.csv  ({n} строк)")


def scalars(name, **vals):
    path = os.path.join(REF, name + ".csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["name", "value"])
        for k, v in vals.items():
            w.writerow([k, f"{float(v):.17g}"])
    print(f"  {name}.csv  ({len(vals)} значений)")


def one_step(p, m, phi_s, phi_w, xs, ws, xw, ww, dt):
    """Один шаг simulate_dyn, разложенный на промежуточные величины."""
    Rtr = np.sqrt(np.maximum(
        p.Rd ** 2 - m / (p.rho_s * p.phi_sed * np.pi * p.Lax), p.r_i ** 2))
    dH = p.Rd - Rtr

    tau_s = np.pi * (Rtr ** 2 - p.r_i ** 2) * p.Lax / p.Qw
    k_s = (p.rho_s - p.rho_w) * p.omega ** 2 / (18 * p.eta_w)
    pref_s = p.Rd ** 2 / max(p.Rd ** 2 - p.r_i ** 2, 1e-12)
    phi_s0 = p.eps_s / (p.eps_w + p.eps_s)
    T_s = T.grade(xs, k_s, tau_s * p.f_clar, pref_s, phi_s, p.phi_ref)
    up_s = np.vstack([phi_s0 * ws, phi_s[:-1]])
    cap_s = (up_s * T_s).sum(1)
    star_s = up_s * (1 - T_s)
    decay_s = np.exp(-dt / np.maximum(tau_s, 1e-9))[:, None]
    phi_s_next = star_s + (phi_s - star_s) * decay_s

    if p.has_oil and p.r_i - p.Ro > 1e-6:
        tau_o = np.full(p.n, np.pi * (p.r_i ** 2 - p.Ro ** 2) * p.Lax) / p.Qo
        k_w = (p.rho_w - p.rho_o) * p.omega ** 2 / (18 * p.eta_o)
        pref_w = p.r_i ** 2 / max(p.r_i ** 2 - p.Ro ** 2, 1e-12)
        phi_w0 = p.eps_wd / (p.eps_o + p.eps_wd)
        T_w = T.grade(xw, k_w, tau_o * p.f_clar, pref_w, phi_w, p.phi_max)
        up_w = np.vstack([phi_w0 * ww, phi_w[:-1]])
        star_w = up_w * (1 - T_w)
        phi_w_next = star_w + (phi_w - star_w) * np.exp(-dt / tau_o)[:, None]
    else:
        tau_o = np.zeros(p.n)
        k_w = pref_w = phi_w0 = 0.0
        T_w = np.zeros((p.n, p.J))
        phi_w_next = np.zeros((p.n, p.J))

    tr = p.phi_sed * p.rho_s * np.pi * (p.Rd ** 2 - Rtr ** 2) * p.u_ax
    dm = p.rho_s * p.Qw * cap_s + np.append(tr[1:], 0.0) - tr

    return dict(Rtr=Rtr, dH=dH, tau_s=tau_s, tau_o=tau_o, k_s=k_s, k_w=k_w,
                pref_s=pref_s, pref_w=pref_w, phi_s0=phi_s0, phi_w0=phi_w0,
                T_s=T_s, T_w=T_w, cap_s=cap_s, tr=tr, dm=dm,
                phi_s_next=phi_s_next, phi_w_next=phi_w_next)


def dump(p, tag, dt=1.0):
    print(f"\n[{tag}]")
    xs, ws = T.psd(p.x50, p.bs, p.J, p.rrsb)
    xw, ww = T.psd(p.d50, p.bw, p.J)
    save_cols(f"psd_solids_{tag}", x=xs, w=ws)
    save_cols(f"psd_drops_{tag}", x=xw, w=ww)

    scalars(f"geometry_{tag}",
            omega=p.omega, beta=p.beta, Lax=p.Lax, L_cone=p.L_cone,
            f_clar=p.f_clar, u_ax=p.u_ax, r_i=p.r_i, Qo=p.Qo, Qw=p.Qw,
            has_oil=float(p.has_oil), dt=dt)

    # --- шаг 1: всё с нуля. Проверяет самый простой случай ---
    m = np.zeros(p.n)
    phi_s = np.zeros((p.n, p.J))
    phi_w = np.zeros((p.n, p.J))
    r = one_step(p, m, phi_s, phi_w, xs, ws, xw, ww, dt)
    save_cols(f"step1_vec_{tag}", Rtr=r["Rtr"], dH=r["dH"], tau_s=r["tau_s"],
              tau_o=r["tau_o"], cap_s=r["cap_s"], tr=r["tr"], dm=r["dm"])
    scalars(f"step1_sca_{tag}", k_s=r["k_s"], k_w=r["k_w"],
            pref_s=r["pref_s"], pref_w=r["pref_w"],
            phi_s0=r["phi_s0"], phi_w0=r["phi_w0"])
    save_mat(f"step1_Ts_{tag}", r["T_s"])
    save_mat(f"step1_phis_next_{tag}", r["phi_s_next"])
    save_mat(f"step1_phiw_next_{tag}", r["phi_w_next"])

    # --- прогон до установившегося режима ---
    m = np.zeros(p.n)
    phi_s = np.zeros((p.n, p.J))
    phi_w = np.zeros((p.n, p.J))
    nstep = int(1500 / dt)
    for _ in range(nstep):
        r = one_step(p, m, phi_s, phi_w, xs, ws, xw, ww, dt)
        m = np.maximum(m + dt * r["dm"], 0.0)
        phi_s, phi_w = r["phi_s_next"], r["phi_w_next"]

    save_mat(f"steady_phis_{tag}", phi_s)
    save_mat(f"steady_phiw_{tag}", phi_w)
    save_cols(f"steady_vec_{tag}", m=m, Rtr=r["Rtr"], dH=r["dH"],
              tau_s=r["tau_s"], cap_s=r["cap_s"], tr=r["tr"], dm=r["dm"])

    # --- шаг из установившегося состояния: главный тест ---
    r2 = one_step(p, m, phi_s, phi_w, xs, ws, xw, ww, dt)
    save_mat(f"steady_step_phis_next_{tag}", r2["phi_s_next"])
    save_cols(f"steady_step_vec_{tag}", cap_s=r2["cap_s"], dm=r2["dm"])

    U = np.mean(p.Rd ** 2 - r["Rtr"] ** 2) / max(p.Rd ** 2 - p.r_i ** 2, 1e-12)
    scalars(f"steady_out_{tag}",
            U=U,
            E_s=1 - phi_s[-1].sum() / r["phi_s0"],
            E_w=(1 - phi_w[-1].sum() / r["phi_w0"]) if r["phi_w0"] else 0.0,
            phi_s_out=phi_s[-1].sum(),
            phi_w_out=phi_w[-1].sum(),
            m_tot=m.sum())

    # --- траектория ---
    res = T.simulate_dyn(p, 1500.0, dt)
    idx = np.arange(0, res["t"].size, 50)
    save_cols(f"trajectory_{tag}", t=res["t"][idx], U=res["U"][idx],
              E_s=res["E_s"][idx], E_w=res["E_w"][idx],
              phi_w=res["phi_w"][idx])


if __name__ == "__main__":
    print("Эталоны динамического варианта в", os.path.abspath(REF))
    dump(T.Params(), "3ph")
    dump(T.two_phase(), "2ph")
    print("\nГотово.")
