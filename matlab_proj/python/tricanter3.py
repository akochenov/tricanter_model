"""
Трёхфазная модель трикантера (нефть / вода / твёрдое).

Обобщение двухфазной модели (gleiss.py). Кольца от оси наружу:
    Ro (слив нефти) < r_i (граница нефть-вода) < Rtr (кек) < Rd (стенка)
Нефть несёт капли воды, они оседают на границу r_i; вода несёт твёрдое,
оно оседает на кек. Состояния — только масса кека, остальное алгебра.

Положение r_i задаётся балансом давлений двух сливов (U-трубка в поле
ω²r) и от кека не зависит. При eps_o = 0 получаем r_i = Rw и модель
в точности совпадает с двухфазной.

Расходы нефти и воды считаются постоянными вдоль машины (капли разбавлены).
"""
import numpy as np
from dataclasses import dataclass

G = 9.81


@dataclass
class Params:
    Rd: float = 0.040        # радиус барабана, м
    Rw: float = 0.0287       # слив воды, м
    Ro: float = 0.028        # слив нефти, м
    Wsc: float = 0.025       # шаг шнека, м
    Lpond: float = 0.17      # длина пруда (объём, время пребывания), м
    alpha: float = 7.0       # полуугол конуса, град
    Lsep: float = 0.0        # путь осаждения; 0 = вывести из геометрии конуса
    n: int = 25
    J: int = 40

    rho_s: float = 2720.0
    rho_w: float = 998.0
    rho_o: float = 870.0
    eta_w: float = 1e-3      # вязкость воды, Па·с
    eta_o: float = 0.03      # вязкость нефти, Па·с

    C: float = 1500.0        # фактор разделения, g
    Q: float = 30 / 3.6e6    # подача, м3/с
    dn: float = 10.0         # диф. скорость шнека, об/мин

    eps_o: float = 0.60      # нефть в подаче
    eps_wd: float = 0.08     # вода, диспергированная в нефти
    eps_w: float = 0.30      # свободная вода
    eps_s: float = 0.02      # твёрдое

    mu_s: float = 0.25       # трение кека
    phi_sed: float = 0.45    # доля твёрдого в осадке
    phi_ref: float = 0.13    # опорная доля для твёрдого в воде
    phi_max: float = 0.64    # плотная упаковка капель
    u_conv: float = 0.57     # эффективность конвейирования

    x50: float = 1.5e-6      # медиана частиц, м
    bs: float = 3.0
    rrsb: bool = False       # форма распределения твёрдого
    d50: float = 3.0e-6      # медиана капель воды, м
    bw: float = 2.0

    @property
    def omega(self):
        return np.sqrt(self.C * G / self.Rd)

    @property
    def beta(self):
        return np.arctan2(self.Wsc, 2 * np.pi * self.Rd)

    @property
    def Lax(self):
        """Осевая длина ячейки."""
        return self.Lpond / self.n

    @property
    def L_cone(self):
        """Вклад конуса в осветление, в эквивалентных метрах цилиндра.
        Захват на осевой срез ~ R(x)^2, отсюда int R^2 dx / Rd^2."""
        t = np.tan(np.radians(self.alpha))
        return (self.Rd ** 3 - self.Rw ** 3) / (3 * t) / self.Rd ** 2

    @property
    def f_clar(self):
        """Во сколько раз осветление сильнее, чем даёт один цилиндр."""
        L = self.Lsep if self.Lsep > 0 else self.Lpond + self.L_cone
        return L / self.Lpond

    @property
    def u_ax(self):
        """Осевая скорость выноса кека шнеком."""
        kappa = np.arctan(self.mu_s) + self.beta
        eff = 1 / (1 + np.tan(self.beta) * np.tan(kappa))
        return self.Wsc * self.dn / 60 * eff * self.u_conv

    @property
    def has_oil(self):
        return self.eps_o + self.eps_wd > 0

    @property
    def r_i(self):
        """Граница нефть-вода из баланса давлений двух сливов."""
        if not self.has_oil:
            return self.Rw
        r2 = (self.rho_w * self.Rw ** 2 - self.rho_o * self.Ro ** 2) \
            / (self.rho_w - self.rho_o)
        return float(np.clip(np.sqrt(max(r2, 0.0)), self.Ro, self.Rd))

    @property
    def Qo(self):
        return self.Q * (self.eps_o + self.eps_wd)

    @property
    def Qw(self):
        return self.Q * (self.eps_w + self.eps_s)


def psd(d50, b, J, rrsb=False):
    """Классы дисперсной фазы и их массовые доли."""
    e = np.logspace(np.log10(d50 / 30), np.log10(d50 * 30), J + 1)
    Q3 = 1 - np.exp(-(e / d50) ** b) if rrsb else 1 / (1 + (d50 / e) ** b)
    w = np.diff(Q3)
    return np.sqrt(e[:-1] * e[1:]), w / w.sum()


def grade(x, k0, tau, pref, phi_cells, phi_ref):
    """Матрица эффективности захвата T[ячейка, класс]."""
    R = np.maximum(0.0, 1 - phi_cells.sum(1) / phi_ref) ** 4.65
    arg = 2 * k0 * x[None, :] ** 2 * (R * tau)[:, None]
    return np.minimum(1.0, pref * (1 - np.exp(-arg)))


def settle(phi0, w, x, k0, tau, pref, phi_ref):
    """Прогонка каскада для одной дисперсной фазы.
    Возвращает захват по ячейкам (в долях расхода) и остаток на выходе."""
    n = tau.size
    phi = phi0 * w
    cap = np.zeros(n)
    for i in range(n):
        R = max(0.0, 1 - phi.sum() / phi_ref) ** 4.65
        T = np.minimum(1.0, pref * (1 - np.exp(-2 * k0 * x ** 2 * R * tau[i])))
        cap[i] = (phi * T).sum()
        phi = phi * (1 - T)
    return cap, phi.sum()


def cascade(p, m, xs, ws, xw, ww):
    """По массе кека возвращает dm/dt и наблюдаемые величины."""
    Rtr = np.sqrt(np.maximum(
        p.Rd ** 2 - m / (p.rho_s * p.phi_sed * np.pi * p.Lax), p.r_i ** 2))
    dH = p.Rd - Rtr

    # твёрдое оседает в воде: кольцо r_i .. Rtr, захват на кек
    tau_s = np.pi * (Rtr ** 2 - p.r_i ** 2) * p.Lax / p.Qw
    k_s = (p.rho_s - p.rho_w) * p.omega ** 2 / (18 * p.eta_w)
    pref_s = p.Rd ** 2 / max(p.Rd ** 2 - p.r_i ** 2, 1e-12)
    phi_s0 = p.eps_s / (p.eps_w + p.eps_s)
    cap_s, phi_s_out = settle(phi_s0, ws, xs, k_s, tau_s * p.f_clar, pref_s, p.phi_ref)

    # капли воды оседают в нефти: кольцо Ro .. r_i, захват на границу
    if p.has_oil and p.r_i - p.Ro > 1e-6:
        tau_o = np.full(p.n, np.pi * (p.r_i ** 2 - p.Ro ** 2) * p.Lax) / p.Qo
        k_w = (p.rho_w - p.rho_o) * p.omega ** 2 / (18 * p.eta_o)
        pref_w = p.r_i ** 2 / max(p.r_i ** 2 - p.Ro ** 2, 1e-12)
        phi_w0 = p.eps_wd / (p.eps_o + p.eps_wd)
        _, phi_w_out = settle(phi_w0, ww, xw, k_w, tau_o * p.f_clar, pref_w, p.phi_max)
    elif p.has_oil:
        phi_w0 = p.eps_wd / (p.eps_o + p.eps_wd)
        phi_w_out = phi_w0                      # нефтяного кольца нет
    else:
        phi_w0, phi_w_out = 0.0, 0.0

    tr = p.phi_sed * p.rho_s * np.pi * (p.Rd ** 2 - Rtr ** 2) * p.u_ax
    dm = p.rho_s * p.Qw * cap_s + np.append(tr[1:], 0.0) - tr

    U = np.mean(p.Rd ** 2 - Rtr ** 2) / max(p.Rd ** 2 - p.r_i ** 2, 1e-12)
    E_s = 1 - phi_s_out / phi_s0 if phi_s0 else 0.0
    E_w = 1 - phi_w_out / phi_w0 if phi_w0 else 0.0
    return dm, U, E_s, E_w, phi_w_out, dH


def simulate_dyn(p, t_end, dt, inputs=None):
    """Концентрации в состояниях: даёт запаздывание, потерянное в квазистатике."""
    t = np.arange(0, t_end + dt, dt)
    live_psd = bool(inputs) and bool(
        {"x50", "bs", "rrsb", "d50", "bw", "J"} & inputs.keys())
    xs, ws = psd(p.x50, p.bs, p.J, p.rrsb)
    xw, ww = psd(p.d50, p.bw, p.J)

    m = np.zeros(p.n)
    phi_s = np.zeros((p.n, p.J))
    phi_w = np.zeros((p.n, p.J))
    U, E_s, E_w, w_oil = (np.zeros(t.size) for _ in range(4))
    dH_hist = np.zeros((t.size, p.n))

    for i in range(t.size):
        if inputs:
            for name, fn in inputs.items():
                setattr(p, name, fn(t[i]))
            if live_psd:
                xs, ws = psd(p.x50, p.bs, p.J, p.rrsb)
                xw, ww = psd(p.d50, p.bw, p.J)

        Rtr = np.sqrt(np.maximum(
            p.Rd ** 2 - m / (p.rho_s * p.phi_sed * np.pi * p.Lax), p.r_i ** 2))
        dH = p.Rd - Rtr

        # твёрдое в воде
        tau_s = np.pi * (Rtr ** 2 - p.r_i ** 2) * p.Lax / p.Qw
        k_s = (p.rho_s - p.rho_w) * p.omega ** 2 / (18 * p.eta_w)
        pref_s = p.Rd ** 2 / max(p.Rd ** 2 - p.r_i ** 2, 1e-12)
        phi_s0 = p.eps_s / (p.eps_w + p.eps_s)
        T_s = grade(xs, k_s, tau_s * p.f_clar, pref_s, phi_s, p.phi_ref)
        decay_s = np.exp(-dt / np.maximum(tau_s, 1e-9))
        cap_s = np.zeros(p.n)
        phi_s_next = np.empty_like(phi_s)
        up = phi_s0 * ws                    # вход первой ячейки
        for c in range(p.n):               # c, не i: i занят внешним циклом по времени
            cap_s[c] = (up * T_s[c]).sum()
            star_c = up * (1 - T_s[c])
            phi_s_next[c] = star_c + (phi_s[c] - star_c) * decay_s[c]
            up = phi_s_next[c]              # уже на текущем шаге, не на прошлом
        phi_s = phi_s_next

        # капли воды в нефти
        if p.has_oil and p.r_i - p.Ro > 1e-6:
            tau_o = np.full(p.n, np.pi * (p.r_i ** 2 - p.Ro ** 2) * p.Lax) / p.Qo
            k_w = (p.rho_w - p.rho_o) * p.omega ** 2 / (18 * p.eta_o)
            pref_w = p.r_i ** 2 / max(p.r_i ** 2 - p.Ro ** 2, 1e-12)
            phi_w0 = p.eps_wd / (p.eps_o + p.eps_wd)
            T_w = grade(xw, k_w, tau_o * p.f_clar, pref_w, phi_w, p.phi_max)
            decay_w = np.exp(-dt / tau_o)
            phi_w_next = np.empty_like(phi_w)
            up = phi_w0 * ww
            for c in range(p.n):
                star_c = up * (1 - T_w[c])
                phi_w_next[c] = star_c + (phi_w[c] - star_c) * decay_w[c]
                up = phi_w_next[c]
            phi_w = phi_w_next
            out_w = phi_w[-1].sum()
        else:
            phi_w0 = p.eps_wd / (p.eps_o + p.eps_wd) if p.has_oil else 0.0
            out_w = phi_w0

        tr = p.phi_sed * p.rho_s * np.pi * (p.Rd ** 2 - Rtr ** 2) * p.u_ax
        m = np.maximum(m + dt * (p.rho_s * p.Qw * cap_s
                                 + np.append(tr[1:], 0.0) - tr), 0.0)

        U[i] = np.mean(p.Rd ** 2 - Rtr ** 2) / max(p.Rd ** 2 - p.r_i ** 2, 1e-12)
        E_s[i] = 1 - phi_s[-1].sum() / phi_s0 if phi_s0 else 0.0
        E_w[i] = 1 - out_w / phi_w0 if phi_w0 else 0.0
        w_oil[i] = out_w
        dH_hist[i] = dH
    return dict(t=t, U=U, E_s=E_s, E_w=E_w, phi_w=w_oil, dH=dH_hist)


def step(t0, before, after):
    """Ступенька: значение before до t0, after после."""
    return lambda t: before if t < t0 else after


def simulate(p, t_end=1500.0, dt=1.0, inputs=None, dynamic=False):
    """Явный Эйлер по массе кека.

    inputs — расписание входов, {имя параметра: функция от времени}, напр.
        inputs={"Q": step(300, 30/3.6e6, 45/3.6e6)}
    dynamic — концентрации тоже состояния (запаздывание по времени
        пребывания жидкости). Обновляются точной экспонентой, поэтому
        шаг не ограничен; при dt >> tau ветвь сходится к квазистатике.
    """
    if dynamic:
        return simulate_dyn(p, t_end, dt, inputs)
    t = np.arange(0, t_end + dt, dt)
    # распределения пересчитываются на каждом шаге, только если их меняют
    live_psd = bool(inputs) and bool(
        {"x50", "bs", "rrsb", "d50", "bw", "J"} & inputs.keys())
    xs, ws = psd(p.x50, p.bs, p.J, p.rrsb)
    xw, ww = psd(p.d50, p.bw, p.J)
    m = np.zeros(p.n)
    U, E_s, E_w, phi_w = (np.zeros(t.size) for _ in range(4))
    dH = np.zeros((t.size, p.n))
    for i in range(t.size):
        if inputs:
            for name, fn in inputs.items():
                setattr(p, name, fn(t[i]))
            if live_psd:
                xs, ws = psd(p.x50, p.bs, p.J, p.rrsb)
                xw, ww = psd(p.d50, p.bw, p.J)
        dm, U[i], E_s[i], E_w[i], phi_w[i], dH[i] = cascade(p, m, xs, ws, xw, ww)
        m = np.maximum(m + dt * dm, 0.0)
    return dict(t=t, U=U, E_s=E_s, E_w=E_w, phi_w=phi_w, dH=dH)


def weir_for_interface(p, r_i):
    """Радиус водяного слива, дающий заданную границу нефть-вода."""
    return np.sqrt((r_i ** 2 * (p.rho_w - p.rho_o) + p.rho_o * p.Ro ** 2) / p.rho_w)


def check_limits(p):
    """Проверка, что режим физически существует."""
    msg = []
    if not p.has_oil:
        return msg
    if p.r_i <= p.Ro + 1e-6:
        msg.append("нефтяное кольцо схлопнулось: вода уходит в нефтяной слив")
    elif p.r_i - p.Ro < 1e-3:
        msg.append("нефтяное кольцо тоньше 1 мм: режим на грани")
    if p.r_i >= p.Rd - 1e-6:
        msg.append("водяное кольцо схлопнулось: нефть уходит в водяной слив")
    return msg


def two_phase(**kw):
    """Двухфазный частный случай: нефти нет, r_i = Rw."""
    p = Params(eps_o=0.0, eps_wd=0.0, eps_w=0.98, eps_s=0.02,
               rho_s=1410, C=250, dn=5, mu_s=0.3, phi_sed=0.5,
               phi_ref=0.5, x50=2.5e-6, bs=2.78, rrsb=True, Rw=0.034)
    for k, v in kw.items():
        setattr(p, k, v)
    return p


def plot(p=None, t_end=1500.0):
    import matplotlib.pyplot as plt
    p = p or Params()
    r = simulate(p, t_end)
    s = np.arange(p.n) * p.Lax

    fig, ax = plt.subplots(1, 4, figsize=(16, 3.6))
    ax[0].plot(r['t'], r['U']);                ax[0].set(xlabel="t, с", ylabel="U")
    ax[1].plot(r['t'], r['E_s'], color="g");   ax[1].set(xlabel="t, с", ylabel="E твёрдого")
    ax[2].plot(r['t'], r['phi_w'] * 100, color="b")
    ax[2].set(xlabel="t, с", ylabel="вода в нефти, %об")
    for i in range(0, r['t'].size, max(1, r['t'].size // 12)):
        ax[3].plot(s, r['dH'][i] * 1000, color=plt.cm.viridis(i / r['t'].size), lw=1)
    ax[3].set(xlabel="осевая координата, м", ylabel="высота осадка, мм")
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    p = Params()
    print(f"r_i = {p.r_i*1000:.1f} мм   (нефть {p.Ro*1000:.0f}..{p.r_i*1000:.1f}, "
          f"вода {p.r_i*1000:.1f}..{p.Rd*1000:.0f})")
    for w in check_limits(p):
        print("!", w)

    r = simulate(p, 900)
    print(f"\nтрёхфаза:  U={r['U'][-1]:.3f}  E твёрдого={r['E_s'][-1]:.3f}  "
          f"вода в нефти {r['phi_w'][-1]*100:.2f}% (было {p.eps_wd/(p.eps_o+p.eps_wd)*100:.1f}%)")

    r2 = simulate(two_phase(), 1500)
    print(f"двухфаза:  U={r2['U'][-1]:.3f}  E={r2['E_s'][-1]:.3f}  "
          f"осадок у входа {r2['dH'][-1, 0]*1000:.1f} мм")
