"""
Двухфазная модель декантера (Gleiss, CES-2017 / CET-2018).

Состояния — только масса кека по ячейкам; концентрации считаются
алгебраически (время пребывания жидкости << времени жизни кека).

Отличия от печатных уравнений статьи:
  1. точная кольцевая геометрия (не развёртка канала);
  2. префактор в T(x) берётся по радиусу барабана, а не по радиусу кека;
  3. T(x) в точной форме, у Gleiss она линеаризована;
  4. путь осаждения 0.21 м (конус) при длине пруда 0.17 м;
  5. u_conv — эффективность конвейирования кека, калибруется.
"""
import numpy as np
from dataclasses import dataclass

G = 9.81


@dataclass
class Params:
    Rd: float = 0.040        # радиус барабана, м
    Rw: float = 0.034        # радиус слива, м
    Wsc: float = 0.025       # шаг шнека, м
    Lpond: float = 0.17      # длина пруда (объём, время пребывания), м
    alpha: float = 7.0       # полуугол конуса, град
    Lsep: float = 0.0        # путь осаждения; 0 = вывести из геометрии конуса
    n: int = 25              # ячеек
    J: int = 40              # классов частиц

    rho_s: float = 1410.0    # плотность твёрдого, кг/м3
    rho_l: float = 998.0
    eta: float = 1e-3        # вязкость жидкости, Па·с

    C: float = 250.0         # фактор разделения, g
    Q: float = 30 / 3.6e6    # подача, м3/с
    phi_in: float = 0.02     # доля твёрдого в подаче
    dn: float = 5.0          # диф. скорость шнека, об/мин

    mu_s: float = 0.3        # трение кека
    phi_sed: float = 0.5     # доля твёрдого в осадке
    phi_ref: float = 0.5     # опорная доля в R(phi)
    u_conv: float = 0.57     # эффективность конвейирования

    x50: float = 2.5e-6      # медиана распределения, м
    b: float = 2.78          # ширина распределения
    rrsb: bool = True        # True: RRSB, False: сигмоида

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


def psd(p):
    """Классы частиц и их массовые доли."""
    e = np.logspace(np.log10(p.x50 / 30), np.log10(p.x50 * 30), p.J + 1)
    Q3 = 1 - np.exp(-(e / p.x50) ** p.b) if p.rrsb else 1 / (1 + (p.x50 / e) ** p.b)
    w = np.diff(Q3)
    return np.sqrt(e[:-1] * e[1:]), w / w.sum()


def cascade(p, m, x, w):
    """По массе кека m возвращает dm/dt, заполнение U, эффективность E,
    высоту осадка dH."""
    Rtr = np.sqrt(np.maximum(
        p.Rd ** 2 - m / (p.rho_s * p.phi_sed * np.pi * p.Lax), p.Rw ** 2))
    dH = p.Rd - Rtr
    tau = np.pi * (Rtr ** 2 - p.Rw ** 2) * p.Lax / p.Q
    pref = p.Rd ** 2 / (p.Rd ** 2 - p.Rw ** 2)

    phi = p.phi_in * w
    sep = np.zeros(p.n)
    for i in range(p.n):
        R = max(0.0, 1 - phi.sum() / p.phi_ref) ** 4.65
        k = (p.rho_s - p.rho_l) * p.omega ** 2 * x ** 2 * R / (18 * p.eta)
        T = np.minimum(1.0, pref * (1 - np.exp(-2 * k * tau[i] * p.f_clar)))
        sep[i] = p.rho_s * p.Q * (phi * T).sum()
        phi = phi * (1 - T)

    tr = p.phi_sed * p.rho_s * np.pi * (p.Rd ** 2 - Rtr ** 2) * p.u_ax
    dm = sep + np.append(tr[1:], 0.0) - tr
    U = np.mean(p.Rd ** 2 - Rtr ** 2) / (p.Rd ** 2 - p.Rw ** 2)
    return dm, U, 1 - phi.sum() / p.phi_in, dH


def simulate(p, t_end=1500.0, dt=1.0, Q_of_t=None):
    """Явный Эйлер по массе кека. Q_of_t(t) — расход, если он меняется."""
    x, w = psd(p)
    t = np.arange(0, t_end + dt, dt)
    m = np.zeros(p.n)
    U, E = np.zeros(t.size), np.zeros(t.size)
    dH = np.zeros((t.size, p.n))
    for i, ti in enumerate(t):
        if Q_of_t is not None:
            p.Q = Q_of_t(ti)
        dm, U[i], E[i], dH[i] = cascade(p, m, x, w)
        m = np.maximum(m + dt * dm, 0.0)
    return t, U, E, dH


def limestone(**kw):
    """Пресет известняка (CET-2018) вместо PVC по умолчанию."""
    p = Params(rho_s=2720, C=1500, dn=10, mu_s=0.25, phi_sed=0.45,
               phi_ref=0.13, x50=1.5e-6, b=3.0, rrsb=False)
    for k, v in kw.items():
        setattr(p, k, v)
    return p


def validate():
    """Сверка с экспериментом Gleiss (CES-2017, PVC)."""
    print("E_sep против эксперимента:")
    for Q, ref in [(30, 0.72), (45, 0.60), (60, 0.51)]:
        _, _, E, _ = simulate(Params(Q=Q / 3.6e6))
        print(f"  {Q} л/ч:  модель {E[-1]:.2f}   опыт {ref:.2f}")

    print("\nвысота осадка у входа:")
    for phi, ref in [(0.02, 0.9), (0.05, 1.9), (0.08, 2.5)]:
        _, _, _, dH = simulate(Params(phi_in=phi))
        print(f"  phi={phi}:  модель {dH[-1, 0] * 1000:.1f} мм   статья {ref} мм")


def plot(p=None, t_end=1500.0):
    """U(t), E(t) и профиль осадка."""
    import matplotlib.pyplot as plt
    p = p or Params()
    t, U, E, dH = simulate(p, t_end)
    s = np.arange(p.n) * p.Lax

    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    ax[0].plot(t, U);                 ax[0].set(xlabel="t, с", ylabel="U")
    ax[1].plot(t, E, color="g");      ax[1].set(xlabel="t, с", ylabel="E_sep")
    for i in range(0, t.size, max(1, t.size // 12)):
        ax[2].plot(s, dH[i] * 1000, color=plt.cm.viridis(i / t.size), lw=1)
    ax[2].set(xlabel="осевая координата, м", ylabel="высота осадка, мм")
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    validate()
