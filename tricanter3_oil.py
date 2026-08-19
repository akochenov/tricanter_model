"""
Трёхфазная модель трикантера (нефть / вода / твёрдое).

Обобщение двухфазной модели (gleiss.py). Кольца от оси наружу:
    Ro (слив нефти) < r_i (граница нефть-вода) < Rtr (кек) < Rd (стенка)

Твёрдое подачи eps_s делится между несущими фазами долей frac_s_in_oil:
она задаётся отдельно, поэтому сумма eps_* остаётся равной 1 и
перебалансировать состав не нужно.

Дисперсных популяций три:
    * твёрдое в воде   (1-frac_s_in_oil) — кольцо r_i..Rtr, несущая вода,
      захват на кек;
    * твёрдое в нефти  (frac_s_in_oil)   — кольцо Ro..r_i, несущая нефть,
      захват на границу r_i; захваченное пересекает границу и достаётся
      водяному кольцу ТОЙ ЖЕ ячейки (обе жидкости текут соосно к сливам),
      дальше идёт по обычному водяному каскаду и попадает в кек;
    * капли воды в нефти (eps_wd)        — кольцо Ro..r_i, захват на r_i.

Осаждение твёрдого в нефти на ~28 раз медленнее, чем в воде
((rho_s-rho_o)/eta_o против (rho_s-rho_w)/eta_w), поэтому именно эта
популяция задаёт содержание мехпримесей в нефтяном продукте. Без неё
модель структурно выдаёт абсолютно чистую нефть при любом режиме.

Состояния — только масса кека, остальное алгебра.

Положение r_i задаётся балансом давлений двух сливов (U-трубка в поле
ω²r) и от кека не зависит. При eps_o = eps_wd = 0 получаем r_i = Rw,
frac_s_in_oil игнорируется, и модель в точности совпадает с двухфазной.
При frac_s_in_oil = 0 — совпадает с прежней трёхфазной версией.

Расходы нефти и воды считаются постоянными вдоль машины (дисперсии
разбавлены; переток твёрдого из нефти в воду расходы не меняет).
Слой на границе раздела (rag layer) не моделируется: перенесённое
твёрдое сразу становится частью водяного кольца.
"""
import numpy as np
from dataclasses import dataclass
import matplotlib
try:                       # на машине без tkinter остаётся бэкенд по умолчанию
    import tkinter as _tk  # noqa: F401
    matplotlib.use('TkAgg')
except Exception:
    pass
import matplotlib.pyplot as plt
G = 9.81

# скалярные наблюдаемые, которые пишутся в историю на каждом шаге
OUTS = ("U", "E_s", "E_so", "E_w", "E_s_tot",
        "phi_s", "phi_so", "phi_w", "bsw", "q_cake")


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
    J: int = 100

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
    eps_s: float = 0.02      # твёрдое, всего в подаче

    # какая доля твёрдого подачи пришла диспергированной в нефти.
    # Эмпирический вход (смачиваемость, история эмульсии), не выводится
    # из состава. 0 => прежняя версия модели.
    frac_s_in_oil: float = 0.25

    # куда попадает твёрдое, перешедшее из нефти в водяное кольцо:
    #   False — домешивается равномерно по сечению кольца; это то же
    #           допущение о перемешивании, на котором стоит вся ячеечная
    #           модель (основной вариант);
    #   True  — плуг: частицы стартуют с r_i и должны пройти всю толщину
    #           кольца, не перемешиваясь. Оценка снизу для захвата.
    inject_at_interface: bool = False

    mu_s: float = 0.25       # трение кека
    phi_sed: float = 0.45    # доля твёрдого в осадке
    phi_ref: float = 0.13    # опорная доля для твёрдого (в воде и в нефти)
    phi_max: float = 0.64    # плотная упаковка капель
    u_conv: float = 0.57     # эффективность конвейирования

    x50: float = 1.5e-6      # медиана частиц (твёрдое в воде), м
    bs: float = 3.0
    rrsb: bool = False       # форма распределения твёрдого
    x50_so: float = 0.0      # медиана твёрдого в нефти; 0 = как x50
    b_so: float = 0.0        # 0 = как bs
    rrsb_so: int = -1        # -1 = как rrsb
    d50: float = 3.0e-6      # медиана капель воды, м
    bw: float = 2.0

    # --- распределение твёрдого в нефти: по умолчанию наследует водяное ---
    @property
    def x50_o(self):
        return self.x50_so if self.x50_so > 0 else self.x50

    @property
    def b_o(self):
        return self.b_so if self.b_so > 0 else self.bs

    @property
    def rrsb_o(self):
        return bool(self.rrsb) if self.rrsb_so < 0 else bool(self.rrsb_so)

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

    # --- деление твёрдого между несущими фазами ---
    @property
    def f_so(self):
        """Доля твёрдого в нефти. Без нефтяной фазы диспергировать не в чем."""
        if not self.has_oil:
            return 0.0
        return float(np.clip(self.frac_s_in_oil, 0.0, 1.0))

    @property
    def eps_s_o(self):
        """Твёрдое, пришедшее внутри нефти (доля подачи)."""
        return self.eps_s * self.f_so

    @property
    def eps_s_w(self):
        """Твёрдое, пришедшее внутри свободной воды (доля подачи)."""
        return self.eps_s * (1.0 - self.f_so)

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
        return self.Q * (self.eps_o + self.eps_wd + self.eps_s_o)

    @property
    def Qw(self):
        return self.Q * (self.eps_w + self.eps_s_w)

    # --- входные доли дисперсных фаз в своих несущих ---
    @property
    def phi_s0(self):
        d = self.eps_w + self.eps_s_w
        return self.eps_s_w / d if d > 0 else 0.0

    @property
    def phi_so0(self):
        d = self.eps_o + self.eps_wd + self.eps_s_o
        return self.eps_s_o / d if d > 0 else 0.0

    @property
    def phi_w0(self):
        d = self.eps_o + self.eps_wd + self.eps_s_o
        return self.eps_wd / d if d > 0 else 0.0

    @property
    def eps_sum(self):
        return self.eps_o + self.eps_wd + self.eps_w + self.eps_s


# ---------------------------------------------------------------- PSD

def psd_edges(d_lo, d_hi, J):
    """Границы классов, логарифмически равномерно."""
    return np.logspace(np.log10(d_lo / 30), np.log10(d_hi * 30), J + 1)


def psd_weights(e, d50, b, rrsb=False):
    """Массовые доли классов на заданной сетке границ."""
    Q3 = 1 - np.exp(-(e / d50) ** b) if rrsb else 1 / (1 + (d50 / e) ** b)
    w = np.diff(Q3)
    return w / w.sum()


def psd(d50, b, J, rrsb=False):
    """Классы дисперсной фазы и их массовые доли."""
    e = psd_edges(d50, d50, J)
    return np.sqrt(e[:-1] * e[1:]), psd_weights(e, d50, b, rrsb)


def psd_solids(p):
    """Общая сетка классов для обеих твёрдых популяций.

    Общая — потому что переток из нефти в воду идёт поклассово (крупные
    уходят из нефти первыми), и индексы классов должны совпадать.
    При x50_so = x50 сетка совпадает с psd(x50, ...) до бита."""
    e = psd_edges(min(p.x50, p.x50_o), max(p.x50, p.x50_o), p.J)
    x = np.sqrt(e[:-1] * e[1:])
    return x, psd_weights(e, p.x50, p.bs, p.rrsb), \
        psd_weights(e, p.x50_o, p.b_o, p.rrsb_o)


# ---------------------------------------------------- осаждение, кинетика

def hindered(phi_sum, phi_ref):
    """Стеснение по Ричардсону-Заки."""
    return max(0.0, 1 - phi_sum / phi_ref) ** 4.65


def grade_cell(x, k0, R, tau, pref):
    """Эффективность захвата в одной ячейке; tau уже умножено на f_clar.
    Порядок умножений сохранён как в исходной settle() — тогда при
    eps_so = 0 результат совпадает со старой версией побитово."""
    return np.minimum(1.0, pref * (1 - np.exp(-2 * k0 * x ** 2 * R * tau)))


def grade(x, k0, tau, pref, phi_cells, phi_ref):
    """Матрица эффективности захвата T[ячейка, класс]."""
    R = np.maximum(0.0, 1 - phi_cells.sum(1) / phi_ref) ** 4.65
    arg = 2 * k0 * x[None, :] ** 2 * (R * tau)[:, None]
    return np.minimum(1.0, pref * (1 - np.exp(-arg)))


class Plug(object):
    """Твёрдое, перешедшее из нефти, при inject_at_interface = True.

    Частицы входят в водяное кольцо на r_i и не перемешиваются:
    r(t) = r0 exp(k t), захват когда радиус дошёл до поверхности кека.
    Внутри класса ведётся один представительный радиус — среднее по массе
    от ln r (траектории мультипликативны, поэтому усреднять надо именно
    логарифм). Это оценка снизу для захвата: в основном варианте те же
    частицы считаются размазанными по всему сечению кольца.

    Радиус ведётся как ln r: перемещение за ячейку тогда просто
    прибавляется, а усреднение по массе линейно и не переполняется."""

    __slots__ = ("on", "lr0", "q", "lr")

    def __init__(self, p, J):
        self.on = bool(p.inject_at_interface)
        self.lr0 = np.log(p.r_i)
        self.q = np.zeros(J)                # масса в долях водяного расхода
        self.lr = np.full(J, self.lr0)      # ln текущего радиуса класса

    def step(self, inj, travel, Rtr_i):
        """inj — приход в эту ячейку, travel = k0 x^2 R tau f_clar (= Δln r).
        Возвращает захваченное на кек в этой ячейке."""
        tot = self.q + inj
        safe = np.maximum(tot, 1e-300)
        self.lr = (self.q * self.lr + inj * self.lr0) / safe
        self.q = tot
        wall = np.log(Rtr_i)
        self.lr = np.minimum(self.lr + travel, wall)
        done = self.lr >= wall - 1e-15
        cap = float(np.where(done, self.q, 0.0).sum())
        self.q = np.where(done, 0.0, self.q)
        self.lr = np.where(done, self.lr0, self.lr)
        return cap

    @property
    def out(self):
        """Непойманное, уходит с водой."""
        return float(self.q.sum()) if self.on else 0.0


def geometry(p, m):
    """Кольца, времена пребывания, кинетика и префакторы по массе кека.

    Общее для квазистатики и динамики. rout в префакторах — фиксированная
    геометрическая граница (Rd для твёрдого в воде, r_i для всего, что
    оседает в нефти), не текущий радиус осадка."""
    Rtr = np.sqrt(np.maximum(
        p.Rd ** 2 - m / (p.rho_s * p.phi_sed * np.pi * p.Lax), p.r_i ** 2))
    dH = p.Rd - Rtr

    # твёрдое в воде: кольцо r_i .. Rtr, захват на кек
    tau_s = np.pi * (Rtr ** 2 - p.r_i ** 2) * p.Lax / p.Qw
    k_s = (p.rho_s - p.rho_w) * p.omega ** 2 / (18 * p.eta_w)
    pref_s = p.Rd ** 2 / max(p.Rd ** 2 - p.r_i ** 2, 1e-12)

    # нефтяное кольцо Ro .. r_i: твёрдое и капли воды, обе на границу r_i
    oil = p.has_oil and p.r_i - p.Ro > 1e-6
    if oil:
        tau_o = np.full(p.n, np.pi * (p.r_i ** 2 - p.Ro ** 2) * p.Lax / p.Qo)
        k_so = (p.rho_s - p.rho_o) * p.omega ** 2 / (18 * p.eta_o)
        k_w = (p.rho_w - p.rho_o) * p.omega ** 2 / (18 * p.eta_o)
        pref_o = p.r_i ** 2 / max(p.r_i ** 2 - p.Ro ** 2, 1e-12)
    else:
        tau_o = np.zeros(p.n)
        k_so = k_w = pref_o = 0.0
    return Rtr, dH, tau_s, tau_o, k_s, k_so, k_w, pref_s, pref_o, oil


def outputs(p, Rtr, dH, phi_s_out, phi_so_out, phi_w_out, q_cake):
    """Наблюдаемые по концентрациям на выходе и геометрии кека.

    Если у популяции нет входа (phi*0 = 0), её эффективность по соглашению
    выдаётся нулём, а не NaN: так же вёл себя E_w в двухфазном режиме."""
    feed_s = p.eps_s                                    # твёрдое подачи, всего
    esc = p.Qw * phi_s_out + p.Qo * phi_so_out          # унос твёрдого, м3/с
    return dict(
        q_cake=q_cake,        # твёрдое в выгрузке кека, м3/с
        U=np.mean(p.Rd ** 2 - Rtr ** 2) / max(p.Rd ** 2 - p.r_i ** 2, 1e-12),
        # E_s — осветление водяного потока относительно его же входа;
        # при заметном frac_s_in_oil может стать отрицательным: нефть
        # досыпает твёрдого в воду. Честная общая величина — E_s_tot.
        E_s=1 - phi_s_out / p.phi_s0 if p.phi_s0 else 0.0,
        E_so=1 - phi_so_out / p.phi_so0 if p.phi_so0 else 0.0,
        E_w=1 - phi_w_out / p.phi_w0 if p.phi_w0 else 0.0,
        E_s_tot=1 - esc / (p.Q * feed_s) if feed_s else 0.0,
        phi_s=phi_s_out,      # твёрдое в воде на выходе, об. доля
        phi_so=phi_so_out,    # твёрдое в нефти на выходе, об. доля
        phi_w=phi_w_out,      # вода в нефти на выходе, об. доля
        bsw=phi_w_out + phi_so_out,
        dH=dH)


# ------------------------------------------------------------- каскад

def cascade(p, m, xs, ws, wso, xw, ww):
    """По массе кека возвращает dm/dt и наблюдаемые величины.

    Один проход по ячейкам: в каждой сначала нефтяное кольцо (твёрдое и
    капли), потом водяное, куда сразу впадает перенесённое твёрдое."""
    Rtr, dH, tau_s, tau_o, k_s, k_so, k_w, pref_s, pref_o, oil = geometry(p, m)
    fc = p.f_clar
    q = p.Qo / p.Qw if oil else 0.0     # пересчёт доли: расход нефти -> воды

    up_s = p.phi_s0 * ws                # вход первой ячейки, водяное кольцо
    up_so = p.phi_so0 * wso             # то же, нефтяное кольцо
    up_w = p.phi_w0 * ww
    zero = np.zeros_like(xs)
    cap_s = np.zeros(p.n)
    plug = Plug(p, p.J)

    for i in range(p.n):
        if oil:
            # твёрдое в нефти -> граница r_i
            R = hindered(up_so.sum(), p.phi_ref)
            T = grade_cell(xs, k_so, R, tau_o[i] * fc, pref_o)
            moved = up_so * T
            up_so = up_so - moved
            # капли воды в нефти -> та же граница
            Rw_ = hindered(up_w.sum(), p.phi_max)
            Tw = grade_cell(xw, k_w, Rw_, tau_o[i] * fc, pref_o)
            up_w = up_w * (1 - Tw)
        else:
            moved = zero
        inj = moved * q                 # пересчёт доли: расход нефти -> воды
        if Rtr[i] - p.r_i < 1e-6:       # водяное кольцо схлопнулось:
            cap_s[i] = inj.sum()        # граница совпала с поверхностью кека
            inj = zero
        else:
            cap_s[i] = 0.0
        # водяное кольцо: свой вход + пришедшее из нефти (домешивается
        # равномерно по сечению кольца — то же допущение, на котором
        # стоит вся ячеечная модель; при inject_at_interface приход
        # обрабатывается отдельно, как непереме́шанный плуг)
        cin = up_s if plug.on else up_s + inj
        R = hindered(cin.sum() + plug.out, p.phi_ref)
        T = grade_cell(xs, k_s, R, tau_s[i] * fc, pref_s)
        cap_s[i] += (cin * T).sum()
        up_s = cin * (1 - T)
        if plug.on:
            cap_s[i] += plug.step(inj, k_s * xs ** 2 * R * tau_s[i] * fc,
                                  Rtr[i])

    tr = p.phi_sed * p.rho_s * np.pi * (p.Rd ** 2 - Rtr ** 2) * p.u_ax
    dm = p.rho_s * p.Qw * cap_s + np.append(tr[1:], 0.0) - tr
    return dm, outputs(p, Rtr, dH, up_s.sum() + plug.out, up_so.sum(),
                       up_w.sum(), tr[0] / p.rho_s)


def simulate_dyn(p, t_end, dt, inputs=None):
    """Концентрации в состояниях: даёт запаздывание, потерянное в квазистатике."""
    t = np.arange(0, t_end + dt, dt)
    live_psd = bool(inputs) and bool(
        {"x50", "bs", "rrsb", "x50_so", "b_so", "rrsb_so",
         "d50", "bw", "J"} & inputs.keys())
    xs, ws, wso = psd_solids(p)
    xw, ww = psd(p.d50, p.bw, p.J)

    m = np.zeros(p.n)
    phi_s = np.zeros((p.n, p.J))        # твёрдое в воде
    phi_so = np.zeros((p.n, p.J))       # твёрдое в нефти
    phi_w = np.zeros((p.n, p.J))        # капли воды в нефти
    hist = {k: np.zeros(t.size) for k in OUTS}
    dH_hist = np.zeros((t.size, p.n))

    for i in range(t.size):
        if inputs:
            for name, fn in inputs.items():
                setattr(p, name, fn(t[i]))
            if live_psd:
                xs, ws, wso = psd_solids(p)
                xw, ww = psd(p.d50, p.bw, p.J)

        Rtr, dH, tau_s, tau_o, k_s, k_so, k_w, pref_s, pref_o, oil = \
            geometry(p, m)
        fc = p.f_clar
        q = p.Qo / p.Qw if oil else 0.0

        T_s = grade(xs, k_s, tau_s * fc, pref_s, phi_s, p.phi_ref)
        decay_s = np.exp(-dt / np.maximum(tau_s, 1e-9))
        if oil:
            T_so = grade(xs, k_so, tau_o * fc, pref_o, phi_so, p.phi_ref)
            T_w = grade(xw, k_w, tau_o * fc, pref_o, phi_w, p.phi_max)
            decay_o = np.exp(-dt / np.maximum(tau_o, 1e-9))


        #plt.plot(phi_so)
        #plt.show()

        cap_s = np.zeros(p.n)
        s_next = np.empty_like(phi_s)
        so_next = np.empty_like(phi_so)
        w_next = np.empty_like(phi_w)
        zero = np.zeros(p.J)
        up_s = p.phi_s0 * ws
        up_so = p.phi_so0 * wso
        up_w = p.phi_w0 * ww
        plug = Plug(p, p.J)             # плуг-фракция считается квазистатически
        for c in range(p.n):           # c, не i: i занят внешним циклом по времени
            if oil:
                moved = up_so * T_so[c]
                star = up_so - moved
                so_next[c] = star + (phi_so[c] - star) * decay_o[c]
                star_w = up_w * (1 - T_w[c])
                w_next[c] = star_w + (phi_w[c] - star_w) * decay_o[c]
                up_so, up_w = so_next[c], w_next[c]
            else:
                moved = zero
            inj = moved * q
            if Rtr[c] - p.r_i < 1e-6:
                cap_s[c] = inj.sum()
                inj = zero
            else:
                cap_s[c] = 0.0
            cin = up_s if plug.on else up_s + inj
            cap_s[c] += (cin * T_s[c]).sum()
            star_c = cin * (1 - T_s[c])
            s_next[c] = star_c + (phi_s[c] - star_c) * decay_s[c]
            up_s = s_next[c]            # уже на текущем шаге, не на прошлом
            if plug.on:
                R = hindered(cin.sum() + plug.out, p.phi_ref)
                cap_s[c] += plug.step(inj, k_s * xs ** 2 * R * tau_s[c] * fc,
                                      Rtr[c])
        phi_s = s_next
        if oil:
            phi_so, phi_w = so_next, w_next
            out_so, out_w = phi_so[-1].sum(), phi_w[-1].sum()
        else:                            # нефтяного кольца нет — проскок
            out_so, out_w = p.phi_so0, p.phi_w0

        tr = p.phi_sed * p.rho_s * np.pi * (p.Rd ** 2 - Rtr ** 2) * p.u_ax
        m = np.maximum(m + dt * (p.rho_s * p.Qw * cap_s
                                 + np.append(tr[1:], 0.0) - tr), 0.0)

        out = outputs(p, Rtr, dH, phi_s[-1].sum() + plug.out, out_so, out_w,
                      tr[0] / p.rho_s)
        for k in OUTS:
            hist[k][i] = out[k]
        dH_hist[i] = dH
    plt.show()
    return dict(t=t, dH=dH_hist, **hist)



def step(t0, before, after):
    """Ступенька: значение before до t0, after после."""
    return lambda t: before if t < t0 else after


def simulate(p, t_end=1500.0, dt=1.0, inputs=None, dynamic=True):
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
        {"x50", "bs", "rrsb", "x50_so", "b_so", "rrsb_so",
         "d50", "bw", "J"} & inputs.keys())
    xs, ws, wso = psd_solids(p)
    xw, ww = psd(p.d50, p.bw, p.J)
    m = np.zeros(p.n)
    hist = {k: np.zeros(t.size) for k in OUTS}
    dH = np.zeros((t.size, p.n))
    for i in range(t.size):
        if inputs:
            for name, fn in inputs.items():
                setattr(p, name, fn(t[i]))
            if live_psd:
                xs, ws, wso = psd_solids(p)
                xw, ww = psd(p.d50, p.bw, p.J)
        dm, out = cascade(p, m, xs, ws, wso, xw, ww)
        for k in OUTS:
            hist[k][i] = out[k]
        dH[i] = out['dH']
        m = np.maximum(m + dt * dm, 0.0)
    return dict(t=t, dH=dH, **hist)


def weir_for_interface(p, r_i):
    """Радиус водяного слива, дающий заданную границу нефть-вода."""
    return np.sqrt((r_i ** 2 * (p.rho_w - p.rho_o) + p.rho_o * p.Ro ** 2) / p.rho_w)


def balance(p, res, i=-1):
    """Невязка баланса твёрдого на срезе i (доля от подачи).

    Подача = унос с водой + унос с нефтью + выгрузка кека. В установившемся
    режиме должна быть ~0; на переходе положительна (машина накапливает)."""
    feed = p.Q * p.eps_s
    esc = p.Qw * res['phi_s'][i] + p.Qo * res['phi_so'][i]
    return (feed - esc - res['q_cake'][i]) / feed if feed else 0.0


def check_limits(p):
    """Проверка, что режим физически существует."""
    msg = []
    if abs(p.eps_sum - 1.0) > 1e-9:
        msg.append(f"доли подачи не дают 1: sum = {p.eps_sum:.4f}")
    if not 0.0 <= p.frac_s_in_oil <= 1.0:
        msg.append("frac_s_in_oil вне [0, 1], зажимается")
    if not p.has_oil:
        if p.frac_s_in_oil > 0:
            msg.append("нефти нет: frac_s_in_oil игнорируется, "
                       "всё твёрдое считается пришедшим с водой")
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
               frac_s_in_oil=0.0,
               rho_s=1410, C=250, dn=5, mu_s=0.3, phi_sed=0.5,
               phi_ref=0.5, x50=2.5e-6, bs=2.78, rrsb=True, Rw=0.034)
    for k, v in kw.items():
        setattr(p, k, v)
    return p


def plot(p=None, t_end=1500.0):
    p = p or Params()
    r = simulate(p, t_end)
    s = np.arange(p.n) * p.Lax

    fig, ax = plt.subplots(1, 5, figsize=(20, 3.6))
    ax[0].plot(r['t'], r['U'])
    ax[0].set(xlabel="t, с", ylabel="U")
    ax[1].plot(r['t'], r['E_s'], color="g", label="в воде")
    ax[1].plot(r['t'], r['E_s_tot'], color="k", ls="--", label="всего")
    ax[1].set(xlabel="t, с", ylabel="E твёрдого")
    ax[1].legend(fontsize=8)
    ax[2].plot(r['t'], r['phi_w'] * 100, color="b")
    ax[2].set(xlabel="t, с", ylabel="вода в нефти, %об")
    ax[3].plot(r['t'], r['phi_so'] * 100, color="r")
    ax[3].set(xlabel="t, с", ylabel="твёрдое в нефти, %об")
    for i in range(0, r['t'].size, max(1, r['t'].size // 12)):
        ax[4].plot(s, r['dH'][i] * 1000, color=plt.cm.viridis(i / r['t'].size), lw=1)
    ax[4].set(xlabel="осевая координата, м", ylabel="высота осадка, мм")
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

    r = simulate(p, 1500)
    print(f"\nтрёхфаза, frac_s_in_oil = {p.frac_s_in_oil}:")
    print(f"  U = {r['U'][-1]:.3f}   E твёрдого всего = {r['E_s_tot'][-1]:.3f}"
          f"   (вода {r['E_s'][-1]:.3f}, нефть {r['E_so'][-1]:.3f})")
    print(f"  вода в нефти {r['phi_w'][-1]*100:.2f}% (было {p.phi_w0*100:.2f}%),"
          f"  твёрдое в нефти {r['phi_so'][-1]*100:.3f}%"
          f" (было {p.phi_so0*100:.3f}%)")
    print(f"  BS&W нефти {r['bsw'][-1]*100:.2f}%,"
          f"  невязка баланса твёрдого {balance(p, r):.1e}")

    # доля твёрдого в нефти задаётся отдельно, сумма eps_* не трогается
    print(f"\nподъём frac_s_in_oil (C = {Params().C:.0f} g):")
    print(" frac     E_so   тв.нефть,%  E_s(вода)  E_s общая   BS&W,%     U"
          "     dH1,мкм   невязка")
    for f in (0.0, 0.1, 0.25, 0.5, 1.0):
        q = Params(frac_s_in_oil=f)
        rr = simulate(q, 1500)
        print(f" {f:5.2f} {rr['E_so'][-1]:8.4f} {rr['phi_so'][-1]*100:10.4f}"
              f" {rr['E_s'][-1]:10.4f} {rr['E_s_tot'][-1]:10.4f}"
              f" {rr['bsw'][-1]*100:9.3f} {rr['U'][-1]:7.4f}"
              f" {rr['dH'][-1,0]*1e6:9.2f} {balance(q, rr):9.1e}")

    # главный эффект: нефть чистится на порядок хуже воды при любом C
    print("\nта же frac_s_in_oil = 0.25 по фактору разделения:")
    print("   C, g     E_so   E_s(вода)  E_s общая  тв.нефть,%  вода в нефти,%")
    for C in (25., 100., 500., 1000., 2000., 3000.):
        q = Params(C=C)
        rr = simulate(q, 1500)
        print(f" {C:6.0f} {rr['E_so'][-1]:8.4f} {rr['E_s'][-1]:10.4f}"
              f" {rr['E_s_tot'][-1]:10.4f} {rr['phi_so'][-1]*100:11.4f}"
              f" {rr['phi_w'][-1]*100:14.3f}")

    # оценка снизу: перенесённое твёрдое не перемешивается по кольцу
    rp = simulate(Params(inject_at_interface=True), 1500)
    print(f"\nплуг-вариант (inject_at_interface): E общая {rp['E_s_tot'][-1]:.3f}"
          f" против {r['E_s_tot'][-1]:.3f} при перемешивании")

    r2 = simulate(two_phase(), 1500)
    print(f"\nдвухфаза:  U={r2['U'][-1]:.3f}  E={r2['E_s'][-1]:.3f}  "
          f"осадок у входа {r2['dH'][-1, 0]*1000:.1f} мм")

    plot(p, 1500)
    plt.show()
    # ------------------------------------------------------------
    # Масса осадка (кека) в каждой ячейке в установившемся режиме
    # ------------------------------------------------------------

    dH_final = r["dH"][-1]          # толщина кека в последний момент времени, м
    Rtr_final = p.Rd - dH_final     # радиус поверхности кека, м

    m_cake = (
        p.rho_s
        * p.phi_sed
        * np.pi
        * (p.Rd**2 - Rtr_final**2)
        * p.Lax
    )                               # кг в каждой ячейке

    cells = np.arange(1, p.n + 1)

    plt.figure(figsize=(10, 5))
    plt.bar(cells, m_cake * 1000)   # перевод кг -> г

    plt.xlabel("Номер ячейки")
    plt.ylabel("Масса кека, г")
    plt.title("Распределение массы осадка по ячейкам")
    plt.xticks(cells)
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.show()

    print("\nМасса кека по ячейкам:")
    for i, m in enumerate(m_cake, start=1):
        print(f"ячейка {i:2d}: {m*1000:.4f} г")

    print(f"\nОбщая масса кека в машине: {m_cake.sum()*1000:.3f} г")


    

    


