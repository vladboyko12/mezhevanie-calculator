"""
Перевод контуров-ломаных в кривые Безье.

Контур делится на участки по острым углам, каждый участок аппроксимируется
кубическими кривыми методом Шнайдера («An Algorithm for Automatically Fitting
Digitized Curves», Graphics Gems, 1990). Почти прямые кривые сохраняются
как отрезки. В TTF кубики переводятся в квадратичные через Cu2QuPen.
"""
import math

import numpy as np

CORNER_DEG = 40      # излом больше этого угла — точка угла
FIT_ERROR = 0.8      # допустимое отклонение кривой от исходного контура, ед.
LINE_TOL = 0.35      # кривую, почти совпадающую с хордой, пишем отрезком
LONG_EDGE = 25       # ребро длиннее — точно прямая, не аппроксимируем
DENSE = 5            # шаг уплотнения точек для проверки ошибки


def _densify(seg):
    out = [seg[0]]
    for a, b in zip(seg, seg[1:]):
        k = max(1, int(np.ceil(np.linalg.norm(b - a) / DENSE)))
        for j in range(1, k + 1):
            out.append(a + (b - a) * j / k)
    return np.array(out)


def _norm(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else v


def _bez(ctrl, t):
    t = np.asarray(t)[:, None]
    p0, p1, p2, p3 = ctrl
    return ((1 - t) ** 3) * p0 + 3 * ((1 - t) ** 2) * t * p1 + 3 * (1 - t) * t ** 2 * p2 + t ** 3 * p3


def _chord_params(pts):
    d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))]
    return d / d[-1] if d[-1] > 0 else d


def _generate(pts, u, t1, t2):
    p0, p3 = pts[0], pts[-1]
    a1 = t1[None, :] * (3 * (1 - u) ** 2 * u)[:, None]
    a2 = t2[None, :] * (3 * (1 - u) * u ** 2)[:, None]
    c00, c01, c11 = (a1 * a1).sum(), (a1 * a2).sum(), (a2 * a2).sum()
    base = _bez([p0, p0, p3, p3], u)
    tmp = pts - base
    x0, x1 = (a1 * tmp).sum(), (a2 * tmp).sum()
    det = c00 * c11 - c01 * c01
    seg = np.linalg.norm(p3 - p0)
    if abs(det) > 1e-12:
        al = (x0 * c11 - x1 * c01) / det
        ar = (c00 * x1 - c01 * x0) / det
    else:
        al = ar = 0
    eps = 1e-6 * seg
    if al < eps or ar < eps:
        al = ar = seg / 3
    return [p0, p0 + t1 * al, p3 + t2 * ar, p3]


def _newton(ctrl, pts, u):
    p0, p1, p2, p3 = ctrl
    d1 = [3 * (p1 - p0), 3 * (p2 - p1), 3 * (p3 - p2)]
    d2 = [2 * (d1[1] - d1[0]), 2 * (d1[2] - d1[1])]
    out = []
    for p, t in zip(pts, u):
        q = _bez(ctrl, [t])[0]
        q1 = (1 - t) ** 2 * d1[0] + 2 * (1 - t) * t * d1[1] + t ** 2 * d1[2]
        q2 = (1 - t) * d2[0] + t * d2[1]
        num = ((q - p) * q1).sum()
        den = (q1 * q1).sum() + ((q - p) * q2).sum()
        out.append(t - num / den if abs(den) > 1e-12 else t)
    return np.clip(np.array(out), 0, 1)


def _max_err(ctrl, pts, u):
    d = np.linalg.norm(_bez(ctrl, u) - pts, axis=1)
    i = int(np.argmax(d))
    return d[i], i


def _fit(pts, t1, t2, err, depth=0):
    if len(pts) == 2:
        d = np.linalg.norm(pts[1] - pts[0]) / 3
        return [[pts[0], pts[0] + t1 * d, pts[1] + t2 * d, pts[1]]]
    u = _chord_params(pts)
    ctrl = _generate(pts, u, t1, t2)
    e, split = _max_err(ctrl, pts, u)
    if e < err:
        return [ctrl]
    if e < err * 4:
        for _ in range(20):
            u = _newton(ctrl, pts, u)
            ctrl = _generate(pts, u, t1, t2)
            e, split = _max_err(ctrl, pts, u)
            if e < err:
                return [ctrl]
    if depth > 30:
        return [ctrl]
    split = min(max(split, 1), len(pts) - 2)
    tc = _norm(pts[split - 1] - pts[split + 1])
    return (_fit(pts[:split + 1], t1, tc, err, depth + 1) +
            _fit(pts[split:], -tc, t2, err, depth + 1))


def _is_straight(ctrl):
    p0, p1, p2, p3 = ctrl
    ch = p3 - p0
    L = np.linalg.norm(ch)
    if L < 1e-9:
        return True
    n = np.array([-ch[1], ch[0]]) / L

    def off(p):
        return abs(((p - p0) * n).sum())

    def inside(p):
        t = ((p - p0) * ch).sum() / (L * L)
        return -0.02 <= t <= 1.02
    return off(p1) < LINE_TOL and off(p2) < LINE_TOL and inside(p1) and inside(p2)


def fit_ring(coords):
    """coords — замкнутый список точек (без повтора первой). Возвращает
    список команд: ('line', p) / ('curve', c1, c2, p), начиная с точки start."""
    pts = np.array(coords, dtype=float)
    # убрать дубликаты
    keep = [0]
    for i in range(1, len(pts)):
        if np.linalg.norm(pts[i] - pts[keep[-1]]) > 1e-6:
            keep.append(i)
    pts = pts[keep]
    if len(pts) > 1 and np.linalg.norm(pts[0] - pts[-1]) < 1e-6:
        pts = pts[:-1]
    n = len(pts)
    if n < 3:
        return pts[0], []

    def turn(i):
        a = _norm(pts[i] - pts[i - 1])
        b = _norm(pts[(i + 1) % n] - pts[i])
        return math.degrees(math.acos(max(-1, min(1, (a * b).sum()))))

    angles = [turn(i) for i in range(n)]
    corners = {i for i in range(n) if angles[i] > CORNER_DEG}
    # длинные рёбра — это настоящие прямые: их концы тоже точки разбиения
    for i in range(n):
        if np.linalg.norm(pts[(i + 1) % n] - pts[i]) > LONG_EDGE:
            corners.add(i)
            corners.add((i + 1) % n)
    corners = sorted(corners)
    if not corners:
        corners = [int(np.argmax(angles))]
        smooth_start = True
    else:
        smooth_start = False
    start = corners[0]
    pts = np.roll(pts, -start, axis=0)
    corners = sorted((c - start) % n for c in corners)
    cmds = []
    bounds = corners + [n]
    for a, b in zip(bounds, bounds[1:]):
        seg = np.array([pts[i % n] for i in range(a, b + 1)])
        if len(seg) < 2:
            continue
        if len(seg) == 2:
            cmds.append(("line", seg[1]))
            continue
        if smooth_start:
            t1 = _norm(seg[1] - seg[-2])
            t2 = -t1
        else:
            t1 = _norm(seg[1] - seg[0])
            t2 = _norm(seg[-2] - seg[-1])
        for c in _fit(_densify(seg), t1, t2, FIT_ERROR):
            if _is_straight(c):
                cmds.append(("line", c[3]))
            else:
                cmds.append(("curve", c[1], c[2], c[3]))
    return pts[0], cmds


def _split(ctrl, t):
    p0, p1, p2, p3 = ctrl
    a, b, c = p0 + (p1 - p0) * t, p1 + (p2 - p1) * t, p2 + (p3 - p2) * t
    d, e = a + (b - a) * t, b + (c - b) * t
    f = d + (e - d) * t
    return [p0, a, d, f], [f, e, c, p3]


def _extrema_t(ctrl):
    """Параметры t, где кривая достигает экстремума по x или y."""
    p0, p1, p2, p3 = ctrl
    ts = []
    for k in (0, 1):
        a = -p0[k] + 3 * p1[k] - 3 * p2[k] + p3[k]
        b = 2 * (p0[k] - 2 * p1[k] + p2[k])
        c = p1[k] - p0[k]
        if abs(a) < 1e-9:
            if abs(b) > 1e-9:
                ts.append(-c / b)
        else:
            disc = b * b - 4 * a * c
            if disc >= 0:
                s = math.sqrt(disc)
                ts += [(-b + s) / (2 * a), (-b - s) / (2 * a)]
    return sorted(t for t in ts if 0.02 < t < 0.98)


def _with_extrema(cmds, start):
    """Опорные точки на экстремумах: габариты контура = габариты кривой."""
    out, cur = [], np.array(start)
    for c in cmds:
        if c[0] == "line":
            out.append(c)
            cur = c[1]
            continue
        ctrl = [cur, c[1], c[2], c[3]]
        prev = 0.0
        for t in _extrema_t(ctrl):
            left, ctrl = _split(ctrl, (t - prev) / (1 - prev))
            prev = t
            out.append(("curve", left[1], left[2], left[3]))
        out.append(("curve", ctrl[1], ctrl[2], ctrl[3]))
        cur = c[3]
    return out


def draw_ring(coords, pen, shift=0.0):
    start, cmds = fit_ring(coords)
    if not cmds:
        return
    cmds = _with_extrema(cmds, start)
    sh = np.array([shift, 0.0])
    q = lambda p: (int(round(p[0] + sh[0])), int(round(p[1] + sh[1])))   # сразу в целые координаты шрифта
    first = cur = q(start)
    pen.moveTo(first)
    for i, c in enumerate(cmds):
        last = i == len(cmds) - 1
        if c[0] == "line":
            pt = q(c[1])
            # последний отрезок замыкает closePath; отрезки нулевой длины не пишем
            if not last and pt != cur:
                pen.lineTo(pt)
                cur = pt
        else:
            p1, p2, p3 = q(c[1]), q(c[2]), q(c[3])
            if p1 == p2 == p3 == cur:
                continue
            pen.curveTo(p1, p2, p3)
            cur = p3
    pen.closePath()
