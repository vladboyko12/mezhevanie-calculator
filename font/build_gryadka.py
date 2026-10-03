#!/usr/bin/env python3
"""
Gryadka (Грядка) — жирный геометрический плакатный гротеск, близкий по
устройству к Villula: на месте прописных — очень широкие буквы, на месте
строчных — узкие капительные, всё одной высоты. Контуры нарисованы с нуля.

Свои изюминки:
  1. Мягкие углы — все углы и срезы чуть скруглены (у Villula они острые).
  2. Листики вместо точек над Ё, Й, ё, й, i, j (фирменная деталь семьи Kapelka).
  3. ss01 — строчные становятся широкими, ss02 — прописные узкими.
  4. Узкие «а» и «е» нарисованы как строчные, остальные — капитель.

Сборка:  pip install fonttools shapely pillow && python3 build_gryadka.py
"""
import math
import os

from shapely.geometry import Point, LineString, Polygon, box
from shapely.ops import unary_union
from shapely import affinity

from fontTools.fontBuilder import FontBuilder
from fontTools.pens.ttGlyphPen import TTGlyphPen
from fontTools.pens.t2CharStringPen import T2CharStringPen
from fontTools.feaLib.builder import addOpenTypeFeaturesFromString

from build_kapelka import ln, arc, P, leaf_poly, draw, polys_of

HERE = os.path.dirname(os.path.abspath(__file__))
FAMILY = "Gryadka"
UPM = 1000
CH = 700
R = 60            # половина толщины штриха (штрих 120)
SB = 48           # боковые отступы
SOFT = 16         # радиус «мягких углов»
T, B, M = CH - R, R, 350   # осевые линии верхней, нижней перекладин и середины


# ---------------------------------------------------------------- примитивы

def rr_ring(x0, y0, x1, y1):
    """Осевая линия «стадиона»/круга, вписанного в прямоугольник."""
    w, h = x1 - x0, y1 - y0
    r = min(w, h) / 2
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    if abs(w - h) < 1:
        core = Point(cx, cy)
    elif w > h:
        core = LineString([(x0 + r, cy), (x1 - r, cy)])
    else:
        core = LineString([(cx, y0 + r), (cx, y1 - r)])
    return core.buffer(r, quad_segs=48).exterior


def oval(x0, y0, x1, y1):
    return list(rr_ring(x0, y0, x1, y1).coords)


def oval_arc(x0, y0, x1, y1, a0, a1):
    ring = rr_ring(x0, y0, x1, y1)
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    n = max(8, int(abs(a1 - a0) / 1.5))
    out = []
    for i in range(n + 1):
        a = math.radians(a0 + (a1 - a0) * i / n)
        ray = LineString([(cx, cy), (cx + 3000 * math.cos(a), cy + 3000 * math.sin(a))])
        hit = ray.intersection(ring)
        if hit.is_empty:
            continue
        pts = [hit] if hit.geom_type == "Point" else list(hit.geoms)
        p = min(pts, key=lambda q: q.distance(Point(cx, cy)))
        out.append((p.x, p.y))
    return out


def rbowl(xl, ytop, ybot, xr):
    """Правая чаша: от (xl,ytop) вправо, скруглённо вниз, обратно к (xl,ybot)."""
    r = min((ytop - ybot) / 2, xr - xl)
    return P(ln((xl, ytop), (xr - r, ytop)),
             arc(xr - r, ytop - r, r, r, 90, 0),
             ln((xr, ytop - r), (xr, ybot + r)),
             arc(xr - r, ybot + r, r, r, 0, -90),
             ln((xr - r, ybot), (xl, ybot)))


def lbowl(xr, ytop, ybot, xl):
    """Левая чаша (зеркало rbowl)."""
    r = min((ytop - ybot) / 2, xr - xl)
    return P(ln((xr, ytop), (xl + r, ytop)),
             arc(xl + r, ytop - r, r, r, 90, 180),
             ln((xl, ytop - r), (xl, ybot + r)),
             arc(xl + r, ybot + r, r, r, 180, 270),
             ln((xl + r, ybot), (xr, ybot)))


def G(strokes, free=(), leaves=(), band=(0, CH)):
    return {"strokes": [list(s) for s in strokes], "free": [list(s) for s in free],
            "leaves": list(leaves), "band": band}


# ---------------------------------------------------------------- буквы
# Каждая буква — функция от внешней ширины w (без боковых отступов).

L = {}


def A_(w):
    xc = w / 2
    return G([ln((0, 0), (xc, CH), (w, 0)), ln((w * 0.2, 230), (w * 0.8, 230))])


def B_(w):
    x1 = w - R
    return G([ln((R, 0), (R, CH)),
              rbowl(R, T, M, x1 - 40), rbowl(R, M, B, x1)])


def Be_(w):
    x1 = w - R
    return G([P(ln((w - 20, T), (R, T)), ln((R, T), (R, 0))), rbowl(R, 400, B, x1)])


def Ge_(w):
    return G([P(ln((R, 0), (R, T)), ln((R, T), (w, T)))])


def De_(w):
    x1 = w - R
    return G([ln((R, -140), (R, B), (x1, B), (x1, -140)),
              ln((R + 50, B), (w * 0.33, T), (x1, T), (x1, B))], band=(-140, CH))


def E_(w):
    return G([P(ln((w, T), (R, T)), ln((R, T), (R, B)), ln((R, B), (w, B))), ln((R, M), (w - 50, M))])


def F_(w):
    return G([P(ln((w, T), (R, T)), ln((R, T), (R, 0))), ln((R, M), (w - 60, M))])


def Zhe_(w):
    xc = w / 2
    return G([ln((xc, 0), (xc, CH)),
              P(ln((0, CH), (xc, M + 10), (0, 0))), P(ln((w, CH), (xc, M + 10), (w, 0)))])


def Ze_(w):
    x1 = w - R
    return G([rbowl(30, T, M, x1 - 30), rbowl(w * 0.3, M, B, x1)])


def I_(w):
    x1 = w - R
    return G([ln((R, 0), (R, CH)), ln((x1, 0), (x1, CH)), ln((R + 30, 0), (x1 - 30, CH))])


def K_(w):
    return G([ln((R, 0), (R, CH)), ln((w - 30, CH), (R + 80, M), (w - 10, 0))])


def El_(w):
    x1 = w - R
    return G([P(ln((10, 0), (w * 0.3, T)), ln((w * 0.3, T), (x1, T)), ln((x1, T), (x1, 0)))])


def M_(w):
    x1 = w - R
    return G([P(ln((R, 0), (R, T)), ln((R, T), (w / 2, 130)), ln((w / 2, 130), (x1, T)), ln((x1, T), (x1, 0)))])


def N_cyr(w):
    x1 = w - R
    return G([ln((R, 0), (R, CH)), ln((x1, 0), (x1, CH)), ln((R, M), (x1, M))])


def O_(w):
    return G([oval(R, B, w - R, T)])


def Pe_(w):
    x1 = w - R
    return G([P(ln((R, 0), (R, T)), ln((R, T), (x1, T)), ln((x1, T), (x1, 0)))])


def Er_(w):
    return G([ln((R, 0), (R, T)), rbowl(R, T, 290, w - R)])


def Es_(w):
    return G([oval_arc(R, B, w - R, T, 38, 322)])


def Te_(w):
    return G([ln((0, T), (w, T)), ln((w / 2, T), (w / 2, 0))])


def U_cyr(w):
    return G([ln((0, CH), (w * 0.47, 290)), ln((w, CH), (w * 0.28, 0))])


def Ef_(w):
    xc = w / 2
    return G([ln((xc, 0), (xc, CH)), oval(R, 150, w - R, 560)])


def Kha_(w):
    return G([ln((0, CH), (w, 0)), ln((w, CH), (0, 0))])


def Tse_(w):
    x1 = w - R - 70
    return G([P(ln((R, CH), (R, B)), ln((R, B), (w - R, B)), ln((w - R, B), (w - R, -140))),
              ln((x1, CH), (x1, B))], band=(-140, CH))


def Che_(w):
    x1 = w - R
    r = 150
    return G([P(ln((R, CH), (R, 300 + r)), arc(R + r, 300 + r, r, r, 180, 270), ln((R + r, 300), (x1, 300))),
              ln((x1, 0), (x1, CH))])


def Sha_(w):
    x1 = w - R
    return G([P(ln((R, CH), (R, B)), ln((R, B), (x1, B)), ln((x1, B), (x1, CH))), ln((w / 2, CH), (w / 2, B))])


def Shcha_(w):
    x1 = w - R - 70
    return G([P(ln((R, CH), (R, B)), ln((R, B), (w - R, B)), ln((w - R, B), (w - R, -140))),
              ln((x1, CH), (x1, B)), ln(((R + x1) / 2, CH), ((R + x1) / 2, B))], band=(-140, CH))


def soft_(w, x0=R):
    return [ln((x0, 0), (x0, CH)), rbowl(x0, 410, B, w - R)]


def Soft_(w):
    return G(soft_(w))


def Hard_(w):
    x0 = w * 0.3
    return G([ln((0, T), (x0, T))] + soft_(w, x0))


def Yeru_(w):
    sw = w * 0.68
    return G(soft_(sw) + [ln((w - R, 0), (w - R, CH))])


def Ee_(w):
    return G([oval_arc(R, B, w - R, T, 142, -142), ln((w * 0.35, M), (w - R, M))])


def Yu_(w):
    return G([ln((R, 0), (R, CH)), ln((R, M), (w * 0.36, M)), oval(w * 0.36, B, w - R, T)])


def Ya_(w):
    x1 = w - R
    return G([ln((x1, 0), (x1, T)), lbowl(x1, T, 290, R), ln((w * 0.48, 290), (0, 0))])


def D_lat(w):
    return G([ln((R, 0), (R, CH)), rbowl(R, T, B, w - R)])


def G_lat(w):
    x1 = w - R
    return G([P(oval_arc(R, B, x1, T, 42, 360), ln((x1, M), (w * 0.55, M)))])


def I_lat(w):
    return G([ln((R, 0), (R, CH))])


def J_lat(w):
    x1 = w - R
    d = x1 - R
    return G([P(ln((x1, CH), (x1, B + d / 2)), oval_arc(R, B, x1, B + d, 0, -180))])


def L_lat(w):
    return G([P(ln((R, CH), (R, B)), ln((R, B), (w, B)))])


def N_lat(w):
    x1 = w - R
    return G([ln((R, 0), (R, CH)), ln((x1, 0), (x1, CH)), ln((R + 30, CH), (x1 - 30, 0))])


def Q_(w):
    return G([oval(R, B, w - R, T), ln((w * 0.6, 200), (w, 0))])


def R_lat(w):
    return G([ln((R, 0), (R, T)), rbowl(R, T, 290, w - R), ln((w * 0.45, 290), (w, 0))])


def S_(w):
    x1 = w - R
    return G([P(oval_arc(R, M, x1, T, 15, 270), oval_arc(R, B, x1, M, 90, -165))])


def U_lat(w):
    x1 = w - R
    return G([P(ln((R, CH), (R, M)), oval_arc(R, B, x1, T, 180, 360), ln((x1, M), (x1, CH)))])


def V_(w):
    return G([ln((0, CH), (w / 2, 0), (w, CH))])


def W_(w):
    return G([ln((0, CH), (w * 0.25, 0), (w / 2, CH), (w * 0.75, 0), (w, CH))])


def Y_lat(w):
    xc = w / 2
    return G([ln((0, CH), (xc, M), (w, CH)), ln((xc, M), (xc, 0))])


def Z_(w):
    return G([ln((30, T), (w - 30, T), (40, B), (w, B))])


# узкие «строчные» формы, которые не совпадают с капителью
def a_low(w):
    x1 = w - R
    return G([oval(R, B, x1, T), ln((x1, 0), (x1, CH))])


def e_low(w):
    x1 = w - R
    return G([P(ln((R, M), (x1, M)), oval_arc(R, B, x1, T, 0, 322))])


# ---------------------------------------------------------------- таблица: символ → (функция, широкая, узкая)

CYR = [
    ("А", "а", A_, 900, 440), ("Б", "б", Be_, 760, 420), ("В", "в", B_, 760, 410),
    ("Г", "г", Ge_, 640, 360), ("Д", "д", De_, 840, 480), ("Е", "е", E_, 680, 380),
    ("Ж", "ж", Zhe_, 1080, 660), ("З", "з", Ze_, 700, 400), ("И", "и", I_, 800, 440),
    ("К", "к", K_, 780, 430), ("Л", "л", El_, 820, 460), ("М", "м", M_, 980, 580),
    ("Н", "н", N_cyr, 780, 440), ("О", "о", O_, 760, 440), ("П", "п", Pe_, 780, 440),
    ("Р", "р", Er_, 740, 410), ("С", "с", Es_, 740, 430), ("Т", "т", Te_, 780, 420),
    ("У", "у", U_cyr, 820, 460), ("Ф", "ф", Ef_, 960, 620), ("Х", "х", Kha_, 860, 460),
    ("Ц", "ц", Tse_, 880, 500), ("Ч", "ч", Che_, 720, 420), ("Ш", "ш", Sha_, 1000, 640),
    ("Щ", "щ", Shcha_, 1080, 720), ("Ъ", "ъ", Hard_, 860, 520), ("Ы", "ы", Yeru_, 1000, 640),
    ("Ь", "ь", Soft_, 720, 410), ("Э", "э", Ee_, 760, 440), ("Ю", "ю", Yu_, 1100, 680),
    ("Я", "я", Ya_, 740, 410),
]
LAT = [
    ("A", "a", A_, 900, 440), ("B", "b", B_, 760, 410), ("C", "c", Es_, 740, 430),
    ("D", "d", D_lat, 780, 430), ("E", "e", E_, 680, 380), ("F", "f", F_, 640, 360),
    ("G", "g", G_lat, 780, 450), ("H", "h", N_cyr, 780, 440), ("I", "i", I_lat, 120, 120),
    ("J", "j", J_lat, 600, 380), ("K", "k", K_, 780, 430), ("L", "l", L_lat, 600, 350),
    ("M", "m", M_, 980, 580), ("N", "n", N_lat, 800, 440), ("O", "o", O_, 760, 440),
    ("P", "p", Er_, 740, 410), ("Q", "q", Q_, 760, 450), ("R", "r", R_lat, 760, 420),
    ("S", "s", S_, 720, 410), ("T", "t", Te_, 780, 420), ("U", "u", U_lat, 780, 440),
    ("V", "v", V_, 860, 460), ("W", "w", W_, 1200, 700), ("X", "x", Kha_, 860, 460),
    ("Y", "y", Y_lat, 840, 480), ("Z", "z", Z_, 720, 400),
]
LOW_SPECIAL = {"а": a_low, "a": a_low, "е": e_low, "e": e_low}


def with_leaves(g, xs, y=815, k=1.6):
    g = dict(g)
    g["leaves"] = g["leaves"] + [(x, y, 40, k) for x in xs]
    return g


def breve(g, w):
    g = dict(g)
    g["free"] = g["free"] + [arc(w / 2, 900, w * 0.2, 90, 200, 340)]
    return g


def build_letters():
    out = {}   # char -> glyph spec
    for up, low, fn, ww, nw in CYR + LAT:
        out[up] = fn(ww)
        out[low] = LOW_SPECIAL.get(low, fn)(nw)
    # диакритика
    out["Ё"] = with_leaves(E_(680), [210, 480])
    out["ё"] = with_leaves(e_low(380), [105, 290], k=1.3)
    out["Й"] = breve(I_(800), 800)
    out["й"] = breve(I_(440), 440)
    out["i"] = with_leaves(I_lat(120), [60], y=825, k=1.35)
    out["j"] = with_leaves(J_lat(380), [320], y=825, k=1.35)
    return out


# ---------------------------------------------------------------- цифры и знаки

def digits():
    w = 470
    x1 = w - R
    d = {}
    d["0"] = G([oval(R, B, x1, T)])
    d["1"] = G([ln((w * 0.62, 0), (w * 0.62, CH)), ln((R, 470), (w * 0.62, CH))])
    d["2"] = G([P(oval_arc(R, 330, x1, T, 175, -15), ln((x1 - 5, 420), (R, B)), ln((R, B), (w, B)))])
    d["3"] = Ze_(w)
    d["4"] = G([ln((w * 0.68, 0), (w * 0.68, CH)), ln((w * 0.68, CH), (20, 230), (w, 230))])
    d["5"] = G([P(ln((x1 + 20, T), (R + 10, T)), ln((R + 10, T), (R, 390))),
                rbowl(R, 390, B, x1)])
    d["6"] = G([oval(R, B, x1, 430), P(ln((R, 245), (R, T - 120)), arc(R + 120, T - 120, 120, 120, 180, 90),
                                         ln((R + 120, T), (x1, T)))])
    d["7"] = G([ln((30, T), (x1, T), (w * 0.3, 0))])
    d["8"] = G([oval(R + 20, M, x1 - 20, T), oval(R, B, x1, M)])
    six = d["6"]
    d["9"] = G([[(w - x, CH - y) for x, y in s] for s in six["strokes"]])
    return {k: (v, w) for k, v in d.items()}


def punct():
    p = {}
    sq = lambda x, y=0: ln((x, y), (x, y + 2 * R))
    p["."] = G([sq(R)])
    p[","] = G([ln((R, 2 * R), (R, 30), (10, -110))], band=(-200, CH))
    p[":"] = G([sq(R), sq(R, 380)])
    p[";"] = G([sq(R, 380), ln((R, 2 * R), (R, 30), (10, -110))], band=(-200, CH))
    p["!"] = G([ln((R, CH), (R, 250)), sq(R)])
    p["?"] = G([P(oval_arc(R, 300, 420, T, 175, -60), ln((300, 330), (240, 330)), ln((240, 330), (240, 250))),
                sq(240)])
    p["-"] = G([ln((0, M), (300, M))])
    p["–"] = G([ln((0, M), (500, M))])
    p["—"] = G([ln((0, M), (900, M))])
    p["'"] = G([ln((R, CH), (R, 470))])
    p['"'] = G([ln((R, CH), (R, 470)), ln((R + 170, CH), (R + 170, 470))])
    p["’"] = G([ln((R + 20, CH), (R + 20, 580), (R - 20, 470))])
    p["‘"] = G([ln((R - 20, 470), (R - 20, 590), (R + 20, CH))])
    p["«"] = G([ln((200, 560), (R, M), (200, 140)), ln((420, 560), (280, M), (420, 140))])
    p["»"] = G([ln((R, 560), (200, M), (R, 140)), ln((280, 560), (420, M), (280, 140))])
    p["("] = G([oval_arc(R, -150, 460, 850, 100, 260)], band=(-200, 900))
    p[")"] = G([oval_arc(-200, -150, 200, 850, 80, -80)], band=(-200, 900))
    p["["] = G([ln((240, 820), (R, 820), (R, -120), (240, -120))], band=(-200, 900))
    p["]"] = G([ln((0, 820), (180, 820), (180, -120), (0, -120))], band=(-200, 900))
    p["/"] = G([ln((0, -100), (420, 800))], band=(-100, 800))
    p["\\"] = G([ln((0, 800), (420, -100))], band=(-100, 800))
    p["+"] = G([ln((0, M), (440, M)), ln((220, 130), (220, 570))])
    p["="] = G([ln((0, 450), (440, 450)), ln((0, 250), (440, 250))])
    p["*"] = G([ln((170, CH), (170, 420)), ln((40, 640), (300, 480)), ln((40, 480), (300, 640))])
    p["%"] = G([oval(R, 400, 300, T), oval(500, B, 740, 300), ln((60, 0), (740, CH))])
    p["_"] = G([ln((0, -120), (560, -120))], band=(-200, CH))
    p["№"] = G([ln((R, 0), (R, CH)), ln((500, 0), (500, CH)), ln((R + 30, CH), (470, 0)),
                     oval(620, 330, 840, T), ln((610, 180), (850, 180))])
    p["₽"] = G([ln((140, 0), (140, T)), rbowl(140, T, 300, 560), ln((0, 300), (140, 300)),
                     ln((0, 150), (400, 150))])
    p["&"] = G([P(ln((700, 0), (200, 420)), oval_arc(140, 380, 520, T, -60, 210), ln((150, 450), (150, 440))),
                oval_arc(R, B, 520, 430, 110, 360)])
    p["@"] = G([oval(260, 190, 560, 510), P(ln((560, 510), (560, 250)), oval_arc(R, -80, 860, T, -20, 300))],
               band=(-200, CH + 100))
    p["#"] = G([ln((200, 0), (280, CH)), ln((460, 0), (540, CH)), ln((40, 470), (720, 470)), ln((20, 230), (700, 230))])
    p["…"] = G([sq(R), sq(R + 220), sq(R + 440)])
    return p


# ---------------------------------------------------------------- отрисовка

def extend_ends(pts, band):
    """Концы, упирающиеся в линию шрифта, продлеваем — потом их срежет «полоса»."""
    lo, hi = band
    pts = list(pts)
    for idx, nb in ((0, 1), (-1, -2)):
        x, y = pts[idx]
        if abs(y - lo) < 1 or abs(y - hi) < 1 or abs(y) < 1 or abs(y - CH) < 1:
            px, py = pts[nb]
            dx, dy = x - px, y - py
            d = math.hypot(dx, dy)
            if d > 0 and abs(dy) > 1e-6:
                pts[idx] = (x + dx / d * 3 * R, y + dy / d * 3 * R)
    return pts


def stroke(pts, band):
    closed = math.dist(pts[0], pts[-1]) < 1 and len(pts) > 3
    if closed:
        poly = Polygon(pts)
        return poly.buffer(R, quad_segs=24).difference(poly.buffer(-R, quad_segs=24))
    pts = extend_ends(pts, band)
    return LineString(pts).buffer(R, cap_style="flat", join_style="round", quad_segs=24)


def render(spec):
    lo, hi = spec["band"]
    parts = [stroke(s, spec["band"]) for s in spec["strokes"]]
    shape = unary_union(parts)
    shape = shape.intersection(box(-2000, lo, 4000, hi))
    # мягкие углы: и выпуклые, и вогнутые
    shape = shape.buffer(-SOFT, join_style="round").buffer(SOFT, join_style="round")
    shape = shape.buffer(SOFT * 0.6, join_style="round").buffer(-SOFT * 0.6, join_style="round")
    extra = [LineString(s).buffer(R * 0.8, cap_style="flat", join_style="round") for s in spec["free"]]
    extra += [leaf_poly(x, y, a, k) for x, y, a, k in spec["leaves"]]
    if extra:
        shape = unary_union([shape] + extra)
    return shape.simplify(0.6, preserve_topology=True)


def gname(ch):
    names = {"0": "zero", "1": "one", "2": "two", "3": "three", "4": "four",
             "5": "five", "6": "six", "7": "seven", "8": "eight", "9": "nine"}
    if ch.isascii() and ch.isalpha():
        return ch
    return names.get(ch, f"uni{ord(ch):04X}")


# ---------------------------------------------------------------- кернинг

def kerning(shapes, target=2 * SB, tighten=0.5, kmin=-140, kmax=50):
    ys = list(range(-200, 960, 20))
    prof = {}
    for n, (shape, shift, adv) in shapes.items():
        minx, _, maxx, _ = shape.bounds
        left, right = [], []
        for y in ys:
            part = shape.intersection(box(minx - 1, y - 10, maxx + 1, y + 10))
            if part.is_empty:
                left.append(None)
                right.append(None)
            else:
                b = part.bounds
                left.append(b[0] + shift)
                right.append(adv - (b[2] + shift))

        def blur(arr):
            out = []
            for i in range(len(arr)):
                win = [v for v in arr[max(0, i - 2): i + 3] if v is not None]
                out.append(min(win) if win else None)
            return out
        prof[n] = (blur(left), blur(right))
    pairs = {}
    names = list(shapes)
    for a in names:
        for b in names:
            gaps = [x + y for x, y in zip(prof[a][1], prof[b][0]) if x is not None and y is not None]
            if len(gaps) < 2:
                continue
            gap = min(gaps)
            k = -(gap - target) * tighten if gap > target else (target - gap) * 0.8
            k = int(round(max(kmin, min(kmax, k)) / 5.0) * 5)
            if abs(k) >= 15:
                pairs[(a, b)] = k
    lines = ["feature kern {"]
    lines += [f"    pos {a} {b} {k};" for (a, b), k in sorted(pairs.items())]
    lines.append("} kern;")
    print("kerning pairs:", len(pairs))
    return "\n".join(lines)


# ---------------------------------------------------------------- сборка

def main():
    specs = build_letters()
    for k, (v, _) in digits().items():
        specs[k] = v
    specs.update(punct())

    shapes, cmap = {}, {}
    for ch, spec in specs.items():
        name = gname(ch)
        shape = render(spec)
        minx, _, maxx, _ = shape.bounds
        shift = -minx + SB
        adv = int(round(maxx - minx + 2 * SB))
        shapes[name] = (shape, shift, adv)
        cmap[ord(ch)] = name
    cmap[0x20] = "space"
    cmap[0xA0] = "uni00A0"

    kern_fea = kerning(shapes)
    order = [".notdef", "space", "uni00A0"] + list(shapes)

    upper = [gname(u) for u, *_ in CYR + LAT] + ["uni0401", "uni0419"]
    lower = [gname(l) for _, l, *_ in CYR + LAT] + ["uni0451", "uni0439"]
    fea = f"""
languagesystem DFLT dflt;
languagesystem latn dflt;
languagesystem cyrl dflt;
@upper = [{' '.join(upper)}];
@lower = [{' '.join(lower)}];
feature ss01 {{ featureNames {{ name "Wide lowercase"; }}; sub @lower by @upper; }} ss01;
feature ss02 {{ featureNames {{ name "Narrow capitals"; }}; sub @upper by @lower; }} ss02;
{kern_fea}
"""

    for fmt in ("ttf", "otf"):
        fb = FontBuilder(UPM, isTTF=(fmt == "ttf"))
        fb.setupGlyphOrder(order)
        fb.setupCharacterMap(cmap)
        metrics = {".notdef": (500, 50), "space": (280, 0), "uni00A0": (280, 0)}
        glyphs = {}

        def newpen(adv):
            return TTGlyphPen(None) if fmt == "ttf" else T2CharStringPen(adv, None)

        def done(pen):
            return pen.glyph() if fmt == "ttf" else pen.getCharString()

        pen = newpen(500)
        for rect in ([(50, 0), (50, 700), (450, 700), (450, 0)], [(110, 60), (390, 60), (390, 640), (110, 640)]):
            pen.moveTo(rect[0])
            for pt in rect[1:]:
                pen.lineTo(pt)
            pen.closePath()
        glyphs[".notdef"] = done(pen)
        glyphs["space"] = done(newpen(280))
        glyphs["uni00A0"] = done(newpen(280))
        for name, (shape, shift, adv) in shapes.items():
            pen = newpen(adv)
            draw(shape, pen, shift)
            glyphs[name] = done(pen)
            metrics[name] = (adv, int(round(shape.bounds[0] + shift)))

        if fmt == "ttf":
            fb.setupGlyf(glyphs)
        else:
            fb.setupCFF(FAMILY + "-Regular", {"FullName": FAMILY + " Regular"}, glyphs, {})
        fb.setupHorizontalMetrics(metrics)
        fb.setupHorizontalHeader(ascent=950, descent=-250)
        fb.setupNameTable({
            "familyName": FAMILY, "styleName": "Regular",
            "uniqueFontIdentifier": f"{FAMILY}-Regular-1.000",
            "fullName": f"{FAMILY} Regular", "psName": f"{FAMILY}-Regular",
            "version": "Version 1.000", "designer": "Kapelka project",
            "description": "Плакатный геометрический гротеск: широкие прописные, узкие капительные строчные, мягкие углы, листики.",
            "licenseDescription": "Free for personal and commercial use.",
        })
        fb.setupOS2(sTypoAscender=800, sTypoDescender=-200, sTypoLineGap=200,
                    usWinAscent=950, usWinDescent=250, sxHeight=CH, sCapHeight=CH,
                    achVendID="KPLK", fsType=0, usWeightClass=700,
                    ulUnicodeRange1=(1 << 0) | (1 << 1) | (1 << 9), ulCodePageRange1=(1 << 0) | (1 << 2))
        fb.setupPost()
        addOpenTypeFeaturesFromString(fb.font, fea)
        out = os.path.join(HERE, f"{FAMILY}-Regular.{fmt}")
        fb.save(out)
        print("saved", out, len(order), "glyphs")
    specimen()


def specimen():
    from PIL import Image, ImageDraw, ImageFont
    path = os.path.join(HERE, f"{FAMILY}-Regular.otf")

    def F(s, **kw):
        return ImageFont.truetype(path, s, layout_engine=ImageFont.Layout.RAQM)
    img = Image.new("RGB", (1900, 1300), (250, 248, 240))
    d = ImageDraw.Draw(img)
    green, ink = (40, 120, 60), (30, 35, 30)
    y = 30
    d.text((60, y), "ГРЯДКА грядка", font=F(150), fill=green); y += 210
    for line, s in [("АБВГДЕЁЖЗИЙКЛМНОПРСТУФХЦЧШЩЪЫЬЭЮЯ", 50), ("абвгдеёжзийклмнопрстуфхцчшщъыьэюя", 62),
                    ("ABCDEFGHIJKLMNOPQRSTUVWXYZ", 50), ("abcdefghijklmnopqrstuvwxyz", 62),
                    ("0123456789 .,:;!?-–— «»()[]/+=*%№@&#…₽", 58)]:
        d.text((60, y), line, font=F(s), fill=ink); y += 105
    y += 20
    d.text((60, y), "Свежие овощи с грядки!", font=F(100), fill=green); y += 140
    d.text((60, y), "ВКУСНО и полезно", font=F(100), fill=ink); y += 140
    d.text((60, y), "Хлеб «Домашний» — 99,90 ₽", font=F(84), fill=ink)
    img.save(os.path.join(HERE, "gryadka-specimen.png"))
    print("saved specimen")


if __name__ == "__main__":
    main()
