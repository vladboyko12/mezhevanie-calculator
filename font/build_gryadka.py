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

from shapely.geometry.polygon import orient
from fontTools.pens.cu2quPen import Cu2QuPen

from build_kapelka import ln, arc, bz, P, leaf_poly, polys_of
from bezier_fit import draw_ring

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


def rbowl(xl, ytop, ybot, xr, xb=None):
    """Правая чаша: от (xl,ytop) вправо, скруглённо вниз, обратно к (xb,ybot)."""
    xb = xl if xb is None else xb
    r = min((ytop - ybot) / 2, xr - max(xl, xb))
    return P(ln((xl, ytop), (xr - r, ytop)),
             arc(xr - r, ytop - r, r, r, 90, 0),
             ln((xr, ytop - r), (xr, ybot + r)),
             arc(xr - r, ybot + r, r, r, 0, -90),
             ln((xr - r, ybot), (xb, ybot)))


def lbowl(xr, ytop, ybot, xl):
    """Левая чаша (зеркало rbowl)."""
    r = min((ytop - ybot) / 2, xr - xl)
    return P(ln((xr, ytop), (xl + r, ytop)),
             arc(xl + r, ytop - r, r, r, 90, 180),
             ln((xl, ytop - r), (xl, ybot + r)),
             arc(xl + r, ybot + r, r, r, 180, 270),
             ln((xl + r, ybot), (xr, ybot)))


def G(strokes, free=(), leaves=(), band=(0, CH), geoms=()):
    return {"strokes": [list(s) for s in strokes], "free": [list(s) for s in free],
            "leaves": list(leaves), "band": band, "geoms": list(geoms)}


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


def soft_(w, x0=None):
    x0 = R if x0 is None else x0
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


def cw(w):
    """Компенсация ширины: у жирного начертания буквы шире, чтобы не забивались просветы."""
    return w + (R - 60) * (1.3 if w >= 600 else 1.9)


def with_leaves(g, xs, y=815, k=1.6):
    g = dict(g)
    k *= (R / 60) ** 0.45
    g["leaves"] = g["leaves"] + [(x, y, 40, k) for x in xs]
    return g


def breve(g, w):
    g = dict(g)
    g["free"] = g["free"] + [arc(w / 2, 900, w * 0.2, 90, 200, 340)]
    return g


def build_letters():
    out = {}   # char -> glyph spec
    for up, low, fn, ww, nw in CYR + LAT:
        out[up] = fn(cw(ww))
        out[low] = LOW_SPECIAL.get(low, fn)(cw(nw))
    # диакритика
    we, ne = cw(680), cw(380)
    out["Ё"] = with_leaves(E_(we), [we * 0.31, we * 0.7])
    out["ё"] = with_leaves(e_low(ne), [ne * 0.28, ne * 0.76], k=1.3)
    out["Й"] = breve(I_(cw(800)), cw(800))
    out["й"] = breve(I_(cw(440)), cw(440))
    out["i"] = with_leaves(I_lat(120), [R], y=825, k=1.35)
    out["j"] = with_leaves(J_lat(cw(380)), [cw(380) - R], y=825, k=1.35)
    return out


# ---------------------------------------------------------------- цифры и знаки

def digits():
    """Табличные цифры: у всех одинаковая ширина — удобно для ценников."""
    w = cw(480)
    x1 = w - R
    d = {}
    d["0"] = G([oval(R, B, x1, T)])
    xs = w * 0.62
    d["1"] = G([ln((xs, 0), (xs, CH)), ln((R + 10, 480), (xs, CH))])
    top2 = oval_arc(R, 300, x1, T, 168, -40)
    d["2"] = G([P(top2, ln(top2[-1], (R, B)), ln((R, B), (w, B)))])
    d["3"] = G([P(ln((30, T), (x1, T), (w * 0.4, 430)), rbowl(w * 0.4, 430, B, x1, xb=30))])
    d["4"] = G([ln((w * 0.68, 0), (w * 0.68, CH)), ln((w * 0.68, CH), (20, 230), (w, 230))])
    d["5"] = G([P(ln((x1 + 20, T), (R, T), (R, 410)), rbowl(R, 410, B, x1, xb=30))])
    r6 = 150
    d["6"] = G([oval(R, B, x1, 450), P(ln((R, 255), (R, T - r6)), arc(R + r6, T - r6, r6, r6, 180, 90),
                                         ln((R + r6, T), (x1, T)))])
    d["7"] = G([ln((30, T), (x1, T), (w * 0.32, 0))])
    d["8"] = G([oval(R + 25, M, x1 - 25, T), oval(R, B, x1, M)])
    d["9"] = G([[(w - x, CH - y) for x, y in s] for s in d["6"]["strokes"]])
    return d, w


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


# ---------------------------------------------------------------- расширенный набор знаков
import unicodedata


def rot180(g, w):
    """Повернуть знак на 180° внутри ширины w (для Ә, Һ, ¡, ¿)."""
    t = lambda s: [(w - x, CH - y) for x, y in s]
    out = dict(g)
    out["strokes"] = [t(s) for s in g["strokes"]]
    out["free"] = [t(s) for s in g["free"]]
    lo, hi = g["band"]
    out["band"] = (CH - hi, CH - lo)
    return out


def addg(g, strokes=(), geoms=(), leaves=(), band=None):
    out = dict(g)
    out["strokes"] = g["strokes"] + [list(s) for s in strokes]
    out["geoms"] = g.get("geoms", []) + list(geoms)
    out["leaves"] = g["leaves"] + list(leaves)
    if band:
        out["band"] = band
    return out


def mstroke(pts, k=0.78):
    return LineString(pts).buffer(R * k, cap_style="flat", join_style="round", quad_segs=16)


def skel_box(g):
    xs = [x for s in g["strokes"] for x, _ in s]
    return min(xs), max(xs)


# ---- диакритические знаки: функции (cx, xr) -> shapely-геометрия
MK = 0.55       # толщина значков относительно основного штриха
MH = 125        # высота зоны значков (по осевой)


def _yb():
    return CH + 40 + R * MK


def _above(geom):
    return geom.intersection(box(-2000, CH + 40, 4000, _yb() + MH + R * MK))


def m_acute(cx, xr):
    yb = _yb()
    return _above(mstroke([(cx - 70, yb - 30), (cx + 80, yb + MH + 30)], MK))


def m_grave(cx, xr):
    yb = _yb()
    return _above(mstroke([(cx + 70, yb - 30), (cx - 80, yb + MH + 30)], MK))


def m_circ(cx, xr):
    yb = _yb()
    return _above(mstroke([(cx - 110, yb - 20), (cx, yb + MH - R * MK * 0.4), (cx + 110, yb - 20)], MK))


def m_caron(cx, xr):
    yb = _yb()
    return _above(mstroke([(cx - 110, yb + MH + 20), (cx, yb + R * MK * 0.4), (cx + 110, yb + MH + 20)], MK))


def m_breve(cx, xr):
    yb = _yb()
    return _above(mstroke(arc(cx, yb + MH, 100, MH - 10, 192, 348), MK))


def m_tilde(cx, xr):
    yb = _yb()
    return _above(mstroke(bz((cx - 115, yb + 15), (cx - 50, yb + MH + 40), (cx + 50, yb - 40), (cx + 115, yb + MH - 15)), MK))


def m_macron(cx, xr):
    yb = _yb()
    return _above(mstroke([(cx - 115, yb + MH / 2), (cx + 115, yb + MH / 2)], MK))


def m_ring(cx, xr):
    yb = _yb()
    ro = MH / 2 + R * MK
    ri = max(22, ro - 2 * R * MK)
    c = Point(cx, yb + MH / 2)
    return c.buffer(ro, quad_segs=24).difference(c.buffer(ri, quad_segs=24))


def m_dot(cx, xr):
    return leaf_poly(cx, _yb() + MH / 2, 40, 1.2 * (R / 60) ** 0.45)


def m_diaer(cx, xr):
    k = 1.15 * (R / 60) ** 0.45
    d = 75 + R * 0.5
    return unary_union([leaf_poly(cx - d, _yb() + MH / 2, 40, k), leaf_poly(cx + d, _yb() + MH / 2, 40, k)])


def m_dacute(cx, xr):
    yb = _yb()
    return _above(unary_union([mstroke([(cx - 120, yb - 30), (cx - 30, yb + MH + 30)], MK * 0.9),
                               mstroke([(cx + 30, yb - 30), (cx + 120, yb + MH + 30)], MK * 0.9)]))


def _below(geom):
    return geom.intersection(box(-2000, -280, 4000, R))


def m_cedilla(cx, xr):
    return _below(mstroke([(cx, R), (cx, -70), (cx - 80, -200)], MK * 1.1))


def m_comma_below(cx, xr):
    return _below(mstroke([(cx + 20, -60), (cx + 20, -120), (cx - 35, -230)], MK * 1.2))


def m_ogonek(cx, xr):
    x = xr - R * 0.3
    return _below(mstroke(bz((x, R), (x - 100, -50), (x - 80, -190), (x + 35, -185)), MK * 1.1))


MARKS = {0x0300: m_grave, 0x0301: m_acute, 0x0302: m_circ, 0x0303: m_tilde, 0x0304: m_macron,
         0x0306: m_breve, 0x0307: m_dot, 0x0308: m_diaer, 0x030A: m_ring, 0x030B: m_dacute,
         0x030C: m_caron, 0x0326: m_comma_below, 0x0327: m_cedilla, 0x0328: m_ogonek}
APOS_CARON = {"ď", "ľ", "Ľ", "ť"}     # словацко-чешский «карон-апостроф»
DOTLESS = {"i": "ı", "j": "ȷ", "і": "ı", "ј": "ȷ"}


def base_table():
    t = {}
    for up, low, fn, ww, nw in CYR + LAT:
        t[up] = (fn, cw(ww))
        t[low] = (LOW_SPECIAL.get(low, fn), cw(nw))
    t["ı"] = (I_lat, 120)
    t["ȷ"] = (J_lat, cw(380))
    t["І"] = (I_lat, 120)
    return t


def specials():
    """Буквы, которые не собираются из основы и знака."""
    sp = {}

    def both(up, low, fn, ww, nw):
        sp[up] = fn(cw(ww))
        sp[low] = fn(cw(nw))

    def AE(w):
        xm = w * 0.46
        return G([ln((0, 0), (xm, CH)), P(ln((w, T), (xm, T)), ln((xm, T), (xm, B)), ln((xm, B), (w, B))),
                  ln((xm, M), (w - 50, M)), ln((xm * 0.3, 230), (xm, 230))])

    def OE(w):
        xm = w * 0.5
        return G([oval(R, B, xm, T), ln((xm, 0), (xm, CH)),
                  ln((xm, T), (w, T)), ln((xm, M), (w - 50, M)), ln((xm, B), (w, B))])

    def Oslash(w):
        return addg(O_(w), strokes=[ln((R * 0.4, 25), (w - R * 0.4, CH - 25))])

    def Dbar(w):
        return addg(D_lat(w), strokes=[ln((R - 130, M), (R + 140, M))])

    def Thorn(w):
        return G([ln((R, 0), (R, CH)), rbowl(R, 560, 150, w - R)])

    def Hbar(w):
        return addg(N_cyr(w), strokes=[ln((R - 90, 545), (w - R + 90, 545))])

    def Lstroke(w):
        return addg(L_lat(w), strokes=[ln((R - 110, 250), (R + 130, 450))])

    def Ldot(w):
        xd = R + (w - R) * 0.55
        return addg(L_lat(w), strokes=[ln((xd, M - R), (xd, M + R))])

    def Eng(w):
        x1 = w - R
        return G([P(ln((R, 0), (R, T)), ln((R, T), (x1, T)), ln((x1, T), (x1, -40)),
                    arc(x1 - 110, -40, 110, 100, 0, -90), ln((x1 - 110, -140), (x1 - 200, -140)))],
                 band=(-300, CH))

    def Tbar(w):
        return addg(Te_(w), strokes=[ln((w / 2 - 140, M), (w / 2 + 140, M))])

    def IJ(w, dots=False):
        j = J_lat(w)
        sh = 2 * R + 110
        g = G([ln((R, 0), (R, CH))] + [[(x + sh, y) for x, y in s] for s in j["strokes"]])
        if dots:
            g = with_leaves(g, [R, sh + w - R], y=825, k=1.2)
        return g

    def Sharp(w):
        x1 = w - R
        return G([P(ln((R, 0), (R, T)), ln((R, T), (x1, T)), ln((x1, T), (w * 0.42, 430))),
                  rbowl(w * 0.42, 430, B, x1, xb=w * 0.3)])

    for up, low, fn, ww, nw in [("Æ", "æ", AE, 1040, 640), ("Œ", "œ", OE, 1100, 680),
                                ("Ø", "ø", Oslash, 760, 440), ("Ð", "ð", Dbar, 780, 430),
                                ("Đ", "đ", Dbar, 780, 430), ("Þ", "þ", Thorn, 700, 400),
                                ("Ħ", "ħ", Hbar, 780, 440), ("Ł", "ł", Lstroke, 600, 350),
                                ("Ŀ", "ŀ", Ldot, 600, 350), ("Ŋ", "ŋ", Eng, 780, 440),
                                ("Ŧ", "ŧ", Tbar, 780, 420), ("ẞ", "ß", Sharp, 700, 430)]:
        both(up, low, fn, ww, nw)
    sp["Ĳ"] = IJ(cw(600))
    sp["ĳ"] = IJ(cw(380), dots=True)
    sp["ı"] = I_lat(120)
    sp["ĸ"] = K_(cw(430))
    sp["İ"] = with_leaves(I_lat(120), [R], y=825, k=1.35)

    # кириллица: украинский, белорусский, казахский, сербский, македонский
    def Ye_ukr(w):
        return G([oval_arc(R, B, w - R, T, 38, 322), ln((R, M), (w * 0.62, M))])

    def Ghe_up(w):
        x1 = w - R
        return G([P(ln((R, 0), (R, T)), ln((R, T), (x1, T)), ln((x1, T), (x1, CH + 150)))], band=(0, CH + 150))

    def Schwa(w):
        return rot180(e_low(w), w)

    def Ghe_bar(w):
        return addg(Ge_(w), strokes=[ln((R - 100, 330), (R + 180, 330))])

    def Ka_tail(w):
        return addg(K_(w), strokes=[ln((w - 50, B), (w + 40, B), (w + 40, -140))], band=(-140, CH))

    def En_tail(w):
        x1 = w - R
        return addg(N_cyr(w), strokes=[ln((x1, B), (x1 + 75, B), (x1 + 75, -140))], band=(-140, CH))

    def O_bar(w):
        return addg(O_(w), strokes=[ln((R, M), (w - R, M))])

    def U_str_bar(w):
        return addg(Y_lat(w), strokes=[ln((w / 2 - 140, 230), (w / 2 + 140, 230))])

    def Shha(w):
        return rot180(Che_(w), w)

    def Dzhe(w):
        x1 = w - R
        return G([ln((R, CH), (R, B), (x1, B), (x1, CH)), ln((w / 2, B), (w / 2, -140))], band=(-140, CH))

    def Lje(w):
        w1 = w * 0.56
        return addg(El_(w1), strokes=[rbowl(w1 - R, 410, B, w - R)])

    def Nje(w):
        w1 = w * 0.55
        return addg(N_cyr(w1), strokes=[rbowl(w1 - R, 410, B, w - R)])

    def Tshe(w, hook=False):
        x1, xs = w - R, w * 0.28
        r = 140
        leg = [ln((xs + 5, 430), (x1 - r, 430)), arc(x1 - r, 430 - r, r, r, 90, 0)]
        if hook:
            leg += [ln((x1, 430 - r), (x1, -30)), arc(x1 - 100, -30, 100, 100, 0, -90),
                    ln((x1 - 100, -130), (x1 - 190, -130))]
        else:
            leg += [ln((x1, 430 - r), (x1, 0))]
        return G([ln((0, T), (w * 0.62, T)), ln((xs, 0), (xs, CH)), P(*leg)],
                 band=(-260, CH) if hook else (0, CH))

    for up, low, fn, ww, nw in [("Є", "є", Ye_ukr, 740, 430), ("Ґ", "ґ", Ghe_up, 640, 360),
                                ("Ә", "ә", Schwa, 760, 420), ("Ғ", "ғ", Ghe_bar, 640, 360),
                                ("Қ", "қ", Ka_tail, 780, 430), ("Ң", "ң", En_tail, 780, 440),
                                ("Ө", "ө", O_bar, 760, 440), ("Ү", "ү", Y_lat, 840, 480),
                                ("Ұ", "ұ", U_str_bar, 840, 480), ("Һ", "һ", Shha, 720, 420),
                                ("Ј", "ј", J_lat, 600, 380), ("Ѕ", "ѕ", S_, 720, 410),
                                ("Џ", "џ", Dzhe, 780, 440), ("Љ", "љ", Lje, 1120, 680),
                                ("Њ", "њ", Nje, 1080, 660), ("Ћ", "ћ", Tshe, 820, 480)]:
        both(up, low, fn, ww, nw)
    sp["Ђ"] = Tshe(cw(820), hook=True)
    sp["ђ"] = Tshe(cw(480), hook=True)
    sp["І"] = I_lat(120)
    sp["і"] = with_leaves(I_lat(120), [R], y=825, k=1.35)
    sp["ј"] = with_leaves(J_lat(cw(380)), [cw(380) - R], y=825, k=1.35)
    dx = 85 + R * 0.55
    sp["Ї"] = with_leaves(I_lat(120), [R - dx, R + dx], y=825, k=0.95)
    sp["ї"] = with_leaves(I_lat(120), [R - dx, R + dx], y=825, k=0.95)
    return sp


EXT_RANGES = list(range(0x00C0, 0x0100)) + list(range(0x0100, 0x0180)) + [0x0218, 0x0219, 0x021A, 0x021B] \
    + list(range(0x0400, 0x0460)) + [0x0490, 0x0491, 0x0492, 0x0493, 0x049A, 0x049B, 0x04A2, 0x04A3,
                                     0x04AE, 0x04AF, 0x04B0, 0x04B1, 0x04BA, 0x04BB, 0x04D8, 0x04D9,
                                     0x04E8, 0x04E9, 0x1E9E]
SKIP = {"×", "÷", "ſ", "ŉ"}


KERN_BASE = {}   # буква с диакритикой -> основа (для классового кернинга)


def extended(existing):
    """Расширенная латиница и кириллица: специальные формы + основа со знаком."""
    out = {}
    sp = specials()
    bases = base_table()
    for cp in EXT_RANGES:
        ch = chr(cp)
        if ch in existing or ch in SKIP or not ch.isalpha():
            continue
        if ch in sp:
            out[ch] = sp[ch]
            continue
        dec = unicodedata.decomposition(ch).split()
        if len(dec) != 2 or dec[0].startswith("<"):
            continue
        base, mark = chr(int(dec[0], 16)), int(dec[1], 16)
        base = DOTLESS.get(base, base)
        if base not in bases or mark not in MARKS:
            continue
        fn, w = bases[base]
        g = fn(w)
        lo, hi = skel_box(g)
        cx = (lo + hi) / 2
        if ch in APOS_CARON:
            xa = hi + 90 if ch in ("ď", "ť") else R + 150
            geom = mstroke([(xa + 15, CH), (xa + 15, CH - 90), (xa - 20, CH - 190)], 0.8)
            geom = geom.intersection(box(-2000, 0, 4000, CH))
        else:
            geom = MARKS[mark](cx, hi)
        out[ch] = addg(g, geoms=[geom])
        KERN_BASE[ch] = base if base in existing else dec_base(ch)
    return out


def dec_base(ch):
    return chr(int(unicodedata.decomposition(ch).split()[0], 16))


def punct_extra():
    p = {}
    sq = lambda x, y=0: ln((x, y), (x, y + 2 * R))
    comma = lambda x, y=0: ln((x, y + 2 * R), (x, y + 30), (x - 50, y - 110))
    turned = lambda x, y: ln((x - 20, y), (x - 20, y + 110), (x + 20, y + 230))
    # улучшенные «?», «&», «@»
    w = cw(440)
    x1, xc = w - R, w / 2
    p["?"] = G([P(oval_arc(R, 330, x1, T, 165, -90), ln((xc, 330), (xc, 250))), sq(xc)])
    p["¿"] = rot180(p["?"], w)
    p["¡"] = rot180(G([ln((R, CH), (R, 250)), sq(R)]), 2 * R)
    w = cw(720)
    p["&"] = G([oval(w * 0.14, 410, w * 0.56, T), ln((w * 0.22, 445), (w, 0)),
                P(oval_arc(R, B, w * 0.72, 480, 118, 360), ln((w * 0.72, 270), (w, 270)))])
    w = cw(880)
    p["@"] = G([oval(w * 0.33, 170, w * 0.63, 470),
                P(ln((w * 0.63, 480), (w * 0.63, 170)), oval_arc(R, -100, w - R, CH + 40, -32, 300))],
               band=(-220, 900))
    # кавычки
    p["„"] = G([comma(R + 50), comma(R + 230)], band=(-250, CH))
    p["‚"] = G([comma(R + 50)], band=(-250, CH))
    p["“"] = G([turned(R + 20, 470), turned(R + 200, 470)])
    p["”"] = G([ln((R + 20, CH), (R + 20, 590), (R - 20, 470)), ln((R + 200, CH), (R + 200, 590), (R + 160, 470))])
    p["‹"] = G([ln((200, 560), (R, M), (200, 140))])
    p["›"] = G([ln((R, 560), (200, M), (R, 140))])
    # скобки и математика
    p["{"] = G([ln((230, 820), (140, 820), (140, 420), (R, M), (140, 280), (140, -120), (230, -120))], band=(-200, 900))
    p["}"] = G([ln((0, 820), (90, 820), (90, 420), (170, M), (90, 280), (90, -120), (0, -120))], band=(-200, 900))
    p["<"] = G([ln((440, 600), (R, M), (440, 100))])
    p[">"] = G([ln((0, 600), (380, M), (0, 100))])
    p["≤"] = G([ln((440, 660), (R, 450), (440, 240)), ln((R, B), (440, B))])
    p["≥"] = G([ln((0, 660), (380, 450), (0, 240)), ln((0, B), (380, B))])
    p["±"] = G([ln((0, 430), (440, 430)), ln((220, 210), (220, 650)), ln((0, B), (440, B))])
    p["×"] = G([ln((40, 140), (400, 560)), ln((40, 560), (400, 140))])
    p["÷"] = G([ln((0, M), (440, M)), sq(220, 490), sq(220, 90)])
    p["−"] = G([ln((0, M), (440, M))])
    p["~"] = G([bz((0, 300), (110, 470), (330, 230), (440, 400))])
    p["≈"] = G([bz((0, 420), (110, 590), (330, 350), (440, 520)), bz((0, 180), (110, 350), (330, 110), (440, 280))])
    p["≠"] = G([ln((0, 450), (440, 450)), ln((0, 250), (440, 250)), ln((110, 60), (330, 640))])
    p["^"] = G([ln((0, 450), (200, CH), (400, 450))])
    p["`"] = G([ln((150, 640), (20, 780))], band=(0, 900))
    p["|"] = G([ln((R, -150), (R, 850))], band=(-150, 850))
    p["¦"] = G([ln((R, -150), (R, 250)), ln((R, 450), (R, 850))], band=(-150, 850))
    p["′"] = G([ln((70, CH), (20, 480))])
    p["″"] = G([ln((70, CH), (20, 480)), ln((230, CH), (180, 480))])
    p["·"] = G([sq(R, M - R)])
    p["•"] = G([], geoms=[Point(150, M).buffer(max(95, R * 1.4), quad_segs=24)])
    # валюты
    w = cw(420)
    p["$"] = addg(S_(w), strokes=[ln((w / 2, -100), (w / 2, 800))], band=(-100, 800))
    w = cw(480)
    p["€"] = G([oval_arc(R, B, w - R, T, 40, 320), ln((-40, 410), (w * 0.6, 410)), ln((-40, 290), (w * 0.6, 290))])
    p["£"] = G([P(oval_arc(w * 0.25, 330, w - R, T, 15, 180), ln((w * 0.25, (330 + T) / 2), (w * 0.25, B))),
                ln((0, B), (w, B)), ln((0, M), (w * 0.62, M))])
    w = cw(700)
    p["¥"] = addg(Y_lat(w), strokes=[ln((w / 2 - 170, 300), (w / 2 + 170, 300)), ln((w / 2 - 170, 150), (w / 2 + 170, 150))])
    w = cw(430)
    p["¢"] = G([oval_arc(R, 110, w - R, 590, 40, 320), ln((w / 2, 0), (w / 2, CH))])
    # прочее
    ro = max(110, R * 1.4 + 40)
    p["°"] = G([], geoms=[Point(ro, 700 - ro).buffer(ro, quad_segs=24).difference(
        Point(ro, 700 - ro).buffer(ro - R * 1.4, quad_segs=24))])
    w = 760
    thin = 0.55
    ring = oval(R, -30, w - R, 730)
    p["©"] = G([ring], band=(-200, 900), geoms=[mstroke(oval_arc(230, 190, 530, 510, 45, 315), thin)])
    p["®"] = G([ring], band=(-200, 900), geoms=[
        mstroke(ln((290, 170), (290, 530)), thin),
        mstroke(rbowl(290, 530, 360, 500), thin),
        mstroke(ln((400, 360), (510, 170)), thin)])
    p["™"] = G([], geoms=[mstroke(ln((0, 680), (280, 680)), thin), mstroke(ln((140, 680), (140, 420)), thin),
                          mstroke(ln((360, 420), (360, 680), (470, 500), (580, 680), (580, 420)), thin)])
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
    extra += spec.get("geoms", [])
    if extra:
        shape = unary_union([shape] + extra)
    return shape


def gname(ch):
    names = {"0": "zero", "1": "one", "2": "two", "3": "three", "4": "four",
             "5": "five", "6": "six", "7": "seven", "8": "eight", "9": "nine"}
    if ch.isascii() and ch.isalpha():
        return ch
    return names.get(ch, f"uni{ord(ch):04X}")


# ---------------------------------------------------------------- кернинг

def kerning(shapes, groups, target, tighten=0.5, kmin=-140, kmax=50):
    """Автокернинг по профилям. shapes — представители классов, groups — {представитель: [члены]}."""
    ys = list(range(-260, 1000, 20))
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
        ra = prof[a][1]
        for b in names:
            gaps = [x + y for x, y in zip(ra, prof[b][0]) if x is not None and y is not None]
            if len(gaps) < 2:
                continue
            gap = min(gaps)
            k = -(gap - target) * tighten if gap > target else (target - gap) * 0.8
            k = int(round(max(kmin, min(kmax, k)) / 5.0) * 5)
            if abs(k) >= 15:
                pairs[(a, b)] = k
    lines = [f"@k_{n} = [{' '.join(groups[n])}];" for n in names]
    lines.append("feature kern {")
    lines += [f"    pos @k_{a} @k_{b} {k};" for (a, b), k in sorted(pairs.items())]
    lines.append("} kern;")
    print("kerning pairs:", len(pairs), "classes:", len(names))
    return "\n".join(lines)


# ---------------------------------------------------------------- сборка

WEIGHTS = [
    # стиль,    полуширина штриха, мягкие углы, отступ, usWeightClass
    ("Light",   32, 10, 44, 300),
    ("Regular", 60, 16, 48, 400),
    ("Black",   88, 20, 52, 900),
]


def set_weight(r, soft, sb):
    global R, T, B, SOFT, SB
    R, SOFT, SB = r, soft, sb
    T, B = CH - R, R


def draw_curves(shape, pen, shift, fmt):
    """Контуры — кривыми Безье. TTF: внешний контур по часовой, OTF — против."""
    sign = -1.0 if fmt == "ttf" else 1.0
    for poly in polys_of(shape):
        poly = orient(poly, sign=sign)
        for ring in [poly.exterior] + list(poly.interiors):
            draw_ring(list(ring.coords)[:-1], pen, shift)


def build_weight(style, r, soft, sb, wclass):
    set_weight(r, soft, sb)
    specs = build_letters()
    dig, dw = digits()
    specs.update(punct())
    specs.update(punct_extra())
    KERN_BASE.clear()
    specs.update(extended(specs))

    shapes, cmap = {}, {}
    for ch, spec in list(specs.items()) + list(dig.items()):
        name = gname(ch)
        shape = render(spec)
        minx, _, maxx, _ = shape.bounds
        if ch in dig:                     # табличные цифры: общая ширина, знак по центру
            adv = int(round(dw + 2 * SB))
            shift = (adv - (maxx - minx)) / 2 - minx
        else:
            shift = -minx + SB
            adv = int(round(maxx - minx + 2 * SB))
        shapes[name] = (shape, shift, adv)
        cmap[ord(ch)] = name
    cmap[0x20] = "space"
    cmap[0xA0] = "uni00A0"
    space = int(round(240 + 1.2 * R))

    groups = {}
    for ch in specs:
        rep = KERN_BASE.get(ch, ch)
        groups.setdefault(gname(rep), []).append(gname(ch))
    reps = {n: shapes[n] for n in groups}
    kern_fea = kerning(reps, groups, target=2 * SB)
    order = [".notdef", "space", "uni00A0"] + list(shapes)

    upper, lower = [], []
    for ch in specs:
        lo = ch.lower()
        if ch != lo and len(lo) == 1 and lo in specs:
            upper.append(gname(ch))
            lower.append(gname(lo))
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
    ps = f"{FAMILY}-{style}"
    for fmt in ("ttf", "otf"):
        fb = FontBuilder(UPM, isTTF=(fmt == "ttf"))
        fb.setupGlyphOrder(order)
        fb.setupCharacterMap(cmap)
        metrics = {".notdef": (500, 50), "space": (space, 0), "uni00A0": (space, 0)}
        glyphs = {}

        def newpen(adv):
            if fmt == "ttf":
                tt = TTGlyphPen(None)
                return tt, Cu2QuPen(tt, max_err=0.6, reverse_direction=False)
            t2 = T2CharStringPen(adv, None)
            return t2, t2

        def done(base):
            return base.glyph() if fmt == "ttf" else base.getCharString()

        base, pen = newpen(500)
        for rect in ([(50, 0), (50, 700), (450, 700), (450, 0)], [(110, 60), (390, 60), (390, 640), (110, 640)]):
            if fmt == "otf":
                rect = rect[::-1]
            pen.moveTo(rect[0])
            for pt in rect[1:]:
                pen.lineTo(pt)
            pen.closePath()
        glyphs[".notdef"] = done(base)
        for sp in ("space", "uni00A0"):
            base, _ = newpen(space)
            glyphs[sp] = done(base)
        for name, (shape, shift, adv) in shapes.items():
            base, pen = newpen(adv)
            draw_curves(shape, pen, shift, fmt)
            glyphs[name] = done(base)
            metrics[name] = (adv, int(round(shape.bounds[0] + shift)))

        if fmt == "ttf":
            fb.setupGlyf(glyphs)
            glyf = fb.font["glyf"]
            for name in shapes:
                g = glyf[name]
                g.recalcBounds(glyf)
                metrics[name] = (metrics[name][0], g.xMin)
        else:
            fb.setupCFF(ps, {"FullName": f"{FAMILY} {style}"}, glyphs, {})
        fb.setupHorizontalMetrics(metrics)
        fb.setupHorizontalHeader(ascent=1000, descent=-300)
        names = {
            "familyName": FAMILY if style == "Regular" else f"{FAMILY} {style}",
            "styleName": "Regular",
            "typographicFamily": FAMILY, "typographicSubfamily": style,
            "uniqueFontIdentifier": f"{ps}-2.000",
            "fullName": f"{FAMILY} {style}", "psName": ps,
            "version": "Version 2.000", "designer": "Kapelka project",
            "description": "Плакатный геометрический гротеск: широкие прописные, узкие капительные строчные, мягкие углы, листики.",
            "licenseDescription": "Free for personal and commercial use.",
        }
        fb.setupNameTable(names)
        fb.setupOS2(sTypoAscender=800, sTypoDescender=-200, sTypoLineGap=300,
                    usWinAscent=1020, usWinDescent=300, sxHeight=CH, sCapHeight=CH,
                    achVendID="KPLK", fsType=0, usWeightClass=wclass,
                    version=4, fsSelection=(0x40 if style == "Regular" else 0) | 0x80,
                    ulUnicodeRange1=(1 << 0) | (1 << 1) | (1 << 9), ulCodePageRange1=(1 << 0) | (1 << 2))
        fb.setupPost()
        addOpenTypeFeaturesFromString(fb.font, fea)
        out = os.path.join(HERE, f"{ps}.{fmt}")
        fb.save(out)
        print("saved", out, len(order), "glyphs")


def main():
    for w in WEIGHTS:
        build_weight(*w)
    specimen()


def specimen():
    from PIL import Image, ImageDraw, ImageFont

    def F(style, s):
        return ImageFont.truetype(os.path.join(HERE, f"{FAMILY}-{style}.otf"), s,
                                  layout_engine=ImageFont.Layout.RAQM)
    img = Image.new("RGB", (1900, 1760), (250, 248, 240))
    d = ImageDraw.Draw(img)
    green, ink, grey = (40, 120, 60), (30, 35, 30), (120, 120, 110)
    lab = F("Regular", 34)
    y = 30
    for style, _, *__ in WEIGHTS:
        d.text((60, y + 30), style, font=lab, fill=grey)
        d.text((300, y), "ГРЯДКА грядка 2025", font=F(style, 120), fill=green if style == "Black" else ink)
        y += 170
    y += 10
    for style, _, *__ in WEIGHTS:
        d.text((60, y + 10), style, font=lab, fill=grey)
        d.text((300, y), "АБВГДЕЁЖЗИЙКЛМНОПРСТУФХЦЧШЩЪЫЬЭЮЯ", font=F(style, 40), fill=ink)
        d.text((300, y + 70), "абвгдеёжзийклмнопрстуфхцчшщъыьэюя  0123456789", font=F(style, 50), fill=ink)
        y += 170
    y += 10
    d.text((60, y), "Цифры", font=lab, fill=grey)
    for i, style in enumerate(("Light", "Regular", "Black")):
        d.text((300, y + i * 95), "0123456789  99,90 ₽", font=F(style, 80), fill=ink)
    y += 300
    d.text((60, y), "СВЕЖИЕ ОВОЩИ", font=F("Black", 110), fill=green)
    d.text((60, y + 140), "с грядки — каждый день", font=F("Light", 80), fill=ink)
    d.text((60, y + 260), "Zażółć gęślą jaźń · Ünïcödé · Їжак · Қазақ · Ђурђевак", font=F("Regular", 58), fill=ink)
    img.save(os.path.join(HERE, "gryadka-specimen.png"))
    print("saved specimen")


if __name__ == "__main__":
    main()
