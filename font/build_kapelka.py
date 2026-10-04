#!/usr/bin/env python3
"""
Kapelka — рукописный «фломастерный» дисплейный шрифт в духе Villula
(монолинейный штрих, скруглённые концы, буквы разной ширины),
но со своими изюминками:

  1. Точки над i, j, ё, Ё и т.п. — маленькие листики.
  2. «Капли чернил» — концы штрихов чуть набухают, как у фломастера.
  3. Живой штрих — толщина и траектория слегка «дышат».
  4. Прыгающая строка — фича calt автоматически чередует три варианта
     каждой буквы (обычный, узкий-приподнятый, широкий-опущенный).

Сборка:  pip install fonttools shapely pillow && python3 build_kapelka.py
Результат: Kapelka-Regular.ttf, Kapelka-Regular.otf, specimen.png
"""
import math
import os
import hashlib

from shapely.geometry import Point, Polygon, MultiPolygon
from shapely.ops import unary_union
from shapely.geometry.polygon import orient
from shapely import affinity

from fontTools.fontBuilder import FontBuilder
from fontTools.pens.ttGlyphPen import TTGlyphPen
from fontTools.pens.t2CharStringPen import T2CharStringPen
from fontTools.feaLib.builder import addOpenTypeFeaturesFromString

HERE = os.path.dirname(os.path.abspath(__file__))
FAMILY = "Kapelka"

UPM = 1000
XH = 500        # высота строчных
CH = 700        # высота прописных
ASC = 740
DSC = -210
R = 54          # половина толщины штриха
SB = 42         # боковые отступы

# ---------------------------------------------------------------- геометрия


def ln(*pts):
    return [tuple(map(float, p)) for p in pts]


def bz(p0, p1, p2, p3, n=48):
    out = []
    for i in range(n + 1):
        t = i / n
        a = (1 - t) ** 3
        b = 3 * (1 - t) ** 2 * t
        c = 3 * (1 - t) * t ** 2
        d = t ** 3
        out.append((a * p0[0] + b * p1[0] + c * p2[0] + d * p3[0],
                    a * p0[1] + b * p1[1] + c * p2[1] + d * p3[1]))
    return out


def arc(cx, cy, rx, ry, a0, a1):
    n = max(12, int(abs(a1 - a0) / 3))
    out = []
    for i in range(n + 1):
        a = math.radians(a0 + (a1 - a0) * i / n)
        out.append((cx + rx * math.cos(a), cy + ry * math.sin(a)))
    return out


def P(*segs):
    """Склеить сегменты в одну ломаную."""
    out = []
    for s in segs:
        if out and math.dist(out[-1], s[0]) < 1e-6:
            out.extend(s[1:])
        else:
            out.extend(s)
    return out


def G(strokes, leaves=()):
    return {"strokes": [list(s) for s in strokes], "leaves": list(leaves)}


def tf(g, sx=1.0, sy=1.0, dx=0.0, dy=0.0):
    """Масштаб + сдвиг скелета (толщина штриха при этом не меняется)."""
    return G([[(x * sx + dx, y * sy + dy) for x, y in s] for s in g["strokes"]],
             [(x * sx + dx, y * sy + dy, a, k) for x, y, a, k in g["leaves"]])


def add(g, *others):
    out = G(g["strokes"], g["leaves"])
    for o in others:
        out["strokes"] += o["strokes"]
        out["leaves"] += o["leaves"]
    return out


def leaf(x, y, angle=35, k=1.0):
    return G([], [(x, y, angle, k)])


def dot(x, y):
    return G([[(x, y), (x + 0.5, y + 0.5)]])


# ---------------------------------------------------------------- латиница: строчные

g = {}

g["o"] = G([arc(200, 250, 200, 250, 90, 450)])
g["a"] = G([arc(185, 250, 185, 250, 20, 375),
            P(ln((375, 500), (375, 70)), bz((375, 70), (375, 0), (410, -5), (450, 20)))])
g["b"] = G([ln((0, ASC), (0, 0)), arc(200, 250, 200, 250, -165, 165)])
g["c"] = G([arc(195, 250, 195, 250, 40, 318)])
g["d"] = G([ln((400, ASC), (400, 0)), arc(200, 250, 200, 250, 15, 345)])
g["e"] = G([P(ln((10, 255), (385, 255)), arc(200, 250, 190, 250, 0, 318))])
g["f"] = G([P(ln((100, 0), (100, 560)), bz((100, 560), (100, 740), (250, 770), (300, 690))),
            ln((0, 480), (250, 480))])
g["g"] = G([arc(185, 250, 185, 250, 0, 360),
            P(ln((375, 500), (375, 0)), bz((375, 0), (375, -240), (50, -260), (20, -120)))])
g["h"] = G([ln((0, ASC), (0, 0)),
            P(bz((0, 300), (40, 500), (350, 560), (350, 300)), ln((350, 300), (350, 0)))])
g["i"] = G([ln((0, 500), (0, 0))], [(10, 690, 40, 1.0)])
g["j"] = G([P(ln((150, 500), (150, -60)), bz((150, -60), (150, -220), (20, -240), (-20, -170)))],
           [(160, 690, 40, 1.0)])
g["k"] = G([ln((0, ASC), (0, 0)), ln((320, 500), (0, 190)), ln((110, 290), (330, 0))])
g["l"] = G([bz((0, ASC), (0, 300), (-5, 0), (100, 15))])
g["m"] = G([ln((0, 500), (0, 0)),
            P(bz((0, 300), (30, 500), (270, 560), (270, 300)), ln((270, 300), (270, 0))),
            P(bz((270, 300), (300, 500), (540, 560), (540, 300)), ln((540, 300), (540, 0)))])
g["n"] = G([ln((0, 500), (0, 0)),
            P(bz((0, 300), (40, 500), (350, 560), (350, 300)), ln((350, 300), (350, 0)))])
g["p"] = G([ln((0, 500), (0, DSC)), arc(200, 250, 200, 250, -165, 165)])
g["q"] = G([P(ln((400, 500), (400, DSC + 30)), bz((400, DSC + 30), (400, DSC - 10), (440, DSC - 10), (470, DSC + 20))),
            arc(200, 250, 200, 250, 15, 345)])
g["r"] = G([ln((0, 500), (0, 0)), bz((0, 310), (40, 480), (170, 540), (270, 460))])
s_low = [bz((330, 420), (290, 530), (30, 530), (30, 390)),
         bz((30, 390), (30, 250), (330, 290), (330, 130)),
         bz((330, 130), (330, -30), (40, -30), (0, 85))]
g["s"] = G([P(*s_low)])
g["t"] = G([P(ln((100, 690), (100, 130)), bz((100, 130), (100, -10), (180, -10), (250, 30))),
            ln((0, 480), (250, 480))])
g["u"] = G([P(ln((0, 500), (0, 190)), bz((0, 190), (0, -40), (350, -40), (350, 190))),
            ln((350, 500), (350, 0))])
g["v"] = G([ln((0, 500), (190, 0), (380, 500))])
g["w"] = G([ln((0, 500), (150, 0), (300, 380), (450, 0), (600, 500))])
g["x"] = G([ln((0, 500), (360, 0)), ln((360, 500), (0, 0))])
g["y"] = G([ln((0, 500), (195, 40)),
            P(ln((390, 500), (130, -150)), bz((130, -150), (105, -215), (60, -225), (20, -200)))])
g["z"] = G([ln((0, 500), (340, 500), (0, 0), (350, 0))])

# ---------------------------------------------------------------- латиница: прописные

g["A"] = G([ln((0, 0), (240, CH), (480, 0)), ln((95, 250), (385, 250))])
g["B"] = G([ln((0, 0), (0, CH)),
            P(ln((0, CH), (160, CH)), bz((160, CH), (330, CH), (330, 390), (160, 390)), ln((160, 390), (0, 390))),
            P(ln((0, 390), (190, 390)), bz((190, 390), (390, 390), (390, 0), (190, 0)), ln((190, 0), (0, 0)))])
g["C"] = G([arc(270, 350, 270, 350, 42, 320)])
g["D"] = G([ln((0, 0), (0, CH)),
            P(ln((0, CH), (150, CH)), bz((150, CH), (470, CH), (470, 0), (150, 0)), ln((150, 0), (0, 0)))])
g["E"] = G([ln((360, CH), (0, CH), (0, 0), (370, 0)), ln((0, 370), (300, 370))])
g["F"] = G([ln((360, CH), (0, CH), (0, 0)), ln((0, 370), (300, 370))])
g["G"] = G([P(arc(275, 350, 275, 350, 45, 360), ln((550, 350), (310, 350)))])
g["H"] = G([ln((0, 0), (0, CH)), ln((440, 0), (440, CH)), ln((0, 370), (440, 370))])
g["I"] = G([ln((0, 0), (0, CH))])
g["J"] = G([P(ln((300, CH), (300, 210)), bz((300, 210), (300, -40), (30, -40), (0, 160)))])
g["K"] = G([ln((0, 0), (0, CH)), ln((380, CH), (0, 250)), ln((130, 390), (410, 0))])
g["L"] = G([ln((0, CH), (0, 0), (350, 0))])
g["M"] = G([ln((0, 0), (30, CH), (270, 170), (510, CH), (540, 0))])
g["N"] = G([ln((0, 0), (0, CH), (440, 0), (440, CH))])
g["O"] = G([arc(295, 350, 295, 350, 90, 450)])
g["P"] = G([ln((0, 0), (0, CH)),
            P(ln((0, CH), (170, CH)), bz((170, CH), (380, CH), (380, 320), (170, 320)), ln((170, 320), (0, 320)))])
g["Q"] = G([arc(295, 350, 295, 350, 90, 450), ln((340, 170), (610, -70))])
g["R"] = add(g["P"], G([ln((170, 320), (410, 0))]))
g["S"] = tf(g["s"], 1.18, 1.4, 0, 0)
g["T"] = G([ln((0, CH), (460, CH)), ln((230, CH), (230, 0))])
g["U"] = G([P(ln((0, CH), (0, 250)), bz((0, 250), (0, -40), (450, -40), (450, 250)), ln((450, 250), (450, CH)))])
g["V"] = G([ln((0, CH), (240, 0), (480, CH))])
g["W"] = G([ln((0, CH), (170, 0), (350, 540), (530, 0), (700, CH))])
g["X"] = G([ln((0, CH), (450, 0)), ln((450, CH), (0, 0))])
g["Y"] = G([ln((0, CH), (230, 350), (460, CH)), ln((230, 350), (230, 0))])
g["Z"] = G([ln((0, CH), (420, CH), (0, 0), (440, 0))])

# ---------------------------------------------------------------- цифры

g["0"] = G([arc(230, 350, 230, 350, 90, 450)])
g["1"] = G([ln((20, 540), (210, CH), (210, 0))])
g["2"] = G([P(bz((20, 540), (50, 740), (390, 750), (380, 520)),
              bz((380, 520), (370, 340), (90, 230), (0, 0)), ln((0, 0), (400, 0)))])
three = [bz((20, 600), (80, 745), (390, 735), (370, 545)),
         bz((370, 545), (355, 410), (230, 390), (160, 390)),
         bz((160, 390), (430, 390), (450, 180), (380, 85)),
         bz((380, 85), (300, -35), (60, -25), (0, 100))]
g["3"] = G([P(*three)])
g["4"] = G([ln((310, 0), (310, CH), (0, 220), (430, 220))])
g["5"] = G([P(ln((390, CH), (60, CH), (35, 400)),
              bz((35, 400), (150, 480), (420, 470), (410, 235)),
              bz((410, 235), (400, -35), (80, -40), (0, 95)))])
g["6"] = G([bz((350, 690), (120, 650), (0, 430), (0, 220)),
            arc(205, 215, 205, 215, 0, 360)])
g["7"] = G([ln((0, CH), (410, CH), (130, 0)), ln((110, 350), (340, 350))])
g["8"] = G([arc(205, 535, 165, 165, -90, 270), arc(205, 190, 205, 190, 90, 450)])
g["9"] = G([arc(205, 485, 205, 215, 0, 360),
            bz((410, 480), (410, 260), (290, 40), (60, 15))])

# ---------------------------------------------------------------- кириллица: прописные

cy = {}
cy["А"] = g["A"]
cy["Б"] = G([P(ln((370, CH), (0, CH), (0, 0), (180, 0)),
               bz((180, 0), (400, 0), (400, 410), (180, 410)), ln((180, 410), (0, 410)))])
cy["В"] = g["B"]
cy["Г"] = G([ln((0, 0), (0, CH), (360, CH))])
cy["Д"] = G([P(ln((0, -140), (0, 0), (500, 0), (500, -140))),
             P(bz((60, 0), (140, 200), (150, 500), (150, CH)), ln((150, CH), (440, CH), (440, 0)))])
cy["Е"] = g["E"]
cy["Ё"] = add(g["E"], leaf(100, 840, 40), leaf(270, 840, 40))
cy["Ж"] = G([ln((300, 0), (300, CH)),
             ln((10, CH), (300, 370), (0, 0)), ln((590, CH), (300, 370), (600, 0))])
cy["З"] = G([P(*three)])
cy["И"] = G([ln((0, CH), (0, 0), (440, CH), (440, 0))])
cy["Й"] = add(cy["И"], G([arc(220, 870, 110, 75, 200, 340)]))
cy["К"] = g["K"]
cy["Л"] = G([P(bz((0, 0), (110, 10), (150, 300), (150, CH)), ln((150, CH), (440, CH), (440, 0)))])
cy["М"] = g["M"]
cy["Н"] = g["H"]
cy["О"] = g["O"]
cy["П"] = G([ln((0, 0), (0, CH), (440, CH), (440, 0))])
cy["Р"] = g["P"]
cy["С"] = g["C"]
cy["Т"] = g["T"]
cy["У"] = G([ln((0, CH), (230, 270)),
             P(ln((470, CH), (190, 60)), bz((190, 60), (160, 0), (90, -15), (30, 20)))])
cy["Ф"] = G([ln((290, CH + 30), (290, -30)), arc(290, 380, 290, 215, 90, 450)])
cy["Х"] = g["X"]
cy["Ц"] = G([ln((0, CH), (0, 0), (450, 0)), ln((450, CH), (450, 0), (520, 0), (520, -140))])
cy["Ч"] = G([P(ln((0, CH), (0, 450)), bz((0, 450), (0, 290), (170, 270), (410, 320))),
             ln((410, CH), (410, 0))])
cy["Ш"] = G([ln((0, CH), (0, 0), (580, 0), (580, CH)), ln((290, CH), (290, 0))])
cy["Щ"] = G([ln((0, CH), (0, 0), (580, 0)), ln((580, CH), (580, 0), (650, 0), (650, -140)),
             ln((290, CH), (290, 0))])
soft = [ln((0, 410), (190, 410)), bz((190, 410), (410, 410), (410, 0), (190, 0)), ln((190, 0), (0, 0))]
cy["Ь"] = G([ln((0, CH), (0, 0)), P(*soft)])
cy["Ъ"] = add(tf(cy["Ь"], 1, 1, 140, 0), G([ln((0, CH), (140, CH))]))
cy["Ы"] = add(cy["Ь"], G([ln((540, CH), (540, 0))]))
cy["Э"] = G([arc(260, 350, 270, 350, 138, -140), ln((130, 360), (520, 360))])
cy["Ю"] = G([ln((0, 0), (0, CH)), ln((0, 360), (160, 360)), arc(400, 350, 230, 350, 90, 450)])
cy["Я"] = G([ln((370, 0), (370, CH)),
             P(ln((370, CH), (200, CH)), bz((200, CH), (-10, CH), (-10, 320), (200, 320)), ln((200, 320), (370, 320))),
             ln((210, 320), (-20, 0))])

# ---------------------------------------------------------------- кириллица: строчные


def tail(x, y=60):
    """Правая ножка с маленьким хвостиком вправо (как у a и q)."""
    return bz((x, y), (x, 0), (x + 30, -5), (x + 70, 20))


soft_l = [ln((0, 290), (150, 290)), bz((150, 290), (320, 290), (320, 0), (150, 0)), ln((150, 0), (0, 0))]

cyl = {}
cyl["а"] = g["a"]
cyl["б"] = G([arc(200, 225, 200, 225, 0, 360),
             bz((25, 330), (40, 620), (190, 690), (390, 735))])
cyl["в"] = G([ln((0, 0), (0, 500)),
             P(ln((0, 500), (140, 500)), bz((140, 500), (285, 500), (285, 275), (140, 275)), ln((140, 275), (0, 275))),
             P(ln((0, 275), (160, 275)), bz((160, 275), (325, 275), (325, 0), (160, 0)), ln((160, 0), (0, 0)))])
cyl["г"] = G([P(ln((0, 0), (0, 500), (220, 500)), bz((220, 500), (250, 500), (275, 490), (290, 465)))])
cyl["д"] = G([P(ln((0, -140), (0, 0), (420, 0), (420, -140))),
             P(bz((45, 0), (110, 140), (125, 330), (125, 500)), ln((125, 500), (355, 500), (355, 0)))])
cyl["е"] = g["e"]
cyl["ё"] = add(g["e"], leaf(105, 675, 40), leaf(285, 675, 40))
cyl["ж"] = G([ln((265, 0), (265, 500)),
             P(bz((20, 500), (90, 440), (180, 260), (265, 255)), bz((265, 255), (170, 250), (60, 90), (0, 0))),
             P(bz((510, 500), (440, 440), (350, 260), (265, 255)), bz((265, 255), (360, 250), (470, 90), (530, 0)))])
cyl["з"] = G([P(bz((10, 410), (50, 530), (300, 525), (295, 385)),
                bz((295, 385), (290, 295), (190, 275), (130, 275)),
                bz((130, 275), (330, 275), (345, 125), (290, 55)),
                bz((290, 55), (220, -25), (45, -20), (0, 75)))])
cyl["и"] = G([P(ln((0, 500), (0, 0), (350, 500), (350, 60)), tail(350))])
cyl["й"] = add(cyl["и"], G([arc(175, 660, 105, 75, 200, 340)]))
cyl["к"] = G([ln((0, 0), (0, 500)), bz((300, 500), (220, 440), (120, 270), (0, 255)),
             P(bz((40, 262), (160, 255), (220, 100), (330, 0)))])
cyl["л"] = G([P(bz((0, 0), (100, 10), (125, 250), (125, 500)), ln((125, 500), (365, 500), (365, 60)), tail(365))])
cyl["м"] = G([ln((0, 0), (35, 500), (235, 130), (435, 500), (470, 0))])
cyl["н"] = G([ln((0, 0), (0, 500)), P(ln((345, 500), (345, 60)), tail(345)), ln((0, 260), (345, 260))])
cyl["о"] = g["o"]
cyl["п"] = G([ln((0, 0), (0, 500)),
             P(bz((0, 440), (40, 500), (90, 500), (150, 500)), ln((150, 500), (250, 500)),
               bz((250, 500), (320, 500), (345, 470), (345, 400)), ln((345, 400), (345, 0)))])
cyl["р"] = g["p"]
cyl["с"] = g["c"]
cyl["т"] = G([ln((0, 500), (420, 500)), ln((210, 500), (210, 0))])
cyl["у"] = g["y"]
cyl["ф"] = G([ln((230, ASC), (230, DSC)), arc(230, 250, 230, 230, 90, 450)])
cyl["х"] = g["x"]
cyl["ц"] = G([ln((0, 500), (0, 0), (355, 0)), ln((355, 500), (355, 0), (415, 0), (415, -135))])
cyl["ч"] = G([P(ln((0, 500), (0, 330)), bz((0, 330), (0, 200), (150, 185), (325, 235))),
             P(ln((325, 500), (325, 60)), tail(325))])
cyl["ш"] = G([ln((0, 500), (0, 0), (510, 0), (510, 500)), ln((255, 500), (255, 0))])
cyl["щ"] = G([ln((0, 500), (0, 0), (510, 0)), ln((510, 500), (510, 0), (570, 0), (570, -135)),
             ln((255, 500), (255, 0))])
cyl["ь"] = G([ln((0, 500), (0, 0)), P(*soft_l)])
cyl["ъ"] = add(tf(cyl["ь"], 1, 1, 115, 0), G([ln((0, 500), (115, 500))]))
cyl["ы"] = add(cyl["ь"], G([ln((440, 500), (440, 0))]))
cyl["э"] = G([arc(200, 250, 200, 250, 140, -140), ln((95, 255), (395, 255))])
cyl["ю"] = G([ln((0, 0), (0, 500)), ln((0, 255), (125, 255)), arc(320, 250, 195, 250, 90, 450)])
cyl["я"] = G([ln((335, 0), (335, 500)),
             P(ln((335, 500), (175, 500)), bz((175, 500), (15, 500), (15, 225), (175, 225)), ln((175, 225), (335, 225))),
             bz((185, 225), (110, 200), (40, 90), (0, 0))])

# ---------------------------------------------------------------- знаки

pun = {}
pun["."] = dot(0, 0)
pun[","] = G([bz((10, 20), (15, -40), (0, -90), (-30, -130))])
pun[":"] = add(dot(0, 0), dot(0, 380))
pun[";"] = add(pun[","], dot(10, 380))
pun["!"] = G([ln((0, CH), (0, 220))], [(0, 10, 90, 0.9)])
pun["?"] = G([P(bz((0, 560), (20, 740), (330, 750), (320, 540)),
                bz((320, 540), (310, 390), (150, 400), (150, 220)))], [(150, 10, 90, 0.9)])
pun["-"] = G([ln((0, 260), (230, 260))])
pun["–"] = G([ln((0, 260), (420, 260))])
pun["—"] = G([ln((0, 260), (760, 260))])
pun["'"] = G([ln((0, CH + 20), (0, 520))])
pun['"'] = G([ln((0, CH + 20), (0, 520)), ln((130, CH + 20), (130, 520))])
pun["’"] = G([bz((20, CH + 20), (25, 640), (10, 580), (-20, 540))])
pun["‘"] = G([bz((0, 540), (-5, 620), (10, 680), (40, CH + 20))])
pun["«"] = G([ln((120, 420), (0, 250), (120, 80)), ln((280, 420), (160, 250), (280, 80))])
pun["»"] = G([ln((0, 420), (120, 250), (0, 80)), ln((160, 420), (280, 250), (160, 80))])
pun["("] = G([bz((170, 780), (-40, 600), (-40, 0), (170, -180))])
pun[")"] = G([bz((0, 780), (210, 600), (210, 0), (0, -180))])
pun["["] = G([ln((150, 780), (0, 780), (0, -180), (150, -180))])
pun["]"] = G([ln((0, 780), (150, 780), (150, -180), (0, -180))])
pun["/"] = G([ln((0, -120), (330, 760))])
pun["\\"] = G([ln((0, 760), (330, -120))])
pun["+"] = G([ln((0, 280), (380, 280)), ln((190, 90), (190, 470))])
pun["="] = G([ln((0, 360), (380, 360)), ln((0, 190), (380, 190))])
pun["*"] = G([ln((150, CH), (150, 440)), ln((30, 640), (270, 500)), ln((30, 500), (270, 640))])
pun["%"] = G([arc(110, 560, 110, 140, 90, 450), arc(470, 140, 110, 140, 90, 450), ln((0, 0), (580, CH))])
pun["_"] = G([ln((0, -120), (480, -120))])
pun["№"] = add(G([ln((0, 0), (0, CH), (400, 0), (400, CH))]),
                    G([arc(590, 480, 110, 140, 90, 450), ln((470, 220), (710, 220))]))
pun["@"] = G([P(arc(330, 300, 130, 160, 0, 360)),
              P(ln((460, 380), (460, 200)), bz((460, 200), (460, 40), (680, 60), (680, 300)),
                arc(340, 300, 340, 380, 0, 330))])
pun["&"] = G([P(bz((520, 0), (300, 250), (60, 460), (130, 620)),
                bz((130, 620), (200, 760), (400, 720), (350, 560)),
                bz((350, 560), (300, 420), (0, 380), (0, 160)),
                bz((0, 160), (0, -40), (330, -60), (480, 260)))])
pun["#"] = G([ln((140, 0), (200, CH)), ln((340, 0), (400, CH)), ln((40, 470), (500, 470)), ln((20, 230), (480, 230))])
pun["…"] = add(dot(0, 0), dot(170, 0), dot(340, 0))
pun["₽"] = G([ln((60, 0), (60, CH)),
                   P(ln((60, CH), (230, CH)), bz((230, CH), (440, CH), (440, 320), (230, 320)), ln((230, 320), (0, 320))),
                   ln((0, 170), (300, 170))])

# ---------------------------------------------------------------- отрисовка


def phase(name, k):
    h = hashlib.md5(f"{name}:{k}".encode()).digest()
    return h[0] / 255 * 2 * math.pi


def resample(pts, step=6.0):
    out = [pts[0]]
    for a, b in zip(pts, pts[1:]):
        d = math.dist(a, b)
        n = max(1, int(d / step))
        for i in range(1, n + 1):
            t = i / n
            out.append((a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t))
    return out


def wobble(pts, ph1, ph2, amp=5.0):
    return [(x + amp * math.sin(y / 95 + ph1), y + amp * 0.7 * math.sin(x / 110 + ph2)) for x, y in pts]


def stroke_poly(pts, ph):
    pts = resample(pts)
    if len(pts) == 1:
        pts = [pts[0], (pts[0][0] + 0.5, pts[0][1] + 0.5)]
    closed = math.dist(pts[0], pts[-1]) < 2
    cum = [0.0]
    for a, b in zip(pts, pts[1:]):
        cum.append(cum[-1] + math.dist(a, b))
    L = cum[-1]
    radii = []
    for s in cum:
        r = R * (1 + 0.045 * math.sin(s / 85 + ph))
        if not closed:
            d = min(s, L - s)
            r *= 1 + 0.13 * math.exp(-(d / 22) ** 2)   # «капля» чернил на конце
        radii.append(r)
    circles = [Point(p).buffer(r, quad_segs=10) for p, r in zip(pts, radii)]
    if len(circles) == 1:
        return circles[0]
    hulls = [circles[i].union(circles[i + 1]).convex_hull for i in range(len(circles) - 1)]
    return unary_union(hulls)


def leaf_poly(x, y, angle, k):
    Lh, W = 82 * k, 40 * k
    top, bot = [], []
    n = 24
    for i in range(n + 1):
        t = i / n
        px = -Lh + 2 * Lh * t
        w = W * math.sin(math.pi * t) ** 0.85 * (1 + 0.25 * (t - 0.5))
        top.append((px, w))
        bot.append((px, -w * 0.8))
    poly = Polygon(top + bot[::-1][1:-1])
    poly = affinity.rotate(poly, angle, origin=(0, 0))
    return affinity.translate(poly, x, y)


def render(name, glyph, sx=1.0, rot=0.0, dy=0.0, seed=""):
    strokes = glyph["strokes"]
    ph1, ph2 = phase(name + seed, 1), phase(name + seed, 2)
    parts = []
    for i, s in enumerate(strokes):
        s = [(x * sx, y + dy) for x, y in s]
        s = wobble(resample(s, 10), ph1, ph2)
        parts.append(stroke_poly(s, phase(name + seed, 10 + i)))
    for x, y, a, k in glyph["leaves"]:
        parts.append(leaf_poly(x * sx, y + dy, a, k))
    shape = unary_union(parts) if parts else None
    if shape is not None and rot:
        c = shape.centroid
        shape = affinity.rotate(shape, rot, origin=(c.x, 250))
    if shape is not None:
        shape = shape.simplify(0.8, preserve_topology=True)
    return shape


def polys_of(shape):
    if shape is None:
        return []
    if isinstance(shape, Polygon):
        return [shape]
    if isinstance(shape, MultiPolygon):
        return list(shape.geoms)
    return [p for p in getattr(shape, "geoms", []) if isinstance(p, Polygon)]


def draw(shape, pen, shift):
    for poly in polys_of(shape):
        poly = orient(poly, sign=-1.0)  # внешний контур по часовой (TrueType)
        rings = [poly.exterior] + list(poly.interiors)
        for ring in rings:
            pts = [(round(x + shift), round(y)) for x, y in list(ring.coords)[:-1]]
            clean = []
            for p in pts:
                if not clean or p != clean[-1]:
                    clean.append(p)
            if len(clean) > 1 and clean[0] == clean[-1]:
                clean.pop()
            if len(clean) < 3:
                continue
            pen.moveTo(clean[0])
            for p in clean[1:]:
                pen.lineTo(p)
            pen.closePath()


# ---------------------------------------------------------------- сборка

def gname(ch):
    cp = ord(ch)
    if ch.isascii() and ch.isalnum():
        return ch if ch.isalpha() else {"0": "zero", "1": "one", "2": "two", "3": "three", "4": "four",
                                        "5": "five", "6": "six", "7": "seven", "8": "eight", "9": "nine"}[ch]
    return f"uni{cp:04X}"


KERN_TARGET = 2 * SB      # «нормальный» просвет, как между двумя вертикальными штрихами
KERN_TIGHTEN = 0.55       # какую долю лишнего просвета убирать
KERN_MIN, KERN_MAX = -130, 50
KERN_STEP = 20            # шаг профиля по высоте
KERN_BLUR = 2             # учитывать соседние строки (±40 ед.) для диагоналей и округлостей


def kerning(shapes, letters):
    """Автокернинг по профилям: для каждой пары меряем минимальный
    просвет между реальными контурами и подтягиваем его к KERN_TARGET.
    Альтернативы букв кернятся тем же классом, что и основная форма."""
    from shapely.geometry import box
    ys = list(range(-260, 960, KERN_STEP))
    bases = [n for n in shapes if "." not in n]
    prof = {}
    for n in bases:
        shape, shift, adv = shapes[n]
        minx, _, maxx, _ = shape.bounds
        left, right = [], []
        for y in ys:
            part = shape.intersection(box(minx - 1, y - KERN_STEP / 2, maxx + 1, y + KERN_STEP / 2))
            if part.is_empty:
                left.append(None)
                right.append(None)
            else:
                b = part.bounds
                left.append(b[0] + shift)
                right.append(adv - (b[2] + shift))
        # «размытие»: у соседних строк берём ближайшую к краю точку
        def blur(arr):
            out = []
            for i in range(len(arr)):
                win = [v for v in arr[max(0, i - KERN_BLUR): i + KERN_BLUR + 1] if v is not None]
                out.append(min(win) if win else None)
            return out
        prof[n] = (blur(left), blur(right))

    pairs = {}
    for a in bases:
        ra = prof[a][1]
        for b in bases:
            lb = prof[b][0]
            gaps = [x + y for x, y in zip(ra, lb) if x is not None and y is not None]
            if len(gaps) < 2:
                continue
            gap = min(gaps)
            if gap > KERN_TARGET:
                k = -(gap - KERN_TARGET) * KERN_TIGHTEN
            else:
                k = (KERN_TARGET - gap) * 0.8
            k = max(KERN_MIN, min(KERN_MAX, k))
            k = int(round(k / 5.0) * 5)
            if abs(k) >= 15:
                pairs[(a, b)] = k

    def cls(n):
        return f"@k_{n}"

    lines = []
    for n in bases:
        members = [n] + ([f"{n}.alt1", f"{n}.alt2"] if n in letters else [])
        lines.append(f"{cls(n)} = [{' '.join(members)}];")
    lines.append("feature kern {")
    for (a, b), k in sorted(pairs.items()):
        lines.append(f"    pos {cls(a)} {cls(b)} {k};")
    lines.append("} kern;")
    print("kerning pairs:", len(pairs))
    return "\n".join(lines)


ALTS = {  # вариант: (ширина, поворот°, сдвиг по вертикали)
    "alt1": (0.92, 2.2, 14),
    "alt2": (1.07, -1.6, -10),
}


def main():
    sources = {}
    sources.update(g)
    sources.update(cy)
    sources.update(cyl)
    sources.update(pun)

    shapes = {}       # glyph name -> (shape, advance)
    cmap = {}
    letters = []

    def place(name, shape, extra=0):
        if shape is None:
            return
        minx, _, maxx, _ = shape.bounds
        shifted = -minx + SB
        adv = int(round(maxx - minx + 2 * SB + extra))
        shapes[name] = (shape, shifted, adv)

    for ch, gl in sources.items():
        name = gname(ch)
        cmap[ord(ch)] = name
        place(name, render(name, gl))
        if ch.isalpha():
            letters.append(name)
            for alt, (sx, rot, dy) in ALTS.items():
                place(f"{name}.{alt}", render(name, gl, sx, rot, dy, seed=alt))

    # пробелы
    cmap[0x20] = "space"
    cmap[0xA0] = "uni00A0"

    order = [".notdef", "space", "uni00A0"] + list(shapes.keys())
    kern_fea = kerning(shapes, set(letters))

    def build(fmt):
        fb = FontBuilder(UPM, isTTF=(fmt == "ttf"))
        fb.setupGlyphOrder(order)
        fb.setupCharacterMap(cmap)
        metrics = {".notdef": (500, 50), "space": (260, 0), "uni00A0": (260, 0)}
        glyphs = {}

        def empty():
            return TTGlyphPen(None).glyph() if fmt == "ttf" else T2CharStringPen(260, None).getCharString()

        def notdef():
            pen = TTGlyphPen(None) if fmt == "ttf" else T2CharStringPen(500, None)
            for rect in ([(50, 0), (50, 700), (450, 700), (450, 0)], [(100, 50), (400, 50), (400, 650), (100, 650)]):
                pen.moveTo(rect[0])
                for p in rect[1:]:
                    pen.lineTo(p)
                pen.closePath()
            return pen.glyph() if fmt == "ttf" else pen.getCharString()

        glyphs[".notdef"] = notdef()
        glyphs["space"] = empty()
        glyphs["uni00A0"] = empty()
        for name, (shape, shift, adv) in shapes.items():
            pen = TTGlyphPen(None) if fmt == "ttf" else T2CharStringPen(adv, None)
            draw(shape, pen, shift)
            glyphs[name] = pen.glyph() if fmt == "ttf" else pen.getCharString()
            minx = shape.bounds[0] + shift
            metrics[name] = (adv, int(round(minx)))

        if fmt == "ttf":
            fb.setupGlyf(glyphs)
        else:
            fb.setupCFF(FAMILY + "-Regular", {"FullName": FAMILY + " Regular"}, glyphs, {})
        from fontTools.misc.timeTools import timestampNow
        y_max = max(sh.bounds[3] for sh, _, _ in shapes.values())
        y_min = min(sh.bounds[1] for sh, _, _ in shapes.values())
        asc = int(math.ceil(max(y_max, 1000)))
        desc = int(math.ceil(max(-y_min, 300)))
        fb.setupHead(fontRevision=1.1, created=timestampNow(), modified=timestampNow())
        fb.setupHorizontalMetrics(metrics)
        fb.setupHorizontalHeader(ascent=asc, descent=-desc, lineGap=0)
        fb.setupNameTable({
            "familyName": FAMILY,
            "styleName": "Regular",
            "uniqueFontIdentifier": f"{FAMILY}-Regular-1.100",
            "fullName": f"{FAMILY} Regular",
            "psName": f"{FAMILY}-Regular",
            "version": "Version 1.100",
            "designer": "Kapelka project",
            "description": "Рукописный дисплейный шрифт: фломастерный штрих, листики вместо точек, прыгающая строка.",
            "copyright": "Copyright 2026 Gryadka Project Authors. All rights reserved.",
            "licenseDescription": "Use rights are granted by the font owner under a separate written license.",
        }, mac=False)
        fb.setupOS2(sTypoAscender=asc, sTypoDescender=-desc, sTypoLineGap=0,
                    usWinAscent=asc, usWinDescent=desc, sxHeight=XH, sCapHeight=CH,
                    achVendID="KPLK", fsType=0, version=4, fsSelection=0x40 | 0x80,
                    ulUnicodeRange1=(1 << 0) | (1 << 1) | (1 << 9),
                    ulCodePageRange1=(1 << 0) | (1 << 2))
        fb.setupPost()
        if fmt == "ttf":                  # сглаживание без хинтинга
            from fontTools.ttLib import newTable
            from fontTools.ttLib.tables import ttProgram
            prep = newTable("prep")
            prep.program = ttProgram.Program()
            prep.program.fromBytecode(bytes([0xB8, 0x01, 0xFF, 0x85, 0xB0, 0x04, 0x8D]))
            fb.font["prep"] = prep
            gasp = newTable("gasp")
            gasp.version = 1
            gasp.gaspRange = {0xFFFF: 0x000F}
            fb.font["gasp"] = gasp

        base = " ".join(letters)
        a1 = " ".join(n + ".alt1" for n in letters)
        a2 = " ".join(n + ".alt2" for n in letters)
        fea = f"""
languagesystem DFLT dflt;
languagesystem latn dflt;
languagesystem cyrl dflt;

@base = [{base}];
@alt1 = [{a1}];
@alt2 = [{a2}];

lookup bounce {{
    sub @alt1 @base' by @alt2;
    sub @base @base' by @alt1;
}} bounce;

feature calt {{ lookup bounce; }} calt;
feature salt {{ sub @base by @alt1; }} salt;
feature ss01 {{ featureNames {{ name "Wide alternates"; }}; sub @base by @alt2; }} ss01;

{kern_fea}
"""
        addOpenTypeFeaturesFromString(fb.font, fea)
        out = os.path.join(HERE, f"{FAMILY}-Regular.{fmt}")
        fb.save(out)
        print("saved", out, len(order), "glyphs")
        return out

    ttf = build("ttf")
    build("otf")
    specimen(ttf)


def specimen(path):
    from PIL import Image, ImageDraw, ImageFont
    kern_test(path)
    W, H = 1800, 1250
    img = Image.new("RGB", (W, H), (250, 248, 240))
    d = ImageDraw.Draw(img)
    green = (40, 120, 60)
    ink = (30, 35, 30)
    try:
        f_big = ImageFont.truetype(path, 150, layout_engine=ImageFont.Layout.RAQM)
        f_mid = ImageFont.truetype(path, 72, layout_engine=ImageFont.Layout.RAQM)
        f_sm = ImageFont.truetype(path, 56, layout_engine=ImageFont.Layout.RAQM)
    except Exception:
        f_big = ImageFont.truetype(path, 150)
        f_mid = ImageFont.truetype(path, 72)
        f_sm = ImageFont.truetype(path, 56)
    y = 30
    d.text((60, y), "Капелька", font=f_big, fill=green); y += 210
    for line in ["АБВГДЕЁЖЗИЙКЛМНОПРСТУФХЦЧШЩЪЫЬЭЮЯ",
                 "абвгдеёжзийклмнопрстуфхцчшщъыьэюя",
                 "ABCDEFGHIJKLMNOPQRSTUVWXYZ",
                 "abcdefghijklmnopqrstuvwxyz",
                 "0123456789 .,:;!?-–— «»()[]/+=*%№@&#…₽"]:
        d.text((60, y), line, font=f_sm, fill=ink); y += 105
    y += 20
    d.text((60, y), "Свежие овощи прямо с грядки!", font=f_mid, fill=green); y += 115
    d.text((60, y), "Вкусно, полезно и с любовью!", font=f_mid, fill=ink); y += 115
    d.text((60, y), "Fresh juicy apples — 99 ₽/kg", font=f_mid, fill=ink)
    out = os.path.join(HERE, "specimen.png")
    img.save(out)
    print("saved", out)


def kern_test(path):
    """Сравнение: без кернинга / с кернингом."""
    from PIL import Image, ImageDraw, ImageFont
    f = ImageFont.truetype(path, 90, layout_engine=ImageFont.Layout.RAQM)
    lab = ImageFont.truetype(path, 40, layout_engine=ImageFont.Layout.RAQM)
    texts = ["ТОВАР «Уют», Год. LT AV", "Гусь, Учёт, Ягода. Yes, To"]
    img = Image.new("RGB", (1800, 680), (250, 248, 240))
    d = ImageDraw.Draw(img)
    y = 20
    for t in texts:
        d.text((60, y), "без кернинга", font=lab, fill=(150, 60, 50)); y += 45
        d.text((60, y), t, font=f, fill=(30, 35, 30), features=["-kern"]); y += 110
        d.text((60, y), "с кернингом", font=lab, fill=(40, 120, 60)); y += 45
        d.text((60, y), t, font=f, fill=(30, 35, 30)); y += 110
    img.save(os.path.join(HERE, "kerning.png"))


if __name__ == "__main__":
    main()
