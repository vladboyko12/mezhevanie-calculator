#!/usr/bin/env python3
"""Лист покрытия: все знаки расширенного набора Gryadka. Результат: gryadka-coverage.png"""
import os
from PIL import Image, ImageDraw, ImageFont

HERE = os.path.dirname(os.path.abspath(__file__))
ROWS = [
    "ÀÁÂÃÄÅÆÇÈÉÊËÌÍÎÏÐÑÒÓÔÕÖØÙÚÛÜÝÞẞ",
    "àáâãäåæçèéêëìíîïðñòóôõöøùúûüýþÿß",
    "ĀĂĄĆĈĊČĎĐĒĔĖĘĚĜĞĠĢĤĦĨĪĬĮİĲĴĶĹĻĽĿŁ",
    "āăąćĉċčďđēĕėęěĝğġģĥħĩīĭįıĳĵķĸĺļľŀł",
    "ŃŅŇŊŌŎŐŒŔŖŘŚŜŞŠŢŤŦŨŪŬŮŰŲŴŶŸŹŻŽȘȚ",
    "ńņňŋōŏőœŕŗřśŝşšţťŧũūŭůűųŵŷźżžșț",
    "ЀЂЃЄЅІЇЈЉЊЋЌЍЎЏҐҒҚҢҮҰҺӘӨ",
    "ѐђѓєѕіїјљњћќѝўџґғқңүұһәө",
    "?¿!¡&@„“”‚‘’«»‹›(){}[]<>≤≥±×÷−~≈≠^`|¦′″·•",
    "$€£¥¢₽°©®™№%#*+=/\\_…",
]


def main():
    styles = ["Light", "Regular", "Black"]
    size, line = 64, 96
    img = Image.new("RGB", (2200, 60 + len(ROWS) * line * len(styles) + 40 * len(ROWS)), (250, 248, 240))
    d = ImageDraw.Draw(img)
    fonts = {s: ImageFont.truetype(os.path.join(HERE, f"Gryadka-{s}.otf"), size,
                                   layout_engine=ImageFont.Layout.RAQM) for s in styles}
    lab = ImageFont.truetype(os.path.join(HERE, "Gryadka-Regular.otf"), 26, layout_engine=ImageFont.Layout.RAQM)
    y = 30
    for row in ROWS:
        for s in styles:
            d.text((20, y + 30), s, font=lab, fill=(130, 130, 120))
            d.text((180, y), row, font=fonts[s], fill=(30, 35, 30))
            y += line
        d.line((20, y + 5, 2180, y + 5), fill=(215, 210, 195), width=2)
        y += 40
    img.save(os.path.join(HERE, "gryadka-coverage.png"))
    print("saved coverage")


if __name__ == "__main__":
    main()
