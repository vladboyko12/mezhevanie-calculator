#!/usr/bin/env python3
"""
Сравнение «Капельки» с образцом другого шрифта (например, скриншотом Villula).

    python3 compare.py reference.png "ВкусВилл — вкусно и полезно"

Сверху — картинка-образец, снизу — тот же текст, набранный Kapelka.
Результат: compare.png
"""
import os
import sys

from PIL import Image, ImageDraw, ImageFont

HERE = os.path.dirname(os.path.abspath(__file__))
FONT = os.path.join(HERE, "Kapelka-Regular.otf")
W = 1800
BG = (250, 248, 240)


def main(ref_path, text):
    ref = Image.open(ref_path).convert("RGB")
    ref = ref.resize((W - 120, int(ref.height * (W - 120) / ref.width)))
    label = ImageFont.truetype(FONT, 44, layout_engine=ImageFont.Layout.RAQM)

    # подбираем кегль, чтобы строка заняла ширину образца
    size = 200
    while size > 20:
        f = ImageFont.truetype(FONT, size, layout_engine=ImageFont.Layout.RAQM)
        if f.getbbox(text)[2] <= W - 120:
            break
        size -= 4
    th = f.getbbox(text)[3] + 40

    img = Image.new("RGB", (W, ref.height + th + 220), BG)
    d = ImageDraw.Draw(img)
    d.text((60, 20), "Образец", font=label, fill=(150, 60, 50))
    img.paste(ref, (60, 85))
    y = 85 + ref.height + 30
    d.text((60, y), "Kapelka", font=label, fill=(40, 120, 60))
    d.text((60, y + 65), text, font=f, fill=(30, 35, 30))
    out = os.path.join(HERE, "compare.png")
    img.save(out)
    print("saved", out)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
