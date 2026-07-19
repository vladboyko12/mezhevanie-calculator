// Тарифы и логика расчёта стоимости межевания.
// Классический скрипт (без ES-модулей): CONFIG и calculate становятся
// глобальными при подключении через <script src="calc.js">. Так страница
// работает и при открытии файла напрямую (file://), и на хостинге.

const CONFIG = {
  base: 8000,            // базовая стоимость: оформление, документы
  perAre: 500,           // цена за одну сотку
  shapeMultiplier: {
    rectangular: 1.0,    // прямоугольный участок
    irregular:   1.2,    // неровный участок (сложнее замер)
  },
  weekendSurcharge: 1000, // надбавка за выезд в выходной
  zones: [
    { id: 'z1', label: 'До 30 км',  price: 1500 },
    { id: 'z2', label: '30–60 км',  price: 3000 },
    { id: 'z3', label: '60–100 км', price: 5000 },
    { id: 'z4', label: '100+ км',   price: 7000 },
  ],
};

function calculate(input) {
  const { ares, shape, zoneId, isWeekend } = input;

  const zone = CONFIG.zones.find(z => z.id === zoneId);
  const multiplier = CONFIG.shapeMultiplier[shape];

  const empty = { survey: 0, travel: 0, weekendFee: 0, total: 0 };

  if (typeof ares !== 'number' || !isFinite(ares) || ares <= 0) {
    return { ...empty, valid: false, error: 'Укажите площадь больше нуля' };
  }
  if (!zone) {
    return { ...empty, valid: false, error: 'Выберите район' };
  }
  if (multiplier === undefined) {
    return { ...empty, valid: false, error: 'Выберите форму участка' };
  }

  const survey = Math.round((CONFIG.base + CONFIG.perAre * ares) * multiplier);
  const weekendFee = isWeekend ? CONFIG.weekendSurcharge : 0;
  const travel = zone.price;
  const total = survey + travel + weekendFee;

  return { survey, travel, weekendFee, total, valid: true };
}
