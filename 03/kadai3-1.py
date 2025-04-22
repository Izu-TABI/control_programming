import numpy as np
import matplotlib.pyplot as plt

# 直流サーボモータのパラメータ値
La = 3.00 * 10 ** -3
Ra = 2.26
Jm = 3.30 * 10 ** -5
Kv = 0.051
Kt = 0.051
Fv = 2.20 * 10 ** -5
Tc = 0.029
Ka = 3.00
As = 0.024
Ka = 3.00
Kp = 1.0




# 解くべき微分方程式1
def f_diadt(ia, om, omr):
  diadt = (Ka * As * Kp * omr - Ra * ia - om * (Kv + Ka * As * Kp)) / La

  return diadt

# 解くべき微分方程式2
def f_domdt(ia, om):
  domdt = (Kt * ia - Fv * om - Tc) / Jm

  return domdt

#ルンゲ・クッタ法
def RungeKutta(ia, om, va, h):
  ki_1 = f_diadt(ia, om, va)
  om_1 = f_domdt(ia, om)

  ki_2 = f_diadt(ia + h * ki_1 / 2, om + h * om_1 / 2, va)
  om_2 = f_domdt(ia + h * ki_1 / 2, om + h * om_1 / 2)

  ki_3 = f_diadt(ia + h * ki_2 / 2, om + h * om_2 / 2, va)
  om_3 = f_domdt(ia + h * ki_2 / 2, om + h * om_2 / 2)

  ki_4 = f_diadt(ia + h * ki_3, om + h * om_3, va)
  om_4 = f_domdt(ia + h * ki_3, om + h * om_3)

  ki = ia + h * (ki_1 + 2 * ki_2 + 2 * ki_3 + ki_4) / 6
  om = om + h * (om_1 + 2 * om_2 + 2 * om_3 + om_4) / 6

  return ki, om

# プロット
def plot_om(x, y1, y2, x_label, y_label):
  fig = plt.figure()

  # プロット準備
  sol = fig.add_subplot(1, 1, 1)
  sol.set_xlabel(x_label, fontsize=15, fontname='serif')
  sol.set_ylabel(y_label, fontsize=15, fontname='serif')
  sol.tick_params(axis='both', length=10, which='major', direction='in')
  sol.tick_params(axis='both', length=5, which='minor', direction='in')
  sol.minorticks_on()
  sol.grid(which='major')
  sol.plot(x, y1, 'r-', markersize=5)
  sol.plot(x, y2, 'b-', markersize=5)
  sol.legend(['h1(t)', 'h2(t)'])

  fig.tight_layout()
  plt.show()
  fig.clf()


# プロット
def plot_om(x, y1, y2, y3, x_label, y_label):
  fig = plt.figure()

  # プロット準備
  sol = fig.add_subplot(1, 1, 1)
  sol.set_xlabel(x_label, fontsize=15, fontname='serif')
  sol.set_ylabel(y_label, fontsize=15, fontname='serif')
  sol.tick_params(axis='both', length=10, which='major', direction='in')
  sol.tick_params(axis='both', length=5, which='minor', direction='in')
  sol.minorticks_on()
  sol.grid(which='major')
  sol.plot(x, y1, 'r-', markersize=5)
  sol.plot(x, y2, 'b-', markersize=5)
  sol.plot(x, y3, 'g-', markersize=5)
  sol.legend(['Wr=100[rad/s]', 'Wr=200[rad/s]', 'Wr=300[rad/s]'])

  fig.tight_layout()
  plt.show()


#メイン処理
if __name__ == "__main__":
  ia = [[0.0], [0.0], [0.0], [0.0]] # 初期値の設定
  om = [[0.0], [0.0], [0.0], [0.0]]
  omr = [100, 200, 300]

  h = 0.5 * 10 ** -3
  t_list = np.arange(0.0, 0.5, h)

  for i in range(3):
    for t in t_list:
      ia_t, om_t = RungeKutta(ia[i][-1], om[i][-1], omr[i], h) #yを求める
      ia[i].append(ia_t) #　求めたyをリストに追加
      om[i].append(om_t) #　求めたyをリストに追加

    ia[i].pop()
    om[i].pop()


  #プロット
  plot_om(t_list, om[0], om[1], om[2], "$t[s]$", "$ω(rad/s)$")
