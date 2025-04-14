import numpy as np
import matplotlib.pyplot as plt

# 解くべき微分方程式1
def f_dh1dt(h1, h2):
  C1 = 5.0
  C2 = 4.0
  R1 = 2.0
  qi = 2.0

  dh1dt = (qi - ((h1 - h2) / R1)) / C1

  return dh1dt

# 解くべき微分方程式2
def f_dh2dt(h1, h2):
  C1 = 5.0
  C2 = 4.0
  R1 = 2.0
  R2 = 10.0

  dh2dt = (((h1 - h2) / R1) - h2 / R2) / C2

  return dh2dt

#ルンゲ・クッタ法
def RungeKutta(h1, h2, h):
  kh1_1 = f_dh1dt(h1, h2)
  kh2_1 = f_dh2dt(h1, h2)

  kh1_2 = f_dh1dt(h1 + h * kh1_1 / 2, h2 + h * kh2_1 / 2)
  kh2_2 = f_dh2dt(h1 + h * kh1_1 / 2, h2 + h * kh2_1 / 2)

  kh1_3 = f_dh1dt(h1 + h * kh1_2 / 2, h2 + h * kh2_2 / 2)
  kh2_3 = f_dh2dt(h1 + h * kh1_2 / 2, h2 + h * kh2_2 / 2)

  kh1_4 = f_dh1dt(h1 + h * kh1_3, h2 + h * kh2_3)
  kh2_4 = f_dh2dt(h1 + h * kh1_3, h2 + h * kh2_3)

  kh1 = h1 + h * (kh1_1 + 2 * kh1_2 + 2 * kh1_3 + kh1_4) / 6
  kh2 = h2 + h * (kh2_1 + 2 * kh2_2 + 2 * kh2_3 + kh2_4) / 6

  return kh1, kh2


# プロット
def plot(x, y1, y2, x_label, y_label):
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




#メイン処理
if __name__ == "__main__":
  h1 = [0.0] # 初期値の設定
  h2 = [0.0]
  h = 0.1
  t_list = np.arange(0.0, 500.0, h)

  for t in t_list:
    h1_t, h2_t = RungeKutta(h1[-1], h2[-1], h) #yを求める
    h1.append(h1_t) #　求めたyをリストに追加
    h2.append(h2_t) #　求めたyをリストに追加

  h1.pop()
  h2.pop()

  #プロット
  plot(t_list, h1, h2, "$t[s]$", "$h1(t),h2(t)$")
