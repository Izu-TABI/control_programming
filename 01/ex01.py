import numpy as np
import matplotlib.pyplot as plt

# 解くべき微分方程式
def f_dydt( y ):
  a = -1.0
  dydt = a * y
  return dydt


#ルンゲ・クッタ法
def RungeKutta(y,h):
  kd1 = f_dydt(y)
  kd2 = f_dydt(y + h * kd1 / 2)
  kd3 = f_dydt(y + h * kd2 / 2)
  kd4 = f_dydt(y + h + kd3)
  ky = y + h * (kd1 + 2 * kd2 + 2 * kd3 + kd4) / 6
  return ky




# プロット

def plot(x, y, x_label, y_label):
  fig = plt.figure()

  # プロット準備
  sol = fig.add_subplot(1, 1, 1)
  sol.set_xlabel(x_label, fontsize=20, fontname='serif')
  sol.set_ylabel(y_label, fontsize=20, fontname='serif')

  sol.tick_params(axis='both', length=10, which='major', direction='in')
  sol.tick_params(axis='both', length=5, which='major', direction='in')
  sol.minorticks_on() # 補助メモリをつける
  sol.grid(which='major')
  sol.plot(x, y, 'b-', markersize=5)

  #スクリーン表示
  fig.tight_layout()
  plt.show()
  fig.clf()


#メイン処理
if __name__ == "__main__":
  y = [1.0] # 初期値の設定
  h = 0.01
  t_list = np.arange(0.0, 10.0, h)

  for t in t_list:
    yt = RungeKutta(y[-1], h) #yを求める
    y.append(yt) #　求めたyをリストに追加

  y.pop() # リストの末尾の要素を削除 要素数1001から1000

  #プロット
  plot(t_list, y, "$t$", "$y(t)$")
