# result.csvを読んでグラフ描画

import csv
import matplotlib.pyplot as plt

x = []
y_true = []
y_maclaurin = []
y_chebyshev = []
error_maclaurin = []
error_chebyshev = []

# CSVからデータを読み出して配列に追加
with open('./result.csv') as file:
    reader = csv.reader(file)
    for row in reader:
        nums = [float(v) for v in row]
        x.append(nums[0])
        y_true.append(nums[1])
        y_maclaurin.append(nums[2])
        y_chebyshev.append(nums[3])
        error_maclaurin.append(nums[4])
        error_chebyshev.append(nums[5])


# --- 描画 --- #
# 近似曲線
fig1 = plt.figure(figsize = (10, 7))  # 横, 縦
ax1 = fig1.add_subplot(211)
ax1.set_title('Approximate curve', fontsize=20)
ax1.step(x, y_true,      label="true",      linestyle="dashdot")
ax1.step(x, y_maclaurin,    label="maclaurin",    linestyle="-", color="green")
ax1.step(x, y_chebyshev, label="chebyshev", linestyle="-", color="orange")
ax1.set_xlim(x[0], x[-1])
ax1.tick_params(labelsize=13)  # 軸目盛の大きさ
ax1.legend(fontsize=15, loc="upper right")
ax1.set_ylabel("sinc(x)", fontsize=15)

# 誤差曲線
ax2 = fig1.add_subplot(212)
ax2.step(x, error_maclaurin,    label="maclaurin",    color="green")
ax2.step(x, error_chebyshev, label="chebyshev", color="orange")
ax2.set_xlim(x[0], x[-1])
ax2.tick_params(labelsize=13)
ax2.legend(fontsize=15, loc="upper right")
ax2.set_ylabel("error",   fontsize=15)
ax2.set_xlabel("x [rad]", fontsize=15)

plt.show()