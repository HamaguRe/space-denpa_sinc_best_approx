# sinc関数の最良近似式を作る

宇宙電波実験室の記事内で使用したプログラムです。
実装する際の参考にしてください。

当該記事：[sinc関数の最良近似式を作る（宇宙電波実験室）](https://space-denpa.jp/2024/09/23/sinc-best-approx/)

## 内容

各ディレクトリにRustで記述した計算プログラムと、グラフ表示用のPythonスクリプトが入っています。

* /telescoping：テレスコーピング法でチェビシェフ多項式の特徴を持った近似式を作るプログラム
* /refine：テレスコーピング法で作った近似式を元にして最良近似式を作るプログラム

## 使い方

各ディレクトリに移動して以下のコマンドを実行してください。

```
$ cargo run && python3 data_plot.py
```