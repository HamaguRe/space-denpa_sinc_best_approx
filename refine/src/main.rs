//! テレスコーピング法で求めた係数を改良してsinc関数の最良近似式を求める

use std::fs::File;
use std::io::{Write, BufWriter};
use nalgebra::{SMatrix, SVector};

use std::f64::consts::PI;

/// N次までの項を計算する（ゼロから数えるので，項の総数はN+1）
/// 
/// 偏差点はN+2個必要で，偏差点の個数は変えられないので最高次数の方を合わせる．
/// 最高次数の係数がゼロになってしまっても問題ない．
const N: usize = 5;

/// グラフ描画の際の分割数
const SPLIT_NUM: usize = 2000;

fn main() {
    // 最適化する近似多項式の係数 [c0, c1, c2, ..., cN]
    let c = SVector::<f64, {N+1}>::from_iterator([
        0.9999937996031747,  // s0
        0.0,  // s1
        -0.16655505952380953,  // s2
        0.0,  // s3
        0.008035714285714285,  // s4
        0.0,  // s5
    ].iter().cloned());

    println!("c (origin): {:?}", c);

    // この係数を改良する
    let mut c_fined = c;

    // xjを配置
    let mut xj = spread_xj([-1.0, 1.0]);

    // 反復計算
    for _ in 0..10 {
        // xjを更新
        let er_closure = |x: f64| er(x, c_fined);  // クロージャを使用してerの引数aとdを固定
        for i in 1..(xj.len() - 1) {
            xj[i] = refine_xj(&er_closure, xj[i]);
        }

        let mut mat = SMatrix::<f64, {N+1}, {N+1}>::zeros();
        let mut e_vec = SVector::<f64, {N+1}>::zeros();
        for i in 0..(N+1) {
            for j in 0..(N+1) {
                mat[(i, j)] = xj[i].powi(j as i32) + xj[i+1].powi(j as i32);
            }
            e_vec[i] = -(er(xj[i], c_fined) + er(xj[i+1], c_fined));
        }

        // 連立方程式を解く
        let delta_c = mat.lu().solve(&e_vec).unwrap();
        c_fined = c_fined + delta_c;

        // 本来0になるはずの係数を強制的に0にする
        for i in 0..=N {
            if (i & 1) == 1 {  // 奇数次の項
                c_fined[i] = 0.0;
            }
        }

        // 収束判定（粒の揃い方を調べる）
        let mut eps_max: f64 = 0.0;
        for i in 0..xj.len() {
            eps_max = eps_max.max( er(xj[i], c_fined).abs() );
        }
        let mut eps_min = 0.0;
        for i in 0..xj.len() {
            eps_min = eps_max.max( er(xj[i], c_fined).abs() / eps_max );
        }
        println!("eps_min: {}", eps_min);
        if eps_min > 0.999999 {  // 9の個数は，厳しくなければ2，少し厳しくしたければ6程度．
            break;
        }
    }

    println!("c (fined):  {:?}", c_fined);

    // 最後にxjを更新（グラフ描画用）
    let er_closure = |x: f64| er(x, c_fined);
    for i in 1..(xj.len() - 1) {
        xj[i] = refine_xj(&er_closure, xj[i]);
    }

    // --- ファイルに書き出し --- //
    // 計算結果の保存先（同一ファイルが存在したら上書き）
    let mut file = BufWriter::new(File::create("./result.csv").unwrap());
    let mut file_xj = BufWriter::new(File::create("./points.csv").unwrap());

    let h = 2.0 / (SPLIT_NUM as f64);  // 計算区間が[-1, 1]の場合
    for i in 0..(SPLIT_NUM+1) {
        let x = h * (i as f64) - 1.0;

        let y_true = sinc_true(x);
        let y_orig = polynominal(x, c);
        let y_fined = polynominal(x, c_fined);
        let error_orig = er(x, c);
        let error_fined = er(x, c_fined);
        file.write(format!(
            "{},{},{},{},{},{}\n",
            x, y_true, y_orig, y_fined, error_orig, error_fined
        ).as_bytes()).unwrap();
    }

    // 偏差点の位置を保存
    for i in 0..xj.len() {
        file_xj.write(format!("{},{}\n", xj[i], er(xj[i], c_fined)).as_bytes()).unwrap();
    }
}

/// sinc関数の真値として扱う
fn sinc_true(x: f64) -> f64 {
    if x == 0.0 {
        1.0
    } else {
        x.sin() / x
    }
}

/// ベクトルの順序を前後入れ替える（ホーナー法で計算するために使う）
fn reverse(a: &mut SVector<f64, {N+1}>) {
    for i in 0..((N + 1) / 2) {
        let tmp1 = a[i];
        let tmp2 = a[N - i];
        a[i] = tmp2;
        a[N - i] = tmp1;
    }
}

/// 多項式を計算
/// 
/// * s: `[s0, s1, s2, ..., sN]`
/// * return: `s0 + s0*x + s1*x^2 + s2*x^3 + ... + sN*x^N`
fn polynominal(x: f64, s: SVector<f64, {N+1}>) -> f64 {
    let mut s = s;
    reverse(&mut s);   // ホーナー法で計算するために順序を反転する

    let mut y = s[0];
    for i in 1..N {  // 一番最後を残して計算
        y = y * x + s[i];
    }
    y = y * x + (s[N] - 1.0);  // 精度改善のノウハウ
    y + 1.0
    
}

/// 誤差関数（絶対誤差基準）
/// 
/// c: 近似式の係数 [c0, c1, c2, ..., cN]
fn er(x: f64, c: SVector<f64, {N+1}>) -> f64 {
    polynominal(x, c) - sinc_true(x)
}

/// ニュートン法で偏差点を改良する
/// 
/// * er(x): 誤差関数
/// * xj: 偏差点のx座標
/// 
/// 参考：近似式のプログラミング, p.69
fn refine_xj<F>(er: &F, xj: f64) -> f64
where F: Fn(f64) -> f64 {
    let mut xj = xj;

    // 一つの山を8分割して偏差点の両側4点を使う
    let mut h = 0.25 / (N as f64);
    for _ in 0..4 {
        // 4次精度の微分
        let s = er(xj + 2.0*h);
        let t = er(xj + h);
        let u = er(xj);
        let v = er(xj - h);
        let w = er(xj - 2.0*h);
        xj -= h * ((s - w) - 8.0*(t - v)) / ((s + w) - 16.0*(t + v) + 30.0*u);

        h *= 0.5;
    }

    xj
}

/// 偏差点を配置する
fn spread_xj(x_range: [f64; 2]) -> [f64; N+2] {
    assert!(x_range[0] < x_range[1]);

    let r = (x_range[1] - x_range[0]) * 0.5;
    let mut xj = [0.0; N+2];  // 偏差点の座標
    for i in 0..(N+2) {
        let p = (i as f64) * PI / ((N + 1) as f64);
        xj[i] = r * p.cos();
    }

    xj
}
