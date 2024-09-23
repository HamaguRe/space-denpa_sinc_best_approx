//! テレスコーピング法でチェビシェフ多項式の特徴を持った近似式を作る
//! 
//! 偶関数の場合には奇数次の項がゼロになるので，N=8かつCUT_TERMS=2とすれば
//! テレスコーピング法で得られる近似式は4次になる．

use std::fs::File;
use std::io::{Write, BufWriter};

/// グラフ描画の際の計算分割数
const SPLIT_NUM: usize = 2000;

/// N次までの項を計算する（ゼロから数えるので，項の総数はN+1）
const N: usize = 8;  // >= 2

/// テレスコーピング法で切り落とす項の数
const CUT_TERMS: usize = 2;

fn main() {
    // テレスコーピング法で係数を求める
    let s = sinc_telescoping(CUT_TERMS);

    // 計算した係数を表示
    println!("let s = [");
    for i in 0..(N+1) {
        println!("    {:?},  // s{}", s[i], i);
    }
    println!("];");

    // --- ファイルに書き出し --- //
    // 計算結果の保存先（同一ファイルが存在したら上書き）
    let mut file = BufWriter::new(File::create("./result.csv").unwrap());

    let h = 2.0 / (SPLIT_NUM as f64);  // 計算区間が[-1, 1]の場合
    for i in 0..(SPLIT_NUM + 1) {
        let x = h * (i as f64) - 1.0;

        // マクローリン展開による近似式
        let y_maclaurin = sinc_maclaurin(x, 6);

        // チェビシェフ近似式
        let y_chebyshev = polynominal(x, s);

        let y_true = x.sin() / x;
        let error_maclaurin = y_maclaurin - y_true;
        let error_chebyshev = y_chebyshev - y_true;
        file.write(format!(
            "{},{},{},{},{},{}\n",
            x, y_true, y_maclaurin, y_chebyshev, error_maclaurin, error_chebyshev
        ).as_bytes()).unwrap();
    }
}

/// 多項式を計算
/// 
/// * s: [s0, s1, s2, ..., sN]
fn polynominal(x: f64, s: [f64; N+1]) -> f64 {
    let mut s = s;
    s.reverse();  // ホーナー法で計算するために順序を反転する

    let mut y = s[0];
    for i in 1..N {
        y = y * x + s[i];
    }
    y = y * x + (s[N] - 1.0);  // 精度改善のノウハウ
    y + 1.0
}

/// n次の項までのマクローリン展開でsinc(x)を計算する（ただし n >= 1）
fn sinc_maclaurin(x: f64, n: i32) -> f64 {
    let z = x * x;
    let mut y = 1.0;
    let mut numer = 1.0;
    let mut denom = 1.0;
    let mut sign_positive = false;
    for i in 1..=n {
        denom *= (i + 1) as f64;
        if (i & 1) == 0 {  // iが偶数
            numer *= z;
            if sign_positive {
                y += numer / denom;
            } else {
                y -= numer / denom;
            }
            sign_positive = !sign_positive;
        }
    }
    y
}

/// 階乗
/// 
/// n!
fn factorial(n: i32) -> f64 {
    assert!(n >= 0);
    let mut y = 1.0;
    for i in 1..=n {
        y *= i as f64;
    }
    y
}

/// 二項係数（binomial coefficients）
/// 
/// nCk
fn binomial(n: i32, k: i32) -> f64 {
    assert!(n >= k);
    factorial(n) / ( factorial(k) * factorial(n - k) )
}

/// 対象とする関数f(x)の近似多項式の係数列aを
/// チェビシェフ多項式Tn(x)の係数列bに変換する．
/// 
/// 任意のf(x)を表す多項式に利用可能（係数しか渡していないため）
/// 
/// * a: [a0, a1, a2, ..., aN]
fn convert_exp2chebyshev(a: [f64; N+1]) -> [f64; N+1] {
    // 多項式の係数をまとめた配列
    // x^0 =   1*T_0(x)
    // x^1 =               1*T_1(x)
    // x^2 = 0.5*T_0(x)               0.5*T_2(x)
    // x^3 =             0.75*T_1(x)              0.25*T_3(x)
    // ...
    let mut coefs = [[0.0; N+1]; N+1];

    coefs[0][0] = 1.0;
    coefs[1][1] = 1.0;
    for n in 2..(N+1) {
        let tmp = 2.0f64.powi(n as i32 - 1);
        for k in 0..=(n / 2) {
            // x^nとT_n(x)の関係式を実装する
            let index = n - 2*k;
            coefs[n][index] = binomial(n as i32, k as i32) / tmp;
            if index == 0 {
                coefs[n][index] *= 0.5;
            }
        }
    }

    // coefsの行ごとにaをかけて同じ列の要素を全て足せばbが求まる
    let mut b = [0.0; N+1];
    for i in 0..N {
        for j in i..N {
            b[i] += a[j] * coefs[j][i];
        }
    }

    b  // b[0]*T_0(x) + b[1]*T_1(x) + b[2]*T_2(x) + b[3]*T_3(x) + ...
}

/// チェビシェフ多項式Tn(x)の係数列bからx^nの多項式係数列aを求める
/// 
/// 任意のf(x)を表す多項式に利用可能（係数しか渡していないため）
/// 
/// * b: [b0, b1, b2, ..., bN]
fn convert_chebyshev2exp(b: [f64; N+1]) -> [f64; N+1] {
    // 多項式の係数をまとめた配列
    // T0(x) =  1
    // T1(x) =       1*x^1
    // T2(x) = -1            2*x^2
    // T3(x) =      -3*x^1           4*x^3
    // ...
    let mut coefs = [[0.0; N+1]; N+1];

    coefs[0][0] = 1.0;
    coefs[1][1] = 1.0;
    for i in 2..(N+1) {
        for j in 0..i {
            // T_n+1(x) = 2*x*T_n(x) - T_n-1(x)
            coefs[i][j+1] = 2.0 * coefs[i-1][j];
            coefs[i][j] -= coefs[i-2][j];
        }
    }

    // coefsの行ごとにbをかけて同じ列の要素を全て足せばsが求まる
    let mut a = [0.0; N+1];
    for i in 0..(N+1) {
        for j in i..(N+1) {
            a[i] += b[j] * coefs[j][i];
        }
    }

    a  // a[0] + a[1]*x + a[2]*x^2 + a[3]*x^3 + ...
}

/// sin(x)/xをマクローリン展開した際の係数列から，
/// テレスコーピング法でチェビシェフ多項式の係数を求める．
/// 
/// * cut_terms: 縮小する項数
fn sinc_telescoping(cut_terms: usize) -> [f64; N+1] {
    assert!(2 * cut_terms < N);

    // sin(x)/xをマクローリン展開した場合の，各項の係数を格納（値がゼロの係数を含む）
    let mut a = [0.0; N+1];
    a[0] = 1.0;
    let mut sign_positive = false;
    for n in 2..(N+1) {
        if (n & 1) == 0 {  // nが偶数のときだけ計算
            let tmp = 1.0 / factorial(n as i32 + 1);
            if sign_positive {
                a[n] = tmp;
            } else {
                a[n] = -tmp;
            }
            sign_positive = !sign_positive;
        }
    }

    let mut b = convert_exp2chebyshev(a);
    // 最高次数から指定した分だけ削る
    if (N & 1) == 0 {  // 最高次数が偶数次の場合
        println!("zero n");
        for i in 0..cut_terms {
            b[N - 2 * i] = 0.0;
        }
    } else {  // 最高次数が奇数次なら二次（偶数次の項まで）削る
        println!("zero n-1");
        for i in 0..cut_terms {
            b[N - 2*i - 1] = 0.0;
        }
    }

    convert_chebyshev2exp(b)
}