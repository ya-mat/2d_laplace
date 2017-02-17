# 2d_laplace

## 内容

2次元Laplace方程式のコードです．`lp_laplace.f90`が本体のようなもので，一重層ポテンシャルと二重層ポテンシャルを計算するための関数があります．main.f90では，正円を区分一定要素を用いた選点法で離散化し，正解が得られる境界値問題を解いて，数値解と正解の相対誤差を比較しています．

## 必要なもの
- fortranコンパイラ
- blas
- lapack

## 使い方
- Makefileを自分の環境に合うように編集する．
    - $(LINKER) fortranコンパイラを指定
    - $(FORTRAN) fortranコンパイラを指定
    - $(LDFLAGS) blas, lapackのpathを指定
- シェルで次のコマンドを打つ
```
make
./a.out
```
- inputの中身を編集することで，未知数の数と円の半径を変更できます．
