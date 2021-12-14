### OpenMP | Introduction with Installation Guide
https://www.geeksforgeeks.org/openmp-introduction-with-installation-guide/

### Commands
export OMP_NUM_THREADS=8

Compile: 
gcc geeksforgeeks.c -o gfg

Execute:
./gfg

### Macでコンパイルする時

clang -Xpreprocessor -fopenmp -lomp -I"$(brew --prefix libomp)/include" -L"$(brew --prefix libomp)/lib" ファイル名.c
export OMP_NUM_THREADS=8

### Parallel Quiz
lmpcc にloginし，行列の乗算と行列式を求めるプログラム (/home/lecture/is/i416/MMMUL.c )を並列化せよ．

a) 並列化版を完成し，行列乗算部分の実行時間を確認せよ．並列度(実行ノード数)を変化させ，並列度に対する実行時間(あるいは高速化率)をグラフに示してレポートに添付せよ．

b) 行列乗算部分は，何回の浮動小数点演算を行っているか？　(ループ計算のための演算などは含まない)
lmpccの理論演算性能に対して，並列化したプログラムは何％の性能を得られたか？　並列度やコア数も考慮し計算せよ．
c) 得られたプログラムは，強いスケーリングか，弱いスケーリングか?　測定による根拠を示しながら説明せよ．

d) 並列化効率を向上させよ．配列の分割(ブロック v.s. サイクリック)や，どのループを並列化させるか(ompディレクティブの挿入場所)，またディレクティブ・オプションの利用がヒントである．(授業で触れていないオプションを使っても良い)

e) プログラムの並列化のポイント，および測定結果とその考察を，5分程度のプレゼンテーションにまとめ，発表せよ．(11月28日の授業時に発表)


### OpenMP Quiz
lmpcc にloginし，モンテカルロ法による円周率計算の例題プログラムを並列化し完成させよ．
プログラムは /home/lecture/is/i416/omp_CYCLIC.c に逐次版が用意してある．

モンテカルロ法による円周率計算の原理:

区間[0-1]の一様乱数を多数回発生させ，  に代入する．全体の発生数に対して， を満たす割合が である．発生回数が増加すると，精度が向上するが，精度向上には多大な処理時間を要する．複数のノードで分割し乱数発生と円内か否かの判定を行うことによって，処理時間を短縮することができる．本アルゴリズムはそれぞれの事象が独立であるので，並列化が容易である．

a) 雛形を参考に並列化版を完成させ，実行を確認せよ．プログラムと実行結果をレポートに記せよ．
b) 実行ノード数を変化させ，実行時間を比較せよ．実行時間のグラフを書きレポートに添付せよ．

c) 並列度に対する速度向上率から，アムダールの法則における並列化可能部分と並列化不可能部分の比率を計算せよ．

d) rand48関数(drand48, srand48など関係する乱数関数)を工夫して，高速化を達成せよ．

以上の実装のポイントと結果を5分程度のプレゼンテーションにまとめ，発表せよ．(11月28日の授業時)


### kagayaki にloginし，行列の乗算の例題プログラムを並列化し完成させよ．
プログラムは /home/lecture/is/i416/mpi_MMMUL.cにひな形が用意してある．


a) 並列化版を完成し，実行時間を確認せよ．並列化版の総和は各自考えよ．(ヒント:授業で説明した関数をうまく使うと良い)

b) 実行ノード数を変化させ，実行時間の推移をOpenMP版と比較しながらグラフに示せ．

c) プログラムの並列化のポイント，および測定結果とその考察を，5分程度のプレゼンテーションにまとめ，発表せよ．(12月18日の授業時に発表)


### MPIによる円周率計算の並列化

kagayaki にloginし，モンテカルロ法による円周率計算の例題プログラムを並列化し完成させよ．
プログラムは /home/lecture/is/i416/mpi_PIrep1.c にひな形が用意してある．

a) 雛形を参考に並列化版を完成させ，実行を確認せよ．
b) 実行ノード数を変化させ，実行時間の推移をOpenMP版と比較しながらグラフに示せ．
c) 実行ノード数を変化させたときの誤差(計算されたPIと真のPIの差)の変化をグラフに示せ。
以上について実装のポイントと結果を5分程度のプレゼンテーションにまとめ，発表せよ．(12月18日の授業時)