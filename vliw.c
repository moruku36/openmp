// OpenMP header
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

### VLIW Quiz
次の命令列を3演算器を持つVLIWシステムに最適化したい．命令列を最適となるように並び替えよ．
また，同時実行される命令の組み合わせを示せ．
条件分岐命令も一つの命令である．但し，トレースはif文を満たす命令列であると仮定する．
命令列は先頭から開始される．


m01; m02; m03; m04;
if( m05 ){ m11; m12; m13; m14; m15; }
else{ m21; m22; }
m31; m32; m33;
if( m34 ){ m41; m42; m43; }
else{ m51; m52; m53; m54; }
m61; m62; m63;
