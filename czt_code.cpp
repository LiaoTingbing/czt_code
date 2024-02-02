// czt_code.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <cmath>
#include <complex>
// #include<mat.h>
#include<string.h>
#include "fftw3.h"
#pragma comment(lib, "libfftw3-3.lib")

#define arrayLength( a ) sizeof( a ) / sizeof( a[0] )
#define arrayLengthPointer( a ) _msize( a ) / sizeof( a[0] )
#define PI 3.14159265358979323846
using namespace std; 



int main()
{
    /* 生成测试时间 */
    //std::cout << "Hello World!\n";
    double f1 = 1.765;
    double alpha1 = 0.005;
    double ampl1 = 1.3;
    double phase1 = 0.4;
    double f2 = 2.345;
    double alpha2 = 0.012;
    double ampl2 = 0.45;
    double phase2 = 1.234;
    // % # Time intervals of the signal :
    int Nsampling = 201;
    double* t_long = new double[Nsampling]();
    for (int i = 0; i < Nsampling; i++) {
        *(t_long + i) = 0.0 + i * (20.0 / (201 - 1)); 
      //  cout << *(t_long + i) << endl;
    }

    double* signal_long = new double[Nsampling]();
    for (int i = 0; i < Nsampling; i++) {
        double t = *(t_long + i);
        double s = ampl1 * cos(2 * PI * f1 * t - phase1) *exp(-alpha1 * t) + ampl2 * cos(2 * PI * f2 * t - phase2)*exp(-alpha2 * t);
        *(signal_long + i) = s ;
        //cout << *(signal_long + i) << endl;
    }

    double* s = signal_long;
    double* t = t_long;

    /*频率区间*/
    f1 = 1;
    f2 = 3;
    int Npoint = 5000;
    double* w = new double[Npoint]();
    for (int i = 0; i < Npoint; i++) {
        *(w + i) =( f1 + i * (f2 - f1) / (Npoint - 1)) *2 * PI;
        //cout << *(w + i) << endl;
    }
    /**/
    int N = int (  _msize(s) / sizeof(s[0]) ) ; 
    int M  =  int ( _msize(w) / sizeof(w[0])) ;
    double dt = *(t + 1) - *(t + 0); 
    double w1 = *(w + 0);
    double dw = *(w + 1) - *(w + 0);
    complex<double> A = exp(-1i * w1 * dt);
    complex<double> W = exp( + 1i * dw * dt);

    //cout << W << endl;
    int L =  int ( pow(2, ceil(log2(N + M + 1))) ) ; 
    
    /* -----------------------Y ------------------------------------*/
    complex<double>* y = new complex<double>[L](); 
    complex<double>* Y = new complex<double>[L]();
    for (int n = 0; n < L; n++) {
        if (n < N) {
            *(y + n) = pow(A, -n) * pow(W, pow(n, 2) / 2.0) * (*(s + n));
        }
    }

    fftw_complex* in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *L );
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * L );
    fftw_plan  p = fftw_plan_dft_1d(L, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (int i = 0; i < L ; i++) {
        in[i][0] = real(*( y+i)); // 实部
        in[i][1] = imag(*(y + i)); // 虚部
    }

    fftw_execute(p); /* repeat as needed */

    for (int i = 0; i < L; i++) {
        *(Y+i) = out[i][0] + 1i* out[i][1]; 
    }

    /*     -------------------------------------------------- v   -----------------------------   */

    complex<double>* v = new complex<double>[L]();
    complex<double>* V = new complex<double>[L]();
    for (int n = 0; n < L; n++) {
        if (n < M) {
            *(v +  n ) = pow(W, -pow(n, 2) / 2);
        }
        if ( (L - N + 1 <= n) && (n < L)  ) {
            *(v + n) = pow(W, -pow(L-n, 2) / 2);
        }
    }

    for (int i = 0; i < L; i++) {
        in[i][0] = real(*(v + i)); // 实部
        in[i][1] = imag(*(v + i)); // 虚部
    }

    fftw_execute(p); /* repeat as needed */

    for (int i = 0; i < L; i++) {
        *(V + i) = out[i][0] + 1i * out[i][1]; // 实部
    }
    
    /*------------计算G----------------------*/
    complex<double>* G = new complex<double>[L]();
    complex<double>* g = new complex<double>[L]();

    for (int i = 0; i < L; i++) {
        *(G + i) = *(V + i) * *(Y + i);
    }

    /* ------------------------------ ifft(G) -------------------------------*/

    p = fftw_plan_dft_1d(L, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    for (int i = 0; i < L; i++) {
        in[i][0] = real(*(G + i)); // 实部
        in[i][1] = imag(*(G + i)); // 虚部
    }

    fftw_execute(p); /* repeat as needed */

    for (int i = 0; i < L; i++) {
        *(g + i) =( out[i][0] + 1i * out[i][1] ) / complex<double>(L , 0)  ; 
    }

    fftw_destroy_plan(p); fftw_free(in); fftw_free(out);

    /*--------------   计算X 结果-------------------------*/
    complex<double>* X = new complex<double>[M]();

    for (int k = 0; k < M; k++) {
        *(X + k) = pow(W, pow(k, 2) / 2) * *(g + k); 
    }
    for (int i = 0; i < M; i++) {
        cout << *(X + i) << endl;
    }

    delete[] t_long , signal_long  , w ,y ,Y , V , v  , G  ,X ;

    // signal_long = ampl1 * cos(2 * pi * f1 * t_long - phase1).*exp(-alpha1 * t_long) + ampl2 * cos(2 * pi * f2 * t_long - phase2).*exp(-alpha2 * t_long);
}
