#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef complex<double> cx;
typedef vector<double> vd;

void fft(vector<cx>& a){
    int n = a.size();
    int L = 31 - __builtin_clz(n);
    static vector<complex<long double>>R(2, 1);
    static vector<cx>rt(2, 1);
    for(int k = 2; k < n; k <<= 1){
        R.resize(n);
        rt.resize(n);
        auto x = polar(1.0L, acos(-1.0L) / k);
        for(int i = k; i < 2 * k; i++){
            rt[i] = R[i] = (i & 1 ? R[i >> 1] * x : R[i >> 1]);
        }
    }

    vector<int>rev(n);
    for(int i = 0; i < n; i++){
        rev[i] = (rev[i >> 1] | (i & 1) << L) >> 1;
        if(i < rev[i]){
            swap(a[i], a[rev[i]]);
        }
    }

    for(int k = 1; k < n; k <<= 1){
        for(int i = 0; i < n; i += 2 * k){
            for(int j = 0; j < k; j++){
                auto x = (double *)& rt[j + k];
                auto y = (double *)& a[i + j + k];
                cx z(x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0]);
                a[i + j + k] = a[i + j] - z;
                a[i + j] += z;
            }
        }
    }
}

template<int M>
vector<int> multiply(const vector<int>& a, const vector<int>& b){
    if(a.empty() || b.empty())return {};
    int sz1 = a.size(), sz2 = b.size();
    vector<int>res(sz1 + sz2 - 1);
    int B = 32 - __builtin_clz(sz1 + sz2 - 1);
    int n = 1 << B, cut = int(sqrt(M));
    
    vector<cx>L(n), R(n), outs(n), outl(n);
    for(int i = 0; i < sz1; i++){
        L[i] = cx((int)a[i] / cut, (int)a[i] % cut);
    }

    for(int i = 0; i < sz2; i++){
        R[i] = cx((int)b[i] / cut, (int)b[i] % cut);
    }

    fft(L), fft(R);
    for(int i = 0; i < n; i++){
        int j = -i & (n - 1);
        outl[j] = (L[i] + conj(L[j])) * R[i] / (2.0 * n);
        outs[j] = (L[i] - conj(L[j])) * R[i] / (2.0 * n) / 1i;
    }

    fft(outl), fft(outs);
    for(int i = 0; i < sz1 + sz2 - 1; i++){
        ll av = ll(real(outl[i]) + .5), cv = ll(imag(outs[i]) + .5);
        ll bv = ll(imag(outl[i]) + .5) + ll(real(outs[i]) + .5);
        res[i] = ((av % M * cut + bv) % M * cut + cv) % M;
    }
    return res;
}

int main(){
    cin.tie(nullptr)->sync_with_stdio(false);
    int n, m;
    cin >> n >> m;
    
    vector<int>a(n);
    for(int i = 0; i < n; i++){
        cin >> a[i];
    }
    
    vector<int>b(m);
    for(int j = 0; j < m; j++){
        cin >> b[j];
    }
    
    auto res = multiply<int(1e9) + 7>(a, b);
    for(int k = 0; k < n + m - 1; k++){
        cout << res[k] << " ";
    }
        
    return 0;
}
