#ifndef _detectorHG_H
#define _detectorHG_H

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <string>
#include <cstring>
#include "BOBHash32.h"
#include "BOBHash64.h"

using namespace std;

#define T 20
#define MAXM 500005
#define MAXL 1000005 

class AlgHG //基于CMsketch的版本
{   
  private:
    int t, M, LowerBound;
    double percentage;
    BOBHash32* bobhash;

    double var[T * MAXM];
    uint64_t id[T * MAXM];
    int counter[T * MAXM]; //low
    int total[T * MAXM]; //high 
  
  public:
    AlgHG(int t, int M, int LowerBound, double percentage) : t(t), M(M), LowerBound(LowerBound), percentage(percentage) { 
        bobhash = new BOBHash32(123);
        memset(counter, 0, sizeof(counter));
        memset(total, 0, sizeof(total));
    } 
    
    bool check_near(double c, double x) {
        if (fabs(x - c) > Range * 2) return false;
        return fabs(x - c) < delta * max(fabs(c), 0.001);
    }

    void insert(const uint64_t& s, const double& x) { //丢一个元素进来
        unsigned long long int H = bobhash -> run((char *)&s, 8);
        unsigned long long int pos = H % M;
        
        int nl = pos * t, nr = (pos + 1) * t, cnt;

        for (int i = nl; i < nr; ++i) if ((cnt = counter[i]) && id[i] == s) {
            ++total[i];
            if (check_near(var[i], x)) {
                var[i] = var[i] * cnt / (cnt + 1) + x / (cnt + 1);
                ++counter[i];
            }
            else --counter[i];
          //  printf("#%.6lf cnt=%d total=%d var=%.6lf\n", x, counter[i], total[i], var[i]);
            return;
        } 
        
        int nx = nl;
        for (int i = nl; i < nr; ++i) if (counter[i] < counter[nx]) nx = i;

        ++total[nx];
        if (!counter[nx]) {
            counter[nx] = 1, var[nx] = x, id[nx] = s;
            return;
        }
        else --counter[nx];

        return;
    }

    pair <bool, double> query(const uint64_t& s) { //询问某个id是否有周期性
        unsigned long long int H = bobhash -> run((char *)&s, 8);
        unsigned long long int pos = H % M;
        
        int nl = pos * t, nr = (pos + 1) * t, cnt;
        
        for (int i = nl; i < nr; ++i) if ((cnt = counter[i]) && id[i] == s) {
            if (total[i] < LowerBound) continue;
            double p = (double)cnt / total[i];
            return make_pair(p >= percentage, var[i]);
        } 

        return make_pair(false, -1);
    }
};
#endif
