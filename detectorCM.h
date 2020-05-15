#ifndef _detectorCM_H
#define _detectorCM_H

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

#define T 6
#define MAXM 500005
#define MAXL 1000005 
#define delta 0.01

class AlgCM //基于CMsketch的版本
{   
  private:
    int t, M, LowerBound;
    double percentage;
    BOBHash32* bobhash[T];

    double var[T * MAXM];
    uint64_t id[T * MAXM];
    int counter[T * MAXM][2]; //low,high 
  
  public:
    Alg(int t, int M, int LowerBound, double percentage) : t(t), M(M), LowerBound(LowerBound), percentage(percentage) { 
        for (int k = 0; k < t; ++k)
            bobhash[k] = new BOBHash32(rand() % 2000);
        memset(counter, 0, sizeof(counter));
    } 

    bool check_near(double c, double x) {
        return fabs(x - c) < delta;
    }

    void insert(const uint64_t& s, const double& x) { //丢一个元素进来
        for (int k = 0; k < t; ++k) {
            unsigned long long int H = bobhash[k] -> run((char *)&s, 8);
            unsigned long long int pos = k * M + H % M;
            
            int l = counter[pos][0];
            if (!l || id[pos] == s && check_near(var[pos], x)) {
                id[pos] = s;
                var[pos] = var[pos] * l / (l + 1) + x / (l + 1);
                ++counter[pos][0];
            }
            else --counter[pos][0];
            
            ++counter[pos][1];
        }

        return;
    }

    pair <bool, double> query(const uint64_t& s) { //询问某个id是否有周期性
        double ans = 0, p;
        int l = 0, h = 2e9;
        for (int k = 0; k < t; ++k) {
            unsigned long long int H = bobhash[k] -> run((char *)&s, 8);
            unsigned long long int pos = k * M + H % M;

            h = min(h, counter[pos][1]);
            if (id[pos] != s) continue;
            if (l < counter[pos][0])
                ans = var[pos], l = counter[pos][0];
        }

        if (counter[pos][1] < LowerBound) return make_pair(false, 0);

        p = (double)l / h;
        return make_pair(p >= percentage, ans);
    }
};
#endif
