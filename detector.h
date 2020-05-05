#ifndef _detector_H
#define _detector_H

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

#define MAXL 1000005 

class Alg //基于hash表的算法
{   
  private:
    int L, LowerBound;
    double percentage;
    BOBHash32* bobhash;
    double var[MAXL];
    uint64_t id[MAXL];
    int counter[MAXL][2]; //low,high 
  
  public:
    int total;
    pair <pair <uint64_t, uint64_t>, double> ans[N];
    
    Alg(int t, int M, int L) : timeStamp(t, M), L(L), total(0) { //前两个参数是bloomfliter的参数, L是hash表的大小
        LowerBound = 8;
        percentage = 0.65;
        bobhash = new BOBHash32(rand() % 5000);
        memset(counter, 0, sizeof(counter));
    } 

    void insert(const uint64_t& s, const double& x) { //丢一个元素进来
        unsigned int H = bobhash -> run((char *)&s, 8);
        unsigned int pos = H % L;
        
        int l = cnt[pos][0];
        if (!l || id[pos] == s && check_near(var[pos], x)) {
            id[pos] = s;
            var[pos] = var[pos] * l / (l + 1) + x / (l + 1);
            ++cnt[pos][0];
        }
        else --cnt[pos][0];
        
        ++cnt[pos][1];

        return;
    }

    pair <bool, double> query(const uint64_t& s) { //询问某个id是否有周期性
        unsigned int H = bobhash -> run((char *)&s, 8);
        unsigned int pos = H % L;

        if (id[pos] != s) return make_pair(false, 0);
        if (cnt[pos][1] < LowerBound) return make_pair(false, 0);
        double p = (double)cnt[pos][0] / cnt[pos][1];
        return make_pair(p > percentage, var[s]);
    }
};
#endif
