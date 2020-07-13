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
    Alg(int L, int LowerBound, double percentage) : L(L), LowerBound(LowerBound), percentage(percentage) { 
        bobhash = new BOBHash32(123);
        memset(counter, 0, sizeof(counter));
    } 

    bool check_near(double c, double x) {
        if (fabs(x - c) > Range * 2) return false;
        return fabs(x - c) < delta * max(fabs(c), 0.001);
    }

    void insert(const uint64_t& s, const double& x) { //丢一个元素进来
       // printf("insert %.6lf\n", x);
        unsigned long long int H = bobhash -> run((char *)&s, 8);
        unsigned long long int pos = H % L;
        
        int l = counter[pos][0];
        if (!l || id[pos] == s && check_near(var[pos], x)) {
            id[pos] = s;
            var[pos] = var[pos] * l / (l + 1) + x / (l + 1);
            ++counter[pos][0];
        }
        else --counter[pos][0];

        ++counter[pos][1];

        return;
    }

    pair <bool, double> query(const uint64_t& s) { //询问某个id是否有周期性
        unsigned long long int H = bobhash -> run((char *)&s, 8);
        unsigned long long int pos = H % L;

        if (id[pos] != s) return make_pair(false, 0);
        if (counter[pos][1] < LowerBound) return make_pair(false, 0);
        double p = (double)counter[pos][0] / counter[pos][1];
        return make_pair(p >= percentage, var[pos]);
    }
};
#endif
