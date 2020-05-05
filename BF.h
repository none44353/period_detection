#ifndef _BF_H
#define _BF_H

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
const double inf = 1e300;

class bloomfliter{
  private:
    int t, M; //t次hash 共t*M个位置
    BOBHash32* bobhash[T];
    double a[T * MAXM];
  public:
    bloomfliter(int t, int M) : t(t), M(M) {
        srand(time(NULL));
        for (int k = 0; k < t; ++k) 
            bobhash[k] = new BOBHash32(rand() % 5000);
        memset(a, 0, sizeof(a));
    }

    void insert(const uint64_t& x, const double& key) {
        for (int k = 0; k < t; ++k) {
            unsigned int H = bobhash[k] -> run((char *)&x, 8);
            unsigned int pos = H % M;

            a[k * M + pos] = key;
        }

        return;
    }

    double query(const uint64_t& x) {
        double res = inf;
        for (int k = 0; k < t; ++k) {
            unsigned int H = bobhash[k] -> run((char *)&x, 8);
            unsigned int pos = H % M;

            res = min(res, a[k * M + pos]);
        }

        return res;
    }
};

#endif
