#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <time.h>
#include <unordered_map>
#include <list>
using namespace std;

#include "ssummary.h"
#include "BF.h"
#include "detector.h"
                
string datapath[60] = {"/usr/share/dataset/CAIDA2018/dataset/130000.dat",
                "/usr/share/dataset/CAIDA2018/dataset/130100.dat",
                "/usr/share/dataset/CAIDA2018/dataset/130200.dat",
                "/usr/share/dataset/CAIDA2018/dataset/130300.dat",
                "/usr/share/dataset/CAIDA2018/dataset/130400.dat",
                "/usr/share/dataset/CAIDA2018/dataset/130500.dat",
                "/usr/share/dataset/CAIDA2018/dataset/130600.dat",
                "/usr/share/dataset/CAIDA2018/dataset/130700.dat",
                "/usr/share/dataset/CAIDA2018/dataset/130800.dat",
                "/usr/share/dataset/CAIDA2018/dataset/130900.dat",
                "/usr/share/dataset/CAIDA2018/dataset/131000.dat",
                "/usr/share/dataset/CAIDA2018/dataset/131100.dat",
                "/usr/share/dataset/CAIDA2018/dataset/131200.dat",
                "/usr/share/dataset/CAIDA2018/dataset/131300.dat",
                "/usr/share/dataset/CAIDA2018/dataset/131400.dat",
                "/usr/share/dataset/CAIDA2018/dataset/131500.dat",
                "/usr/share/dataset/CAIDA2018/dataset/131600.dat",
                "/usr/share/dataset/CAIDA2018/dataset/131700.dat",
                "/usr/share/dataset/CAIDA2018/dataset/131800.dat",
                "/usr/share/dataset/CAIDA2018/dataset/131900.dat",
                "/usr/share/dataset/CAIDA2018/dataset/132000.dat",
                "/usr/share/dataset/CAIDA2018/dataset/132100.dat",
                "/usr/share/dataset/CAIDA2018/dataset/132200.dat",
                "/usr/share/dataset/CAIDA2018/dataset/132300.dat",
                "/usr/share/dataset/CAIDA2018/dataset/132400.dat",
                "/usr/share/dataset/CAIDA2018/dataset/132500.dat",
                "/usr/share/dataset/CAIDA2018/dataset/132600.dat",
                "/usr/share/dataset/CAIDA2018/dataset/132700.dat",
                "/usr/share/dataset/CAIDA2018/dataset/132800.dat",
                "/usr/share/dataset/CAIDA2018/dataset/132900.dat",
                "/usr/share/dataset/CAIDA2018/dataset/133000.dat",
                "/usr/share/dataset/CAIDA2018/dataset/133100.dat",
                "/usr/share/dataset/CAIDA2018/dataset/133200.dat",
                "/usr/share/dataset/CAIDA2018/dataset/133300.dat",
                "/usr/share/dataset/CAIDA2018/dataset/133400.dat",
                "/usr/share/dataset/CAIDA2018/dataset/133500.dat",
                "/usr/share/dataset/CAIDA2018/dataset/133600.dat",
                "/usr/share/dataset/CAIDA2018/dataset/133700.dat",
                "/usr/share/dataset/CAIDA2018/dataset/133800.dat",
                "/usr/share/dataset/CAIDA2018/dataset/133900.dat",
                "/usr/share/dataset/CAIDA2018/dataset/134000.dat",
                "/usr/share/dataset/CAIDA2018/dataset/134100.dat",
                "/usr/share/dataset/CAIDA2018/dataset/134200.dat",
                "/usr/share/dataset/CAIDA2018/dataset/134300.dat",
                "/usr/share/dataset/CAIDA2018/dataset/134400.dat",
                "/usr/share/dataset/CAIDA2018/dataset/134500.dat",
                "/usr/share/dataset/CAIDA2018/dataset/134600.dat",
                "/usr/share/dataset/CAIDA2018/dataset/134700.dat",
                "/usr/share/dataset/CAIDA2018/dataset/134800.dat",
                "/usr/share/dataset/CAIDA2018/dataset/134900.dat",
                "/usr/share/dataset/CAIDA2018/dataset/135000.dat",
                "/usr/share/dataset/CAIDA2018/dataset/135100.dat",
                "/usr/share/dataset/CAIDA2018/dataset/135200.dat",
                "/usr/share/dataset/CAIDA2018/dataset/135300.dat",
                "/usr/share/dataset/CAIDA2018/dataset/135400.dat",
                "/usr/share/dataset/CAIDA2018/dataset/135500.dat",
                "/usr/share/dataset/CAIDA2018/dataset/135600.dat", 
                "/usr/share/dataset/CAIDA2018/dataset/135700.dat",
                "/usr/share/dataset/CAIDA2018/dataset/135800.dat",
                "/usr/share/dataset/CAIDA2018/dataset/135900.dat"};

ifstream fin;

const int M = 2e9;
const double unitTime = 0.01;

/*
string datapath[50] = {"../formatted00.dat",
                "../formatted01.dat",
                "../formatted02.dat",
                "../formatted03.dat",
                "../formatted04.dat"};*/

/*pair <uint64_t, double> Read()
{   
    static bool isfirstRead = true;
    static int curFinID = 0;
    static double offset = 0;
    static double lastT = 0;

	double t; uint64_t s;
    if (isfirstRead) {
        isfirstRead = false;
        fin.open(datapath[curFinID], std :: ios :: binary);
    }

	if (fin.eof()) {
        fin.close();
        fin.open(datapath[++curFinID], std :: ios :: binary);
    }
    fin.read((char*)&t, sizeof(double));
    fin.read((char*)&s, sizeof(uint64_t));

    
    t += offset;
    if(t < lastT) {
        offset += lastT - t;
        t += (lastT - t);
    }
    lastT = t;

    t *= 1e308;
	return make_pair(s, t); //假设时间戳的单位是s，现在的效果是每1s有50W个item流过
}*/
pair <uint64_t, double> Read()
{   
    static bool isfirstRead = true;
    static int curFinID = 0;
    static double offset = 0;
    static double lastT = 0;

	double t; uint64_t s, _s;
    if (isfirstRead) {
        isfirstRead = false;
        fin.open(datapath[curFinID], std :: ios :: binary);
    }

	if (fin.eof()) {
        fin.close();
        fin.open(datapath[++curFinID], std :: ios :: binary);
        if (curFinID > 60) {
            fin.close();
            exit(0);
        }
    }
    fin.read((char*)&s, sizeof(uint64_t)); //srcip(4)+dstip(4) 
    fin.read((char*)&_s, 5);//srcport destport protcol

    fin.read((char*)&t, sizeof(double));

    
    t += offset;
    if(t < lastT) {
        offset += lastT - t;
        t += (lastT - t);
    }
    lastT = t;

	return make_pair(s, t * 2); 
}

map <uint64_t, double> startTime_predict;
map <uint64_t, double> endTime_predict;
map <uint64_t, int> counter;

int n[2], precision[2];
double aae[2];

void check(double p, double x, int type) {
    if (fabs(p - x) < delta) {
        ++precision[type];
        aae[type] += fabs(p - x);
    }
}

int mem;


bool checkRegular(double*batchValue, int batchCnt)
{
    if(batchCnt<=4)
        return false;
#if 0
    sort(batchValue, batchValue+batchCnt);
    double middle = batchValue[batchCnt/2];
#else
    double middle = 0;
    for(int i=0;i<batchCnt;i++)
        middle+=batchValue[i];
    middle/=batchCnt;
#endif
    int reg=0;
    for(int i=0;i<batchCnt;i++)
    {
        if(fabs(middle-batchValue[i])<0.2*middle)
            reg++;
    }
    if((double)reg/batchCnt>0.8)
    {
        //for(int i=0;i<batchCnt;i++)
        //    cerr<<1000*batchValue[i]<<" ";
        //cerr<<endl;
        return true;
    }
    else return false;
}

int main() {//包含元素个数
    freopen("result.txt", "w", stdout);
    cerr << "150KB" << endl;
    //Alg* durationTimeDetector = new Alg(6400, 10, 0.6); //L * 24 / 1024 = 150KB
    //Alg* intervalDetector = new Alg(6400, 10, 0.6);  //L * 24 / 1024 = 150KB
    //bloomfliter* timeStamp = new bloomfliter(2, 64000); //t * M * 8 / 1024 = 100KB
    //bloomfliter* startTime = new bloomfliter(2, 64000); //t * M * 8 / 1024 = 100KB
    unordered_map<uint64_t, list<double>*> durationTime;      
    unordered_map<uint64_t, list<double>*> batchStartTime;      
    unordered_map<uint64_t, list<uint32_t>*> batchCnt;      
    unordered_map<uint64_t, double>timeStamp;
    unordered_map<uint64_t, double>startTime;
    unordered_map<uint64_t, uint32_t>currentBatchCnt;

    for (int i = 0; i < M + 1; ++i) {
        auto e = Read(); 

        uint64_t id = e.first;
        double curTime = e.second;
        double lastTime = 0;
        bool flag=false;
        if(timeStamp.find(id)==timeStamp.end())
            flag=true;
        else
            lastTime = timeStamp[id];

        if(flag)
        {
            durationTime[id] = new list<double>;
            batchStartTime[id] = new list<double>;
            batchCnt[id] = new list<uint32_t>;
            startTime[id]=curTime;
            currentBatchCnt[id]=1;
        }
        else if (curTime - lastTime > unitTime)
        {
            durationTime[id]->push_back(timeStamp[id]-startTime[id]);
            batchStartTime[id]->push_back(startTime[id]);
            batchCnt[id]->push_back( currentBatchCnt[id]);
            startTime[id]=curTime;
            currentBatchCnt[id]=1;
        }
        else
            currentBatchCnt[id]+=1;
        timeStamp[id]=curTime;
        if(i==100000000)
        {
            cerr<<i<<endl;
            break;
        }
    }
    
    //unordered_map<uint64_t, list<double>*> durationTime;      
    //unordered_map<uint64_t, list<double>*> batchStartTime;      
    //unordered_map<uint64_t, list<uint32_t>*> batchCnt;      
    uint32_t totalIdCnt=0;
    uint32_t targetIdCnt=0;
    uint32_t totalItemCnt=0;
    uint32_t targetItemCnt=0;
    for(auto it=durationTime.begin();it!=durationTime.end();it++)
    {
        uint64_t id = it->first;
        list<double>&dT = *(it->second);
        list<double>&bT = *(batchStartTime[id]);
        list<uint32_t>&bC = *(batchCnt[id]);
        if(dT.size()==0)
            continue;
        int itemCnt=0;
        for(auto jt=bC.begin();jt!=bC.end();jt++)
            itemCnt+=*jt;
        totalIdCnt++;
        totalItemCnt+=itemCnt;
        if(itemCnt<=20)
            continue;
        
        auto idt = dT.begin();auto ibt = bT.begin();auto ibc = bC.begin();
        
        int batchCnt=0;
        double flagTime=-1;
        double*batchValue = new double[bC.size()];
        while(idt!=dT.end())
        {
            if(*ibc>=2)
            {
#if 0
                if(flagTime==-1)
                    flagTime=*ibt;
                else
                {
                    batchValue[batchCnt]=*ibt-flagTime;
                    flagTime=*ibt;
                    batchCnt++;
                }
#else
                batchValue[batchCnt]=*ibc;
                batchCnt++;
#endif           
            }
            ++idt;++ibt;++ibc;
        }
        if(!checkRegular(batchValue, batchCnt))
            continue;
        targetIdCnt++;
        targetItemCnt+=itemCnt;
        delete [] batchValue;
    }
    fprintf(stderr,"id ratio:%f, item ratio:%f\n",
        (double)targetIdCnt/totalIdCnt, (double)targetItemCnt/totalItemCnt);
    fin.close();
	return 0;
}
