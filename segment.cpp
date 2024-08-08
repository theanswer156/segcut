#include "segment.h"
#include <random>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <cmath>
#include <time.h>
#include <algorithm>
using namespace std;
segment::segment(const int& count)
{
    segcount = count;
    initial();
    setSrcData();
    getDesData();
    setCoeficient();
    setCoeficientBPrime();
    setCoeficientG();

     computeArcLength1();
     computeArcLength2();
     getIsoTime1();
     getIsoTime2();
    computeTimeSeq();
    computeCrossTime();
//    NewtonIterator();
    // violenceSolution();
     computeIsoPoint();
    computeSelfCrossPoint();
    findSmallRectangle();
}
segment::segment()
{
    initial();
    setSrcData();
    getDesData();

    setCoeficient();
    setCoeficientBPrime();
    setCoeficientG();
    computeArcLength1();
    computeArcLength2();
    getIsoTime1();
    getIsoTime2();
    computeTimeSeq();
    computeCrossTime();
 // NewtonIterator();
 // violenceSolution();
    computeIsoPoint();
    computeSelfCrossPoint();
    findSmallRectangle();
}


segment::~segment()
{
}
//      初始化所有有关变量的维数
void segment::initial()
{
    srcdata.resize(2);
    desdata.resize(2);
    segdata.resize(2);
    segtime.resize(2);
    timingseq.resize(2);
//  crosstime.resize(2);
    crosspoint.resize(2);
    coefficient.resize(2);
    coefBPrime.resize(2);
    coefficientG.resize(2);
    selfcrosstime.resize(2);
//  selfcrosspoint.resize(2);
    minRectangle.resize(2);

}

void segment::setSrcData()
{
     srand(time(0));
     for(int i = 0;i<=1;++i){
         for (int j = 0; j < 4; ++j)
         {
             srcdata[i].emplace_back(rand()/4e1,rand()/4e1);
         }
     }
    
//   srcdata[0].emplace_back(Point(354.275,547.75));
//   srcdata[0].emplace_back(Point(631.375,27.2));
//   srcdata[0].emplace_back(Point(805.75,647.125));
//   srcdata[0].emplace_back(Point(112.525,252.725));

//   srcdata[1].emplace_back(Point(0.0, 200.0));
//   srcdata[1].emplace_back(Point(500.0, 500.0));
//   srcdata[1].emplace_back(Point(0.0, 500.0));
//   srcdata[1].emplace_back(Point(274.0, 400.0));
}

void segment::getDesData()
{
    size_t size = srcdata[0].size();
    for (double t = 0; t < 1.0000; t += precis)
    {
        Point resultPoint0 = {0, 0};
        Point resultPoint1 = {0, 0};
        vector<double> bernstein1(size, 1), bernstein2(size, 1);
        for (int i = 1; i <= size - 1; ++i)
        {
            bernstein1[i] *= bernstein1[i - 1] * t;
            bernstein2[size - i - 1] *= bernstein2[size - i] * (1 - t);
        }
        for (int j = 0; j < size; ++j)
        {
            resultPoint0 += srcdata[0][j] * pascaTri[j + 1] * bernstein1[j] * bernstein2[j];
            resultPoint1 += srcdata[1][j] * pascaTri[j + 1] * bernstein1[j] * bernstein2[j];
        }
        desdata[0].emplace_back(resultPoint0);
        desdata[1].emplace_back(resultPoint1);
    }
    
}

void segment::setCoeficient()
{
    if (coefficient.empty()) return;
    //     coefficient[0]-[3]
    //     贝塞尔曲线的四个系数   t in [0,1]
    coefficient[0].emplace_back(-srcdata[0][0] + 3 * srcdata[0][1] - 3 * srcdata[0][2] + srcdata[0][3]);
    coefficient[0].emplace_back(3 * srcdata[0][0] - 6 * srcdata[0][1] + 3 * srcdata[0][2]);
    coefficient[0].emplace_back(-3 * srcdata[0][0] + 3 * srcdata[0][1]);
    coefficient[0].emplace_back(1*srcdata[0][0]);

    coefficient[1].emplace_back(-srcdata[1][0] + 3 * srcdata[1][1] - 3 * srcdata[1][2] + srcdata[1][3]);
    coefficient[1].emplace_back(3 * srcdata[1][0] - 6 * srcdata[1][1] + 3 * srcdata[1][2]);
    coefficient[1].emplace_back(-3 * srcdata[1][0] + 3 * srcdata[1][1]);
    coefficient[1].emplace_back(1 * srcdata[1][0]);
}
void segment::setCoeficientBPrime()
{
    if(coefBPrime.empty()) return;
    coefBPrime[0].emplace_back(3 * coefficient[0][0]);
    coefBPrime[0].emplace_back(2 * coefficient[0][1]);
    coefBPrime[0].emplace_back(1 * coefficient[0][2]);

    coefBPrime[1].emplace_back(3 * coefficient[1][0]);
    coefBPrime[1].emplace_back(2 * coefficient[1][1]);
    coefBPrime[1].emplace_back(1 * coefficient[1][2]);
}
void segment::setCoeficientG()
{
    if(coefficientG.empty()) return;
    coefficientG[0].emplace_back(coefBPrime[0][0].x*coefBPrime[0][0].x+coefBPrime[0][0].y*coefBPrime[0][0].y);
    coefficientG[0].emplace_back(4*(coefBPrime[0][0].x*coefBPrime[0][1].x+coefBPrime[0][0].y*coefBPrime[0][1].y));
    coefficientG[0].emplace_back(2*(coefBPrime[0][0].x*coefBPrime[0][2].x+coefBPrime[0][0].y*coefBPrime[0][2].y)
                        +4*(coefBPrime[0][1].x*coefBPrime[0][1].x+coefBPrime[0][1].y*coefBPrime[0][1].y));
    coefficientG[0].emplace_back(4*(coefBPrime[0][1].x*coefBPrime[0][2].x+coefBPrime[0][1].y*coefBPrime[0][2].y));
    coefficientG[0].emplace_back(coefBPrime[0][2].x*coefBPrime[0][2].x+coefBPrime[0][2].y*coefBPrime[0][2].y);

    coefficientG[1].emplace_back(coefBPrime[1][0].x*coefBPrime[1][0].x+coefBPrime[1][0].y*coefBPrime[1][0].y);
    coefficientG[1].emplace_back(4*(coefBPrime[1][0].x*coefBPrime[1][1].x+coefBPrime[1][0].y*coefBPrime[1][1].y));
    coefficientG[1].emplace_back(2*(coefBPrime[1][0].x*coefBPrime[1][2].x+coefBPrime[1][0].y*coefBPrime[1][2].y)
                        +4*(coefBPrime[1][1].x*coefBPrime[1][1].x+coefBPrime[1][1].y*coefBPrime[1][1].y));
    coefficientG[1].emplace_back(4*(coefBPrime[1][1].x*coefBPrime[1][2].x+coefBPrime[1][1].y*coefBPrime[1][2].y));
    coefficientG[1].emplace_back(coefBPrime[1][2].x*coefBPrime[1][2].x+coefBPrime[1][2].y*coefBPrime[1][2].y);
}



void segment::computeArcLength1()
{
    int n = 10;
    int m = 10;
    bool exitLoops = false;

    double tol = 1e-1;

    vector<vector<double>> romberg(n+1,vector<double>(m+1,0));
    double h = 1;
    romberg[0][0] = h*(compute_f1(0)+compute_f1(h))/2;
    for(int i = 1;i<=n;++i){
        h/=2;
        double rest = 0;

        //!     计算零阶龙贝格积分
        for(int j = 1;j<=pow(2,i-1);++j){
            rest+=compute_f1((2*j-1)*h);
        }
        romberg[i][0] =0.5*romberg[i-1][0]+rest*h;
        if(rest<tol){
            arclength[0]=romberg[i][0];
            break;
        }
        //!     计算M阶龙贝格积分
        for(int j = 1;j<=i;++j){
            double error = pow(pow(4,j)-1,-1)*(romberg[i][j-1]-romberg[i-1][j-1]);
            romberg[j-1][i] = error;
            romberg[i][j] = romberg[i][j-1]+error;
            if(abs(error)<tol){
                exitLoops = true;
                arclength[0]=romberg[i][j];
                break;
            }
        }
        if(exitLoops) break;
    }
    // for(int i = 0;i<romberg.size();++i){
    //     if(romberg[i][0]<1e-5) break;
    //     std::copy(romberg[i].begin(),romberg[i].end(),ostream_iterator<double>(std::cout," "));
    // }

}

void segment::computeArcLength2()
{
    int n = 10;
    int m = 10;
    bool exitLoops = false;

    double tol = 1e-1;

    vector<vector<double>> romberg(n+1,vector<double>(m+1,0));
    double h = 1;
    romberg[0][0] = h*(compute_f2(0)+compute_f2(h))/2;
    for(int i = 1;i<=n;++i){
        h/=2;
        double rest = 0;

        //!     计算零阶龙贝格积分
        for(int j = 1;j<=pow(2,i-1);++j){
            rest+=compute_f2((2*j-1)*h);
        }
        romberg[i][0] =0.5*romberg[i-1][0]+rest*h;
        if(rest<tol){
            arclength[1]=romberg[i][0];
            break;
        }
        //!     计算M阶龙贝格积分
        for(int j = 1;j<=i;++j){
            double error = pow(pow(4,j)-1,-1)*(romberg[i][j-1]-romberg[i-1][j-1]);
            romberg[j-1][i] = error;
            romberg[i][j] = romberg[i][j-1]+error;
            if(abs(error)<tol){
                exitLoops = true;
                arclength[1]=romberg[i][j];
                break;
            }
        }
        if(exitLoops) break;
    }
    // for(int i = 0;i<romberg.size();++i){
    //     if(romberg[i][0]<1e-5) break;
    //     std::copy(romberg[i].begin(),romberg[i].end(),ostream_iterator<double>(std::cout," "));
    // }

}


void segment::getIsoTime1()
{
    double tol = 1e-1;
    segtime[0].emplace_back(0);
    for(int i = 1;i<segcount-1e-5;++i){
        double left = segtime[0].at(segtime[0].size()-1);
        double right = 1;
        double subArcLength = static_cast<double>(i *1.0/segcount)*arclength[0];
        double biCutPrircis = 1e-5;
        while(right-left > biCutPrircis){
            double mid = (left+right)/2;
            double F_mid = computeSubArcLength1(mid);
            //std::cout<<std::endl;
            //std::cout<<"Compute the "<<i<<"th segment"<<std::endl;
            //std::cout<<"The arclength[0] from 0 to "<<std::setprecision(6)<<mid<<" is "
            //        <<std::setprecision(10)<<F_mid<<std::endl;
            if(abs(F_mid-subArcLength)<tol){
                segtime[0].emplace_back(mid);
                break;
            }
            if(F_mid>subArcLength){
                right = mid;
            }else{
                left = mid;
            }
        }
//!     segtime[0].emplace_back((left+right)/2);
//!     默认不会到这一步
    }
    segtime[0].emplace_back(1);
}
void segment::getIsoTime2()
{
    double tol = 1e-1;
    segtime[1].emplace_back(0);
    for(int i = 1;i<segcount-1e-5;++i){
        double left = segtime[1].at(segtime[1].size()-1);
        double right = 1;
        double subArcLength = static_cast<double>(i *1.0/segcount)*arclength[1];
        double biCutPrircis = 1e-5;
        while(right-left > biCutPrircis){
            double mid = (left+right)/2;

            double F_mid = computeSubArcLength2(mid);
            //std::cout<<std::endl;
            //std::cout<<"Compute the "<<i<<"th segment"<<std::endl;
            //std::cout<<"The arclength[1] from 0 to "<<std::setprecision(6)<<mid<<" is "
            //        <<std::setprecision(10)<<F_mid<<std::endl;
            if(abs(F_mid-subArcLength)<tol){
                segtime[1].emplace_back(mid);
                break;
            }
            if(F_mid>subArcLength){
                right = mid;
            }else{
                left = mid;
            }
        }
    }
    segtime[1].emplace_back(1);
}

//  牛顿迭代法计算等分点
void segment::NewtonIterator()
{
    int MaxIterator = 1e2;
    double tol = 1e-4;
    double damping = 1;
    segtime[0].emplace_back(0);
    for(double i = 1;i<segcount;++i){
        double subarclength = (i/segcount)*arclength[0];
        double x0 = segtime[0].at(segtime[0].size()-1)+1.00/segcount;
        for(int j = 0;j<MaxIterator;++j){
            double f0 = computeSubArcLength1(x0)-subarclength;
            if(abs(f0)<epslion){
                //printf("abs(f0): %lf <epslion,stop this loop",f0);
                //printf("The %lf th segtime[0] is %lf",i,x0);
                segtime[0].emplace_back(x0);
                break;
            }
            double df = compute_gradf1(x0);
            if(abs(df)<delta){
                //printf("abs(df): %lf <epslion,stop this loop",df);
                //printf("The %lf th segtime[0] is %lf",i,x0);
                segtime[0].emplace_back(x0);
                break;
            }
            double x1=x0 - damping*(f0/df);
            if(abs(x1-x0)<tol){
                //printf("x0: %lf, x1: %lf .abs(x1-x0)<tol,stop this loop",x0,x1);
                //printf("The %lf th segtime[0] is %lf",i,x0);
                segtime[0].emplace_back(x0);
                break;
            }
            if(j == MaxIterator-1){
                //printf("reaching the MaxIterator , stop this loop");
                //printf("The %lf th segtime[0] is %lf",i,x0);
                segtime[0].emplace_back(x0);
                break;
            }
            x0 = x1;
        }

    }
    segtime[0].emplace_back(1);
}

//      计算两条曲线在区域[begin1,end1]✖[begin2,end2]里可能的交点
void segment::NewtonIterator(const double &begin1, const double &end1, const double &begin2, const double &end2)
{   

    std::random_device rd; // 用于获得随机种子
    std::mt19937 gen(rd()); // 生成随机数引擎

    std::uniform_real_distribution<double> dis(0.0, 1.0); // 定义双精度浮点数的均匀分布
    int maxiteratortime = 2e3;
    double eps = 1;
    double tol = 1e-4;
    //  设置初始值为两个时间区间的中间点  即可行域的中心
    double s0 = (begin1+end1)/2;
    double t0 = (begin2+end2)/3;
    for(int k = 0;k<maxiteratortime;++k){
        if((s0 >end1 || s0 < begin1) || (t0 >end2 || t0 < begin2)){
            s0 = dis(gen);
            t0 = dis(gen);
            continue;
        }
        Point point1 = computePoint1(s0);
        Point point2 = computePoint2(t0);

        double f1 = pow(point1.x-point2.x,2);
        double f2 = pow(point1.y-point2.y,2);
        double f = f1+f2;
        if(abs(f)<eps){
            //printf("F is close to zero , stop this iteration.");
            // map.emplace(make_pair(s0,f));
            crosstime.emplace_back(s0);
            filterCrossPoint(s0,t0);
            s0 = dis(gen);
            t0 = dis(gen);
            continue;
        }

        double gradb1x = (point1.x-point2.x)*(coefBPrime[0][0].x*s0*s0+coefBPrime[0][1].x*s0+coefBPrime[0][2].x);
        double gradb1y = (point1.y-point2.y)*(coefBPrime[0][0].y*s0*s0+coefBPrime[0][1].y*s0+coefBPrime[0][2].y);
        double gradb2x = -(point1.x-point2.x)*(coefBPrime[1][0].x*t0*t0+coefBPrime[1][1].x*t0+coefBPrime[1][2].x);
        double gradb2y = -(point1.y-point2.y)*(coefBPrime[1][0].y*t0*t0+coefBPrime[1][1].y*t0+coefBPrime[1][2].y);
        double del = gradb1x*gradb2y-gradb1y*gradb2x;
        if(abs(del)<eps){
            //printf("DF is close to zero , stop this iteration. We get Nothing");
            s0 = dis(gen);
            t0 = dis(gen);
            continue;
        }
        double deltas = (gradb2y*f1-gradb1y*f2)/del;
        double deltat = (-gradb2x*f1+gradb1x*f2)/del;

        double s1 = s0-Adapting(deltas);
        double t1 = t0-Adapting(deltat);
        if(abs(s1-s0)+abs(t0-t1)<tol){
            //printf("(s1,t1) is close to (s0,t0), stop this iteration");
             crosstime.emplace_back(s1);
//            map.emplace(make_pair(s0,t0));
             filterCrossPoint(s1,t1);
            continue;
        }
        s0 = s1;
        t0 = t1;
    }
    //printf("Reaching the maxiterator time , there is no cross-point betweeen two curve in area (%lf,%lf) times (%lf,%lf).",begin1,end1,begin2,end2);
}

double segment::Adapting(double num)
{
    if(abs(num)<0.05) return num;
    while(abs(num)>0.05){
        num/=10;
    }
    return num;

}

void segment::filterCrossPoint(double s,double t)
{
    if(map.empty()){
        map.emplace(make_pair(s,t));
        return;
    }
    double f = dist(computePoint1(s),computePoint2(t));
    bool compared = false;
    auto it = map.begin();
    while (it != map.end())
    {   
        double time1 = (*it).first;
        double time2 = (*it).second;
        double f_new = dist(computePoint1(time1),computePoint2(time2));
        if( (abs(time1-s)+abs(time2-t)) < 1e-2){
            compared = true;
            if(f_new > f){
                map.emplace(make_pair(s,t));
                map.erase(time1);
                break;
            }
        }
        it++;
    }
    if(!compared){
        map.emplace(make_pair(s,t));
    }
    return;
}

void segment::computeIsoPoint()
{
    for(auto &time:segtime[0]){
    segdata[0].emplace_back(computePoint1(time));
    }
    for(auto &time:segtime[1]){
    segdata[1].emplace_back(computePoint2(time));
    }
}
//      计算贝塞尔曲线自相交的交点情况
void segment::computeSelfCrossPoint()
{
    for(int i = 0;i<=1;++i){
        if((coefficient[i].at(0) == Point{0,0})|| (coefficient[i].at(1) == Point{0,0})
                                            ||(coefficient[i].at(2) == Point{0,0})) return;
        double judgement = coefficient[i].at(0).x*coefficient[i].at(1).y-coefficient[i].at(0).y*coefficient[i].at(1).x;
        if(abs(judgement)<1e-10) return;
        double k1 = (coefficient[i].at(2).x*coefficient[i].at(1).y-coefficient[i].at(2).y*coefficient[i].at(1).x)/judgement;
        double k2 = (coefficient[i].at(0).x*coefficient[i].at(2).y-coefficient[i].at(0).y*coefficient[i].at(2).x)/judgement;
        if(k1>0 || k2>0) return;
        double delta = -4*k1-3*k2*k2;
        //!     当delta判别式等于0的时候  time1==time2  我们把这种情况剔除
        if(delta<=0) return;
        double time1 = (-k2+sqrt(delta))/2;
        double time2 = (-k2-sqrt(delta))/2;

        if(time1>=0 && time1 <=1 && time2>=0 && time2<=1){
            selfcrosstime[i] = time1;
        }
        if(selfcrosstime[i]>=0){
            printf("The bezier curve self-cross at time %llf \n",selfcrosstime[i]);
        }else{
            printf("The bezier curve is not self-cross \n");
        }

    }
    
}

void segment::violenceSolution()
{   
    double precise = 1e-3;
    for(double t1 = 0;t1<1.0;){
        Point point1 = computePoint1(t1);
        t1+=precise;
        Point point2 = computePoint1(t1);
        for(double t2 = 0;t2<1.0;){
            Point point3 = computePoint2(t2);
            t2+=precise;
            Point point4 = computePoint2(t2);
            if(isRectangleIntersecting(point1,point2,point3,point4)){
                crosstime.emplace_back(t1);
            }
        }
    }
}

void segment::computeCrossTime()
{
    //! 计算两个贝塞尔曲线交点  还是只能采用包络框的方式 由computTimeSeq函数计算
    //! 得到的两条曲线单调性变化的位置最后由二分法计算   可以得到所有的交点
    //! 所以我们要做的就是能够快速找到有交的包络框
    //! 由于timeseq1、timeseq2的定义  我们可以知道在两个时间段内曲线要么是严格凹的
    //! 要么是严格凸的  

    //      先做退出去的判断   这里不能对crosstime做判断  因为我们是要做递归调用 
    //      采用类似于剪枝的方法去做的
    //      其实可以暴力的使用这种方法去计算交点  比如将0-1以0.001分为1000份
    //      进行一百万次计算就可以精确到0.001  很是不错
    //      如果通过先暴力  平均分为10份   总共分成100个区域  结合曲线的单调性来
    //      计算  将在区间内严格单调的
    //      不是  我直接根据计算而来的两个贝塞尔曲线单调性变化的数据timeseq将
    //      两个贝塞尔曲线直接分成这些曲线段  然后用牛顿迭代法去找不就是了  这样岂不是更快
    //      这样最多要执行25次牛顿迭代  但是牛顿迭代是二次收敛的   怕什么
    //      这种对于两个线段有多个交点的情况还是搜索不到   但是这种情况最多也只有两个交点了
    //      可以搜索两次  就好了  而且是对于有交点的情况下再搜索一次就好了
    //      由于其曲线在区域内是严格凹或凸的   所以如果有交点的话  尾部在迭代一次就好了
    // for(int i = 0;i<timingseq[0].size()-1;++i){
    //     Point point1 = computePoint1(timingseq[0].at(i));
    //     Point point2 = computePoint1(timingseq[0].at(i+1));
    //     for(int j = 0;j<timingseq[1].size()-1;++j){
    //         Point point3 = computePoint2(timingseq[1].at(j));
    //         Point point4 = computePoint2(timingseq[1].at(j+1));
    //         if(isRectangleIntersecting(point1,point2,point3,point4)){
    //             NewtonIterator(0,1,0,1);
    //         }
    //     }
    // }
    NewtonIterator(0,1,0,1);
}

void segment::computeTimeSeq()
{   //  计算 x 参数方程单调性变化的两个时间点  我们只要单调性变化的时间点  
    //  恒为常数或者一直递增或递减可以不存入  这些情况都没有单调性的变化
    //  但是这样一直用if 确实不太好
    for(int i = 0;i<=1;++i){
        timingseq[i].emplace_back(0);
        if(true){
        if(coefBPrime[0].at(0).x != 0){
            double delta = coefBPrime[i].at(1).x*coefBPrime[i].at(1).x-4*coefBPrime[i].at(0).x*coefBPrime[i].at(2).x;
            if(delta>=0){
                double time1 = (-coefBPrime[i].at(1).x+sqrt(delta))/(2*coefBPrime[i].at(0).x);
                double time2 = (-coefBPrime[i].at(1).x-sqrt(delta))/(2*coefBPrime[i].at(0).x);
                if(time1>0 && time1 <1){
                    timingseq[i].emplace_back(time1);
                }
                if(time2>0 && time2 <1){
                    timingseq[i].emplace_back(time2);
                }
            }
        }else if(coefBPrime[i].at(1).x != 0){
            double time1 = -coefBPrime[i].at(2).x/coefBPrime[i].at(1).x;
            if(time1>0 && time1 <1){
                timingseq[i].emplace_back(time1);
            }
        }
        }
        //  计算 y 参数方程单调性变化的两个时间点
        if(true){
            if(coefBPrime[i].at(0).y != 0){
                double delta = coefBPrime[i].at(1).y*coefBPrime[i].at(1).y-4*coefBPrime[i].at(0).y*coefBPrime[i].at(2).y;
                if(delta>=0){
                    double time1 = (-coefBPrime[i].at(1).y+sqrt(delta))/(2*coefBPrime[i].at(0).y);
                    double time2 = (-coefBPrime[i].at(1).y-sqrt(delta))/(2*coefBPrime[i].at(0).y);
                    if(time1>0 && time1 <1){
                        timingseq[i].emplace_back(time1);
                    }
                    if(time2>0 && time2 <1){
                        timingseq[i].emplace_back(time2);
                    }
                }
            }else if(coefBPrime[i].at(1).y != 0){
                double time1 = -coefBPrime[i].at(2).y/coefBPrime[i].at(1).y;
                if(time1>0 && time1 <1){
                    timingseq[i].emplace_back(time1);
                }
            }
        }
    //      对单调性变化的时间进行逆排序   从1到0的逆序排序
    timingseq[i].emplace_back(1);
    std::sort(timingseq[i].begin(),timingseq[i].end());
    }
    
}

bool segment::isRectangleIntersecting(const Point &p1, const Point &p2, const Point &p3, const Point &p4)
{
    if((p1 == p2) || (p3 == p4) ) return false;
    double x1_min = std::min(p1.x, p2.x);
    double y1_min = std::min(p1.y, p2.y);
    double x1_max = std::max(p1.x, p2.x);
    double y1_max = std::max(p1.y, p2.y);

    double x2_min = std::min(p3.x, p4.x);
    double y2_min = std::min(p3.y, p4.y);
    double x2_max = std::max(p3.x, p4.x);
    double y2_max = std::max(p3.y, p4.y);

    // 检查矩形是否在x轴上相交
    bool intersect_x = x1_min < x2_max && x2_min < x1_max;
    // 检查矩形是否在y轴上相交
    bool intersect_y = y1_min < y2_max && y2_min < y1_max;

    return intersect_x && intersect_y;
}
void segment::mergeBounds(const Rectangle& rect, double& minX, double& maxX, double& minY, double& maxY)
{
    minX = std::min(minX, min(rect.LeftUp.x,rect.RightDown.x));
    maxX = std::max(maxX, max(rect.RightDown.x,rect.LeftUp.x));
    minY = std::min(minY, min(rect.RightDown.y,rect.LeftUp.y));
    maxY = std::max(maxY, max(rect.LeftUp.y,rect.RightDown.y));
}
void segment::findSmallRectangle()
{
    if (!timingseq[0].empty())      //      计算第一条曲线的外包络框
    {
        Point p1 = computePoint1(timingseq[0][0]);
        Point p2 = computePoint1(timingseq[0][1]);
        double minX = p1.x, maxX = p2.x, minY = p1.y, maxY = p2.y;
        for (int i = 0; i < timingseq[0].size()-1; ++i)
        {
            Point p3 = computePoint1(timingseq[0][i]);
            Point p4 = computePoint1(timingseq[0][i+1]);
            Rectangle rect{ p3,p4 };
            mergeBounds(rect,minX,maxX,minY,maxY);
        }
        minRectangle[0].LeftUp.x = minX;
        minRectangle[0].LeftUp.y = maxY;
        minRectangle[0].RightDown.x = maxX;
        minRectangle[0].RightDown.y = minY;
    }
    if (!timingseq[1].empty())      //      计算第二条曲线的外包络框
    {
        Point p1 = computePoint2(timingseq[1][0]);
        Point p2 = computePoint2(timingseq[1][1]);
        double minX = p1.x, maxX = p2.x, minY = p1.y, maxY = p2.y;
        for (int i = 0; i < timingseq[1].size()-1; ++i)
        {
            Point p3 = computePoint2(timingseq[1][i]);
            Point p4 = computePoint2(timingseq[1][i + 1]);
            Rectangle rect{ p3,p4 };
            mergeBounds(rect, minX, maxX, minY, maxY);
        }
        minRectangle[1].LeftUp.x = minX;
        minRectangle[1].LeftUp.y = maxY;
        minRectangle[1].RightDown.x = maxX;
        minRectangle[1].RightDown.y = minY;
    }

}



double segment::compute_f1(const double &t)
{
    if(t<0) return 1e-9;
    double result = 0;
    double xvalue = coefBPrime[0].at(0).x*t*t+coefBPrime[0][1].x*t+coefBPrime[0][2].x;
    double yvalue = coefBPrime[0][0].y*t*t+coefBPrime[0][1].y*t+coefBPrime[0][2].y;
    result = pow(xvalue,2)+pow(yvalue,2);


    result = sqrt(abs(result));
//    printf("f(%lf)=%lf",t,result);
    return result;
}
double segment::compute_f2(const double &t)
{
    if(t<0) return 1e-9;
    double result = 0;
    double xvalue = coefBPrime[1][0].x*t*t+coefBPrime[1][1].x*t+coefBPrime[1][2].x;
    double yvalue = coefBPrime[1][0].y*t*t+coefBPrime[1][1].y*t+coefBPrime[1][2].y;
    result = pow(xvalue,2)+pow(yvalue,2);


    result = sqrt(abs(result));
//    printf("f(%lf)=%lf",t,result);
    return result;
}

double segment::compute_gradf1(const double &t)
{
    double xvalue = coefBPrime[0].at(0).x*t*t+2*coefBPrime[0].at(1).x*t+coefBPrime[0].at(2).x;
    double xsign = xvalue<0?-1:1;
    xvalue = abs(xvalue)*xsign*2*(coefBPrime[0].at(0).x*t+coefBPrime[0].at(1).x);
    double yvalue = coefBPrime[0].at(0).y*t*t+2*coefBPrime[0].at(1).y*t+coefBPrime[0].at(2).y;
    double ysign = yvalue<0?-1:1;
    yvalue = abs(yvalue)*ysign*2*(coefBPrime[0].at(0).y*t+coefBPrime[0].at(1).y);
    return xvalue+yvalue;
}
double segment::compute_gradf2(const double &t)
{
    double xvalue = coefBPrime[1].at(0).x*t*t+2*coefBPrime[1].at(1).x*t+coefBPrime[1].at(2).x;
    double xsign = xvalue<0?-1:1;
    xvalue = abs(xvalue)*xsign*2*(coefBPrime[1].at(0).x*t+coefBPrime[1].at(1).x);
    double yvalue = coefBPrime[1].at(0).y*t*t+2*coefBPrime[1].at(1).y*t+coefBPrime[1].at(2).y;
    double ysign = yvalue<0?-1:1;
    yvalue = abs(yvalue)*ysign*2*(coefBPrime[1].at(0).y*t+coefBPrime[1].at(1).y);
    return xvalue+yvalue;
}

double segment::computeSubArcLength1(const double &t)
{
    int n = 10;
    int m = 10;
    bool exitLoops = false;

    double tol = 1e-1;

    double result = 0;
    vector<vector<double>> romberg(n+1,vector<double>(m+1,0));
    double h = t;
    romberg[0][0] = h*(compute_f1(0)+compute_f1(h))/2;
    for(int i = 1;i<=n;++i){
        h/=2;
        double rest = 0;

        //!     计算零阶龙贝格积分
        for(int j = 1;j<=pow(2,i-1);++j){
            rest+=compute_f1((2*j-1)*h);
        }
        romberg[i][0] =0.5*romberg[i-1][0]+rest*h;
        if(rest<tol){
            result = romberg[i][0];
            break;
        }
        //!     计算M阶龙贝格积分
        for(int j = 1;j<=i;++j){
            double error = pow(pow(4,j)-1,-1)*(romberg[i][j-1]-romberg[i-1][j-1]);
            romberg[j-1][i] = error;
            romberg[i][j] = romberg[i][j-1]+error;
            if(abs(error)<tol){
                exitLoops = true;
                result = romberg[i][j];
                break;
            }
        }
        if(exitLoops) break;
    }
    //std::cout<<"\n";
    //printf("The first curve arclength from 0 to %lf is %lf",t,result);
//    for(int i = 0;i<romberg.size();++i){
//        if(romberg[i][0]<1e-5) break;
//        std::cout<<romberg[i];
//    }


    return result;
}

double segment::computeSubArcLength2(const double &t)
{

    int n = 10;
    int m = 10;
    bool exitLoops = false;

    double tol = 1e-1;

    double result = 0;
    vector<vector<double>> romberg(n+1,vector<double>(m+1,0));
    double h = t;
    romberg[0][0] = h*(compute_f2(0)+compute_f2(h))/2;
    for(int i = 1;i<=n;++i){
        h/=2;
        double rest = 0;

        //!     计算零阶龙贝格积分
        for(int j = 1;j<=pow(2,i-1);++j){
            rest+=compute_f2((2*j-1)*h);
        }
        romberg[i][0] =0.5*romberg[i-1][0]+rest*h;
        if(rest<tol){
            result = romberg[i][0];
            break;
        }
        //!     计算M阶龙贝格积分
        for(int j = 1;j<=i;++j){
            double error = pow(pow(4,j)-1,-1)*(romberg[i][j-1]-romberg[i-1][j-1]);
            romberg[j-1][i] = error;
            romberg[i][j] = romberg[i][j-1]+error;
            if(abs(error)<tol){
                exitLoops = true;
                result = romberg[i][j];
                break;
            }
        }
        if(exitLoops) break;
    }
    //std::cout<<"\n";
    //printf("The second curve arclength from 0 to %lf is %lf",t,result);
//    for(int i = 0;i<romberg.size();++i){
//        if(romberg[i][0]<1e-5) break;
//        std::cout<<romberg[i];
//    }


    return result;
}

double segment::computeSubArcLength1(const double &begin, const double &end)
{

    int n = 10;
    int m = 10;
    bool exitLoops = false;

    double tol = 1e-1;

    double result = 0;
    vector<vector<double>> romberg(n+1,vector<double>(m+1,0));
    double t1 = begin;
    double t2 = end;
    double h = t2-t1;
    romberg[0][0] = h*(compute_f1(t1)+compute_f1(t2))/2;
    for(int i = 1;i<=n;++i){
        h/=2;
        double rest = 0;

        //!     计算零阶龙贝格积分
        for(int j = 1;j<=pow(2,i-1);++j){
            rest+=compute_f1(t1+(2*j-1)*h);
        }
        romberg[i][0] =0.5*romberg[i-1][0]+rest*h;
        if(rest<tol){
            result = romberg[i][0];
            break;
        }
        //!     计算M阶龙贝格积分
        for(int j = 1;j<=i;++j){
            double error = pow(pow(4,j)-1,-1)*(romberg[i][j-1]-romberg[i-1][j-1]);
            romberg[j-1][i] = error;
            romberg[i][j] = romberg[i][j-1]+error;
            if(abs(error)<tol){
                exitLoops = true;
                result = romberg[i][j];
                break;
            }
        }
        if(exitLoops) break;
    }
    //std::cout<<"\n";
    //printf("The first arclength from %lf to %lf is %lf",t1,t2,result);
    // for(int i = 0;i<romberg.size();++i){
    //     if(romberg[i][0]<1e-5) break;
    //     std::copy(romberg[i].begin(),romberg[i].end(),std::ostream_iterator<double>(std::cout," "));

    // }
    return result;
}


double segment::computeSubArcLength2(const double &begin, const double &end)
{

    int n = 10;
    int m = 10;
    bool exitLoops = false;

    double tol = 1e-1;

    double result = 0;
    vector<vector<double>> romberg(n+1,vector<double>(m+1,0));
    double t1 = begin;
    double t2 = end;
    double h = t2-t1;
    romberg[0][0] = h*(compute_f2(t1)+compute_f2(t2))/2;
    for(int i = 1;i<=n;++i){
        h/=2;
        double rest = 0;

        //!     计算零阶龙贝格积分
        for(int j = 1;j<=pow(2,i-1);++j){
            rest+=compute_f2(t1+(2*j-1)*h);
        }
        romberg[i][0] =0.5*romberg[i-1][0]+rest*h;
        if(rest<tol){
            result = romberg[i][0];
            break;
        }
        //!     计算M阶龙贝格积分
        for(int j = 1;j<=i;++j){
            double error = pow(pow(4,j)-1,-1)*(romberg[i][j-1]-romberg[i-1][j-1]);
            romberg[j-1][i] = error;
            romberg[i][j] = romberg[i][j-1]+error;
            if(abs(error)<tol){
                exitLoops = true;
                result = romberg[i][j];
                break;
            }
        }
        if(exitLoops) break;
    }
    //std::cout<<"\n";
    //printf("The arclength from %lf to %lf is %lf",t1,t2,result);
    // for(int i = 0;i<romberg.size();++i){
    //     if(romberg[i][0]<1e-5) break;
    //     std::copy(romberg[i].begin(),romberg[i].end(),std::ostream_iterator<double>(std::cout," "));

    // }
    return result;
}

double segment::dist(const Point &p1, const Point &p2)
{
    return sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2));
}

Point segment::computePoint1(double &t)
{
    Point resultPoint;

    int size = srcdata[0].size();
    vector<double> bernstein1(size,1),bernstein2(size,1);
    for(int i = 1;i<=size-1;++i){
        bernstein1[i]*=bernstein1[i-1]*t;
        bernstein2[size-i-1]*=bernstein2[size-i]*(1-t);
    }
    for(int i = 0;i<size;++i){
        resultPoint+=srcdata[0][i]*pascaTri[i+1]*bernstein1[i]*bernstein2[i];
    }
    return resultPoint;
}

Point segment::computePoint2(double &t)
{
    Point resultPoint;

    int size = srcdata[1].size();
    vector<double> bernstein1(size,1),bernstein2(size,1);
    for(int i = 1;i<=size-1;++i){
        bernstein1[i]*=bernstein1[i-1]*t;
        bernstein2[size-i-1]*=bernstein2[size-i]*(1-t);
    }
    for(int i = 0;i<size;++i){
        resultPoint+=srcdata[1][i]*pascaTri[i+1]*bernstein1[i]*bernstein2[i];
    }

    return resultPoint;
}

vector<vector<Point>> segment::outSegData()
{
    return segdata;
}

vector<vector<Point>> segment::outSrcData()
{
    return srcdata;
}

vector<vector<Point>> segment::outDesData()
{
    return desdata;
}

vector<Point> segment::outSelfCrossPoint()
{   
    if(selfcrosstime.empty()) return vector<Point>();
    for(int i = 0;i<selfcrosstime.size();++i){
        if(i%2 == 0){
            selfcrosspoint.emplace_back(computePoint1(selfcrosstime[i]));
        }else{
            selfcrosspoint.emplace_back(computePoint2(selfcrosstime[i]));
        }
    }
    return selfcrosspoint;
}

vector<Point> segment::outCrossPoint()
{
    if(map.empty()) return vector<Point>();
    for(auto pair:map){
        double time = pair.first;
        Point point = computePoint1(time);
        crosspoint.emplace_back(point);
    }
    return crosspoint;
}

vector<Rectangle> segment::outMinRectangle()
{
    return minRectangle;
}
