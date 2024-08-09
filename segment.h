#ifndef SEGMENT_H
#define SEGMENT_H
#include <vector>
#include <math.h>
#include <map>
using namespace std;
struct TimePair
{
    double begin;
    double end;
    TimePair(double begin = 0,double end = 0):begin(begin),end(end){}

    friend bool operator< (const TimePair& pair1,const TimePair& pair2){
        return (pair1.begin > pair2.begin) && (pair1.end < pair2.end);
    }  
    friend TimePair operator-(const TimePair& pair1,const TimePair& pair2){
        return {pair1.begin-pair2.begin,pair1.end-pair2.end};
    }
    friend TimePair operator+(const TimePair& pair1,const TimePair& pair2){
        return {pair1.begin+pair2.begin,pair1.end+pair2.end};
    }   
     
    
};



struct Point{
    double x;
    double y;
    //      在这里可以直接定义点的模长以及叉乘了
    Point(double x = 0.0, double y = 0.0) : x(x), y(y) {}
    Point operator-()const{
        return Point(-x,-y);
    }
    double dist(const Point& p1,const Point& p2){
        return sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2));
    }
    friend Point operator*(const Point& p,double scalar){
        return {p.x*scalar,p.y*scalar};
    }
    friend Point operator*(double scalar,const Point& p){
        return p*scalar;
    }
    friend Point operator*(const Point& p1,const Point& p2){
        return {p1.x*p2.x,p1.y*p2.y};
    }
    friend Point operator+(double scalar,const Point& p){
        return {p.x+scalar,p.y+scalar};
    }
    friend Point operator+(const Point& p,double scalar){
        return scalar+p;
    }
    friend Point operator+(const Point& p1,const Point& p2){
        return {p1.x+p2.x,p1.y+p2.y};
    }
    friend Point operator-(const Point& p1,const Point& p2){
        return {p1.x-p2.x,p1.y-p2.y};
    }
    bool operator==(const Point& other) const {
        return (x == other.x) && (y == other.y);
    }
    Point& operator+=(const Point& p1){
        this->x += p1.x;
        this->y += p1.y;
        return *this;
    }
};

struct Rectangle
{
    Point LeftUp;
    Point RightDown;
};

struct  TreeNode
{
    double Begin;
    double End;
    TreeNode(double begin =  0, double end = 0) :Begin(begin), End(end){}
    ~TreeNode() {}
};

class segment{
public:
    segment(const int& count);
    segment();
    ~segment();
    void initial();                         //      初始化所有的私有成员变量是二维的
    void setSrcData();                      //      随机设置贝塞尔曲线的控制点
    void getDesData();                      //      利用矩阵方法计算贝塞尔曲线
    void setCoeficient();                   //      设置贝塞尔曲线关于时间t的参数方程的参数
    void setCoeficientG();                  //      设置函数 g 的参数
    void setCoeficientBPrime();

    void computeArcLength1();                //      计算贝塞尔曲线的弧长
    void computeArcLength2();                //      计算贝塞尔曲线的弧长

    void getIsoTime1();                      //      计算贝塞尔曲线等距点的时间值
    void getIsoTime2();                      //      计算贝塞尔曲线等距点的时间值
    void NewtonIterator();                   //      牛顿迭代方法计算等距点的时间值
    void NewtonIterator(const double& begin1,const double& end1,const double& begin2,const double& end2);
                                             //      计算两条曲线在区域[begin1,end1]✖[begin2,end2]里可能的交点
    void TreeTraverse(TreeNode* leftNode,TreeNode* rightNode);
    double Adapting(double num);
    void filterCrossPoint(double s,double t);
    void filterCrossPoint();                 //     过滤所有的点
    // double armijoLineSearch(const double s,const double t,const double& alpha0);               //      线性搜索最佳步长
    void computeIsoPoint();
    void computeSelfCrossPoint();             //      计算贝塞尔曲线自相交的交点
    void violenceSolution();                  //      暴力方法找交点
    void computeCrossTime();                 //      计算两个贝塞尔曲线交点
    void computeTimeSeq();                  
    bool isRectangleIntersecting(const Point& p1, const Point& p2, const Point& p3, const Point& p4);// 
    void mergeBounds(const Rectangle& rect, double& minX, double& maxX, double& minY, double& maxY);
    void findSmallRectangle();               //      寻找最小包络框
    double compute_f1(const double &t);      //      计算 f1 对应的函数值
    double compute_f2(const double &t);      //      计算 f2 对应的函数值
    double compute_gradf1(const double &t);  //      计算 f1 的导数
    double compute_gradf2(const double &t);  //      计算 f2 的导数
    double computeSubArcLength1(const double &begin);           //计算从0到t时间的曲线弧长
    double computeSubArcLength1(const double &begin,const double &end);//计算从begin到end时间的曲线弧长
    double computeSubArcLength2(const double &begin);           //计算从0到t时间的曲线弧长
    double computeSubArcLength2(const double &begin,const double &end);//计算从begin到end时间的曲线弧长
    double dist(const Point& p1,const Point& p2);
    Point computePoint1(double &t);
    Point computePoint2(double &t);         //    C++名字查找先于类型检查
    vector<vector<Point>> outSegData();     //    提供公共接口访问私有数据
    vector<vector<Point>> outSrcData();     //    提供公共接口访问私有数据
    vector<vector<Point>> outDesData();     //    提供公共接口访问私有数据
    vector<Point> outSelfCrossPoint();      //    提供公共接口访问曲线自相交的点
    vector<Point> outCrossPoint();          //    提供公共接口访问曲线相交的点
    vector<Rectangle> outMinRectangle();    //    返回曲线的外包络框
private:
    vector<vector<Point>> srcdata;      //  两条曲线的四个控制点

    vector<vector<Point>> desdata;
    vector<vector<Point>> segdata;      //  两条曲线的等距分割点

    vector<vector<double>> segtime;     //  两条曲线的等距分割的时间点

    vector<vector<double>> timingseq;   //  根据coefficientBPrime计算四个参数方程单调性改变的时间点

    vector<double> selfcrosstime;
    vector<Point> selfcrosspoint;       //  分别记录两条曲线自相交的点


    vector<double> crosstime;           //  记录两条贝塞尔曲线相交点的时间值                                
    map<double,double> map;             //  保存相交的时间点以及对应的距离值
    vector<Point> crosspoint;           //  记录两条曲线相交点

    vector<vector<Point>> coefficient;
    vector<vector<Point>> coefBPrime;
    vector<vector<double>> coefficientG;
    vector<int> pascaTri{0,1,3,3,1};
        
    vector<Rectangle> minRectangle;
    bool corsslabel = false;

    vector<double> arclength{0.0,0.0};
    double precis = 1e-2;
    double delta = 1e-5;
    double epslion = 1e-5;
    double tolerenceerror = 1e-7;

    int segcount = 50;



};
#endif
