#include "SCI.h"
#include <algorithm>



double myLog2(double v) {
    if (v == 0.0) {
        return 0.0;
    }
    else {
        return log2(v);
    }
}
double log2fac(int n) {
    double sum = 0;
    for (int i = 2; i <= n; i++) {
        sum += myLog2(i);
    }
    return sum;
}

double stirling(int n, int k) {
    return (double)k * myLog2(n) - log2fac(k);
}

double log2nChoosek(int n, int k) {
    if (k > n || k == 0) {
        return 0;
    }
    else {
        return log2fac(n) - log2fac(k) - log2fac(n - k);
    }
}
double binaryRegretPrecal(int M) {
    int p = 10;
    if (M < 1)
        return 0.0;
    double sum = 1.0;
    double b = 1.0;
    int bound = (int)ceil(2.0 + sqrt(2.0 * M * p * log((double)10.0)));
    for (int i = 1; i <= bound; i++) {
        b = (M - i + 1) * (b / M);
        sum += b;
    }
    return sum;
}
double regretPrecal(int M, int K) {
    if (K < 1) {
        return 0.0;
    }
    else if (K == 1) {
        return 1.0;
    }
    else {
        double sum = binaryRegretPrecal(M);
        double old_sum = 1.0;
        if (K > 2) {
            for (int j = 3; j <= K; j++) {
                double new_sum = sum + (M * old_sum) / ((double)j - 2.0);
                old_sum = sum;
                sum = new_sum;
            }
        }
        return sum;
    }
}
double regret(int M, int K) {
    if (K > 100) {
        double alpha = (double)K / (double)M;
        double ca = 0.5 + 0.5 * sqrt(1.0 + 4.0 / alpha);
        double logReg = (double)M * (log(alpha) + (alpha + 2.0) * log(ca) - 1.0 / ca) - 0.5 * log(ca + 2.0 / alpha);
        return logReg / log((double)2.0);
    }
    else {
        double costs = regretPrecal(M, K);
        if (costs <= 0.0)
            return 0.0;
        else
            return myLog2(costs);
    }
}
double entropy(std::map<int, int>& counts, int sum) {
    double res = 0.0;
    for (std::map<int, int>::iterator it = counts.begin(); it != counts.end(); ++it) {
        int c = it->second;
        double frac = ((double)c) / ((double)sum);
        res -= frac * myLog2(frac);
    }
    return res;
}
double entropy(std::vector<int>& counts) {
    int sum = counts[counts.size() - 1];
    if (sum == 0) {
        return 0.0;
    }
    else {
        double res = 0.0;
        for (size_t i = 0; i < (counts.size() - 1); i++) {
            int c = counts[i];
            double frac = ((double)c) / ((double)sum);
            res -= frac * myLog2(frac);
        }
        return res;
    }
}

double SC(std::vector<int>& x) {
    int n = x.size();
    std::map<int, int> counts;
    for (int i = 0; i < n; i++) {
        counts[x[i]]++;
    }
    double hx = (double)n * entropy(counts, n);
    double rx = regret(n, counts.size());
    return hx + rx;
}

double conditionalSC(std::vector<int>& x, std::vector<int>& y) {
    int n = x.size();
    std::map<int, std::map<int, int> > condCounts;
    std::map<int, int> counts;
    std::map<int, int> totals;
    for (int i = 0; i < n; i++) {
        counts[x[i]]++;
        (condCounts[y[i]])[x[i]]++;
        totals[y[i]]++;
    }
    int domX = counts.size();
    // compute sc(x|y)
    double scxgy = 0.0;
    for (std::map<int, std::map<int, int> >::iterator it = condCounts.begin(); it != condCounts.end(); ++it) {
        int key = it->first;
        std::map<int, int> currCounts = it->second;
        int currTotal = totals[key];
        scxgy += (double)currTotal * entropy(currCounts, currTotal) + regret(currTotal, domX);
    }
    return scxgy;
}
double conditionalEntropy(const std::vector<int>& x, const std::vector<int>& y) {
    int n = x.size();
    std::map<int, std::map<int, int> > condCounts;
    std::map<int, int> counts;
    std::map<int, int> totals;
    for (int i = 0; i < n; i++) {
        counts[x[i]]++;
        (condCounts[y[i]])[x[i]]++;
        totals[y[i]]++;
    }
    // compute H(x|y)
    double hxgy = 0.0;
    for (std::map<int, std::map<int, int> >::iterator it = condCounts.begin(); it != condCounts.end(); ++it) {
        int key = it->first;
        std::map<int, int> currCounts = it->second;
        int currTotal = totals[key];
        hxgy += ((double)currTotal / (double)n) * entropy(currCounts, currTotal);
    }
    return hxgy;
}

std::vector<int> getNiceCategories(std::vector<int>& v) {
    std::map<int, int> converter;
    std::vector<int> niceV;
    int currentCategory = 0;
    for (size_t i = 0; i < v.size(); i++) {
        int current = v[i];
        std::map<int, int>::iterator it = converter.find(current);
        if (it != converter.end()) {
            niceV.push_back(it->second);
        }
        else {
            niceV.push_back(currentCategory);
            converter[current] = currentCategory;
            currentCategory++;
        }
    }
    niceV.push_back(currentCategory);
    return niceV;
}
using std::vector;

std::vector<int> matrixToVector(std::vector<std::vector<int>>& yy) {
    int n = yy.size();          //行数
    int m = yy[0].size();       //列数

    // create y as single vector
    std::vector<int> y(n, 0);
    int maxDom = 1;
    for (int j = 0; j < m; j++) {
        std::set<int> setj;
        for (int i = 0; i < n; i++) {
            int curr = yy[i][j];
            setj.insert(curr);
            if (curr > maxDom) {
                maxDom = curr;
            }
        }
    }
    maxDom++;
    int currentF = 1;
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            int curr = (yy[i][j] + 1) * currentF;
            y[i] += curr;
        }
        currentF *= (maxDom + 1);
    }
    std::vector<int> yN = getNiceCategories(y);
    return yN;
}

std::vector<int> joinVectors(std::vector<int> yy, std::vector<int>& zz) {
    int n = yy.size();

    // create y as single vector
    int maxDom = 1;
    std::set<int> setj;
    for (int i = 0; i < n; i++) {
        int curr = yy[i];
        setj.insert(curr);
        if (curr > maxDom) {
            maxDom = curr;
        }
    }
    for (int i = 0; i < n; i++) {
        int curr = zz[i];
        setj.insert(curr);
        if (curr > maxDom) {
            maxDom = curr;
        }
    }
    maxDom++;
    int currentF = maxDom + 1;
    for (int i = 0; i < n; i++) {
        int curr = (zz[i] + 1) * currentF;
        yy[i] += curr;
    }
    std::vector<int> yN = getNiceCategories(yy);
    return yN;
}




double indepAsymNML(std::vector<vector<int>>& xEXP, std::vector<vector<int>>& yEXP, std::vector<std::vector<int>>& zEXP) {
    // convert input

    vector<vector<int>> xx(xEXP);
    vector<vector<int>> yy(yEXP);
    vector<vector<int>> zz(zEXP);

    //int n = yy.nrow();

    std::vector<int> xN = matrixToVector(xx);
    std::vector<int> yN = matrixToVector(yy);
    std::vector<int> zN = matrixToVector(zz);
    xN.pop_back();
    yN.pop_back();
    zN.pop_back();

    // join z*, y
    std::vector<int> yzN = joinVectors(yN, zN);
    yzN.pop_back();

    // additional model costs
    double result = 0.0;
    //double modelCosts = 0.0;
    // calculate costs
    result = conditionalSC(xN, zN) - conditionalSC(xN, yzN);

    //result += modelCosts;

     // return result
    return result;
}

double indepAsymNML(std::vector<vector<int>>& xEXP, std::vector<vector<int>>& yEXP) {
    // convert input

    vector<vector<int>> xx(xEXP);
    vector<vector<int>> yy(yEXP);

    std::vector<int> xN = matrixToVector(xx);
    std::vector<int> yN = matrixToVector(yy);
    xN.pop_back();
    yN.pop_back();

    // additional model costs
    double result = 0.0;
    //double modelCosts = 0.0;
    // calculate costs
    result = SC(xN) - conditionalSC(xN, yN);

    //result += modelCosts;

     // return result
    return result;
}



//与matalb的接口函数
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    double* x;
    double* y;
    double* z;
    double* CI;
    double cxgy, cygx;



    x = mxGetPr(prhs[0]); //得到输入矩阵的第一个元素的指针
    y = mxGetPr(prhs[1]);
    z = mxGetPr(prhs[2]);

    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    CI = mxGetPr(plhs[0]);


    int mrows = mxGetM(prhs[0]);   //获得矩阵的行
    int ncols = mxGetN(prhs[2]);   //获得矩阵的列

    vector<vector<int>> X(mrows, vector<int>(1));
    vector<vector<int>> Y(mrows, vector<int>(1));
    vector<vector<int>> Z;

    //将mxArray转化为vector
    for (int i = 0; i < mrows; i++) {
       X[i][0] = x[i]; //拷贝矩阵的元素到vector of vector
       Y[i][0] = y[i];
        
    }

    if (ncols != 0) {
        
        Z.resize(mrows);  //初始化
        for (int i = 0; i < mrows; i++) {
            Z[i].resize(ncols);
        }

        for (int i = 0; i < mrows; i++) {
            for (int j = 0; j < ncols; j++) {
                Z[i][j] = z[j * mrows + i]; //拷贝矩阵的元素到vector of vector
            }
        }

        cxgy = indepAsymNML(X, Y, Z);
        cygx = indepAsymNML(Y, X, Z);
    }
    else {
        cxgy = indepAsymNML(X, Y);
        cygx = indepAsymNML(Y, X);
    }

    *CI = (cxgy > cygx) ? cxgy : cygx;
    
}

