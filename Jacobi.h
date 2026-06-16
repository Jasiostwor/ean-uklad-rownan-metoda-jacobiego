#ifndef JACOBI_H
#define JACOBI_H

#include <vector>
#include <string>
#include "interval.h"

using namespace std;
using namespace interval_arithmetic;

typedef vector<vector<long double>> MatrixNormal;
typedef vector<long double> VectorNormal;

typedef vector<vector<Interval<long double>>> MatrixInterval;
typedef vector<Interval<long double>> VectorInterval;

struct JacobiResult {
    bool success;
    string message;
    string iterationLog;
    int st;
    int it;
    vector<string> finalX;
};

// Funkcja pomocnicza:
string DoubleToString(long double val);
string IntervalToString(Interval<long double> val);

// Główne funkcje obliczeniowe
JacobiResult SolveJacobiNormal(const MatrixNormal& A, const VectorNormal& B, const VectorNormal& X0, int maxIter, long double eps);
JacobiResult SolveJacobiInterval(const MatrixInterval& A, const VectorInterval& B, const VectorInterval& X0, int maxIter, long double eps);

#endif
