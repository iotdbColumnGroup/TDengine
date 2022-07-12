#ifndef _TD_TSDB_TSDB_GROUPAPPRENGINE_H
#define _TD_TSDB_TSDB_GROUPAPPRENGINE_H

# include "stdlib.h"

#define FLOAT_NULL  0x7FF00000
#define PI 3.1415926

#define true	1
#define false	0

typedef struct {
    long epsilon;
    int n;
    long start;
    long sigma;
}Feature;

typedef struct {
    int *columns; // list of cid
    int columnLen; // number of current columns
    int maxGain;
    int maxGainPosition;
    int timeSeriesLength;
    long interval;
    float *values;
    int valuesLen;
    int *gains;
    Feature **features;
    long *rangeMap;
    int rangeMapLen;
    int **rangeMapCols;
    int *rangeMapColsLen;
    int groupNum;
}ApprColumnGroup;


void groupingAppr(long *timestamps, float **dataColumns, int numOfRows, int numOfCols, int *columnSchemaIdx);

#endif //_TD_TSDB_TSDB_GROUPAPPRENGINE_H



