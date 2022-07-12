#ifndef _TD_TSDB_TSDB_GROUPENGINE_H
#define _TD_TSDB_TSDB_GROUPENGINE_H

# include "stdlib.h"

#define FLOAT_NULL  0x7FF00000
#define true	1
#define false	0

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
}ColumnGroup;

void grouping(long *timestamps, float **dataColumns, int numOfRows, int numOfCols, int *columnSchemaIdx);

#endif //_TD_TSDB_TSDB_GROUPENGINE_H



