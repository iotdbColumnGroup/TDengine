#include "tsdbGroupApprEngine.h"
#include "stdio.h"
#include "math.h"
#include "time.h"

#define APPR_GROUP_FILE "/workspaces/TDengine/results/group_result_appr.csv"
#define APPR_TIME_FILE "/workspaces/TDengine/results/time_exact_appr.csv"

static void outputTimeToFile(const char *filename, double t);
static void outputGroupingToFile(ApprColumnGroup **columnGroups, const char *filename, int groupNum, int *columnSchemaIdx);
static void CreateApprColumnGroup(ApprColumnGroup *columnGroup, int cid, int groupNum, float *values, int valuesLen, long *timestamps);
static void updateTimestampInfo(ApprColumnGroup *columnGroup, long *timestamps);
static void destroyApprColumnGroup(ApprColumnGroup *columnGroup);
static void outputResults(ApprColumnGroup **columnGroups, int groupNum, int *columnSchemaIdx);
static int computeGain(ApprColumnGroup *g1, ApprColumnGroup *g2);
static int estimateOverlap(ApprColumnGroup *g1, ApprColumnGroup *g2);
static float computeProbability(Feature *f1, Feature *f2);
static long lcm(long number1, long number2);

static int isFloatNull(float *val) {
    return *(__uint32_t *)val == FLOAT_NULL;
}

static void createFeature(Feature *feature, long epsilon, int n, long start, long sigma) {
    feature->epsilon = epsilon;
    feature->n = n;
    feature->start = start;
    feature->sigma = sigma;
}

static void CreateApprColumnGroup(ApprColumnGroup *columnGroup, int cid, int groupNum, float *values, int valuesLen, long *timestamps) {
    columnGroup->columns = (int *)malloc(sizeof(int) * groupNum);
    columnGroup->columns[0] = cid;
    columnGroup->columnLen = 1;
    columnGroup->maxGain = -1;
    columnGroup->maxGainPosition = -1;
    columnGroup->gains = (int *)malloc(sizeof(int) * groupNum);
    columnGroup->valuesLen = valuesLen;
    columnGroup->values = (float *)malloc(sizeof(float) * valuesLen);
    for (int i = 0; i < valuesLen; i++) {
        columnGroup->values[i] = values[i];
    }
    columnGroup->features = (Feature **)malloc(sizeof(Feature *) * valuesLen);
    columnGroup->rangeMap = (long *)malloc(sizeof(long) * valuesLen);
    columnGroup->rangeMapCols = (int **)malloc(sizeof(int) * valuesLen);
    columnGroup->rangeMapColsLen = (int *)malloc(sizeof(int) *valuesLen);
    columnGroup->rangeMapLen = 2;
    columnGroup->groupNum = groupNum;
    updateTimestampInfo(columnGroup, timestamps);
}

static void updateTimestampInfo(ApprColumnGroup *columnGroup, long *timestamps) {
    // count number of not-nan values
    columnGroup->timeSeriesLength = 0;
    for (int i = 0; i < columnGroup->valuesLen; i++) {
        if (!isFloatNull(&columnGroup->values[i])) {
            columnGroup->timeSeriesLength++;
        }
    }

    // also update features
    int intervalGran = 1;
    long intervalSum = 0;
    int intervalCount = 0;
    int intervalMax = 1000;
    int timestampNumber = columnGroup->valuesLen;

    for (int ts = 1; ts < timestampNumber; ts++) {
        long interval_ = timestamps[ts] - timestamps[ts-1];
        if (interval_ < intervalMax) {
            intervalSum += interval_;
            intervalCount++;
        }
    }
    long interval = intervalSum / intervalCount;
    columnGroup->interval = interval / intervalGran * intervalGran;

    long start = timestamps[0];
    long sigma = 0;
    long sigmaSum = 0;

    long offset = 0;
    for (int ts = 0; ts < timestampNumber; ts++) {
        sigma = timestamps[ts] - start - ts * interval - offset;
        sigma = sigma > 0 ? sigma : -sigma; // abslong(sigma)
        if (sigma > 10 * interval) {
            long tmp = timestamps[ts] - start - ts * interval;
            tmp = tmp > 0 ? tmp : -tmp;
            sigma = tmp % interval;
            offset = tmp / interval * interval;
        }
        sigmaSum += sigma;
    }
    sigma = sigmaSum / timestampNumber;
    columnGroup->features[0] = (Feature *)malloc(sizeof(Feature));
    createFeature(columnGroup->features[0], interval, columnGroup->timeSeriesLength, start, sigma);
    columnGroup->rangeMap[0] = start;
    columnGroup->rangeMap[1] = start + interval * columnGroup->timeSeriesLength;
    columnGroup->rangeMapCols[0] = (int *)malloc(columnGroup->groupNum * sizeof(int));
    columnGroup->rangeMapCols[0][0] = columnGroup->columns[0];
    columnGroup->rangeMapColsLen[0] = 1;
}

static void destroyApprColumnGroup(ApprColumnGroup *columnGroup) {
    free(columnGroup->gains);
    free(columnGroup->columns);
    free(columnGroup->values);
//    for (int i = 0; i < columnGroup->columnLen; i++) {
//        free(columnGroup->features[i]);
//    }
    for (int i = 0; i < columnGroup->rangeMapLen - 1; i++) {
        free(columnGroup->rangeMapCols[i]);
    }
    free(columnGroup->features);
    free(columnGroup->rangeMap);
    free(columnGroup->rangeMapCols);
    free(columnGroup->rangeMapColsLen);
    free(columnGroup);
}

static void mergeApprColumnGroup(ApprColumnGroup *g1, ApprColumnGroup *g2) {
    /** merge g2 values into g1 and update **/
    // merge columns, features
    for (int i = 0; i < g2->columnLen; i++) {
        g1->columns[i + g1->columnLen] = g2->columns[i];
        g1->features[i + g1->columnLen] = g2->features[i];
    }
    g1->columnLen = g1->columnLen + g2->columnLen;

    // update length
    int overlap = estimateOverlap(g1, g2);
    int timeSeriesLength = g1->timeSeriesLength + g2->timeSeriesLength - overlap;
    
    // merge range Map
    long *map1;
    long *map2;
    int **mapCols1;
    int **mapCols2;
    int *mapCols1Len;
    int *mapCols2Len;

    int map1Len;
    int map2Len;
    long *newMap = (long *)malloc(sizeof(long) * g1->valuesLen);
    int **newMapCols = (int **)malloc(sizeof(int) * g1->valuesLen);
    int *newMapColsLen = (int *)malloc(sizeof(int) * g1->valuesLen);

    int newMapPtr = 0;


    if (g1->rangeMapLen < g2->rangeMapLen) {
        map1 = g1->rangeMap;
        map2 = g2->rangeMap;
        mapCols1 = g1->rangeMapCols;
        mapCols2 = g2->rangeMapCols;
        map1Len = g1->rangeMapLen;
        map2Len = g2->rangeMapLen;
        mapCols1Len = g1->rangeMapColsLen;
        mapCols2Len = g2->rangeMapColsLen;
    } else {
        map1 = g2->rangeMap;
        map2 = g1->rangeMap;
        mapCols1 = g2->rangeMapCols;
        mapCols2 = g1->rangeMapCols;
        map1Len = g2->rangeMapLen;
        map2Len = g1->rangeMapLen;
        mapCols1Len = g2->rangeMapColsLen;
        mapCols2Len = g1->rangeMapColsLen;
    }

    int i = 0;
    int j = 0;
    int lastPoint = map1[0] == map2[0] ? 3 : 1; //1: map1, 2: map2, 3:equal

    newMap[newMapPtr] = map1[i];
    newMapPtr++;

    while (1) {
        if (lastPoint == 1) {
            i++;
            if (i >= map1Len) {
                break;
            }
            if (map1[i] < map2[j]) {
                newMap[newMapPtr] = map1[i];
                newMapCols[newMapPtr-1] = (int *)malloc(g1->groupNum * sizeof(int));
                for (int k = 0; k < mapCols1Len[i-1]; k++) {
                    newMapCols[newMapPtr-1][k] = mapCols1[i-1][k];
                }
                newMapColsLen[newMapPtr-1] = mapCols1Len[i-1];
                newMapPtr ++;
            } else{
                if (map1[i] == map2[j])
                    lastPoint = 3;
                else{
                    lastPoint = 2;
                }
                newMap[newMapPtr] = map2[j];
                if (j == 0) {
                    newMapCols[newMapPtr-1] = (int *)malloc(g1->groupNum * sizeof(int));
                    for (int k = 0; k < mapCols1Len[i-1]; k++) {
                        newMapCols[newMapPtr-1][k] = mapCols1[i-1][k];
                    }
                    newMapColsLen[newMapPtr-1] = mapCols1Len[i-1];
                } else {
                    // merge columns
                    newMapCols[newMapPtr-1] = (int *)malloc(g1->groupNum * sizeof(int));
                    for (int k = 0; k < mapCols1Len[i-1]; k++) {
                        newMapCols[newMapPtr-1][k] = mapCols1[i-1][k];
                    }
                    for (int k = 0; k < mapCols2Len[j-1]; k++) {
                        newMapCols[newMapPtr-1][k+mapCols1Len[i-1]] = mapCols2[j-1][k];
                    }
                    newMapColsLen[newMapPtr-1] = mapCols1Len[i-1] + mapCols2Len[j-1];
                }
                newMapPtr ++;
            }
        } else if (lastPoint == 2) {
            j++;
            if (j >= map2Len) {
                break;
            }
            if (map2[j] < map1[i]) {
                newMap[newMapPtr] = map2[j];
                newMapCols[newMapPtr-1] = (int *)malloc(g1->groupNum * sizeof(int));
                for (int k = 0; k < mapCols2Len[j-1]; k++) {
                    newMapCols[newMapPtr-1][k] = mapCols2[j-1][k];
                }
                newMapColsLen[newMapPtr-1] = mapCols2Len[j-1];
                newMapPtr ++;
            } else{
                if (map1[i] == map2[j])
                    lastPoint = 3;
                else{
                    lastPoint = 1;
                }
                newMap[newMapPtr] = map1[i];
                if (i == 0) {
                    newMapCols[newMapPtr-1] = (int *)malloc(g1->groupNum * sizeof(int));
                    for (int k = 0; k < mapCols2Len[j-1]; k++) {
                        newMapCols[newMapPtr-1][k] = mapCols2[j-1][k];
                    }
                    newMapColsLen[newMapPtr-1] = mapCols2Len[j-1];
                } else {
                    // merge columns
                    newMapCols[newMapPtr-1] = (int *)malloc(g1->groupNum * sizeof(int));
                    for (int k = 0; k < mapCols2Len[j-1]; k++) {
                        newMapCols[newMapPtr-1][k] = mapCols2[j-1][k];
                    }
                    for (int k = 0; k < mapCols1Len[i-1]; k++) {
                        newMapCols[newMapPtr-1][k+mapCols2Len[j-1]] = mapCols1[i-1][k];
                    }
                    newMapColsLen[newMapPtr-1] = mapCols1Len[i-1] + mapCols2Len[j-1];
                }
                newMapPtr ++;
            }
        } else {
            //last Point = 3
            i++;
            j++;
            if ((i >= map1Len) || (j >= map2Len)) {
                break;
            }
            if (map1[i] < map2[j]) {
                lastPoint = 1;
                newMap[newMapPtr] = map1[i];
            } else if (map1[i] == map2[j]) {
                newMap[newMapPtr] = map1[i];
            } else {
                lastPoint = 3;
                newMap[newMapPtr] = map2[j];
            }
            // merge columns
            newMapCols[newMapPtr-1] = (int *)malloc(g1->groupNum * sizeof(int));
            for (int k = 0; k < mapCols1Len[i-1]; k++) {
                newMapCols[newMapPtr-1][k] = mapCols1[i-1][k];
            }
            for (int k = 0; k < mapCols2Len[j-1]; k++) {
                newMapCols[newMapPtr-1][k+mapCols1Len[i-1]] = mapCols2[j-1][k];
            }
            newMapColsLen[newMapPtr-1] = mapCols1Len[i-1] + mapCols2Len[j-1];

            newMapPtr ++;
        }

    }

    if (i >= map1Len) {
        if (lastPoint != 1) {
            j++;
        }
        while (j < map2Len) {
            newMap[newMapPtr] = map2[j];
            newMapCols[newMapPtr-1] = (int *)malloc(g1->groupNum * sizeof(int));
            for (int k = 0; k < mapCols2Len[j-1]; k++) {
                newMapCols[newMapPtr-1][k] = mapCols2[j-1][k];
            }
            newMapColsLen[newMapPtr-1] = mapCols2Len[j-1];
            newMapPtr ++;
            j++;
        }
    } else if (j >= map2Len) {
        if (lastPoint != 2) {
            i++;
        }
        while (i < map1Len) {
            newMap[newMapPtr] = map1[i];
            newMapCols[newMapPtr-1] = (int *)malloc(g1->groupNum * sizeof(int));
            for (int k = 0; k < mapCols1Len[i-1]; k++) {
                newMapCols[newMapPtr-1][k] = mapCols1[i-1][k];
            }
            newMapColsLen[newMapPtr-1] = mapCols1Len[i-1];
            newMapPtr ++;
            i++;
        }
    }

    // destroy
    for (int idx = 0; idx < g1->rangeMapLen - 1; idx++) {
        free(g1->rangeMapCols[idx]);
    }
    free(g1->rangeMap);
    free(g1->rangeMapCols);
    free(g1->rangeMapColsLen);
    destroyApprColumnGroup(g2);

    // assign value
    g1->rangeMap = newMap;
    g1->rangeMapCols = newMapCols;
    g1->rangeMapColsLen = newMapColsLen;
    g1->rangeMapLen = newMapPtr;

    g1->timeSeriesLength = timeSeriesLength;
}

static int computeGain(ApprColumnGroup *g1, ApprColumnGroup *g2) {
    int overlap = estimateOverlap(g1, g2);
    // col - col
    int gain;
    if (g1->columnLen == 1 && g2->columnLen == 1) {
        gain = overlap * sizeof(long) - 2 * (g1->timeSeriesLength + g2->timeSeriesLength - overlap);
        return gain;
    }
    if (g1->columnLen == 1) {
        // col - group
        int m_s_a = g1->timeSeriesLength;
        int n_g_a = g2->columnLen;
        int m_g_a = g2->timeSeriesLength;
        gain = overlap * sizeof(long) + n_g_a * m_g_a - (n_g_a + 1) * (m_s_a + m_g_a - overlap);
        return gain;
    }
    if (g2->columnLen == 1) {
        // group - col
        int m_s_a = g2->timeSeriesLength;
        int n_g_a = g1->columnLen;
        int m_g_a = g1->timeSeriesLength;
        gain = overlap * sizeof(long) + n_g_a * m_g_a - (n_g_a + 1) * (m_s_a + m_g_a - overlap);
        return gain;
    }
    // group - group
    int n_g_a = g1->columnLen;
    int m_g_a = g1->timeSeriesLength;
    int n_g_b = g2->columnLen;
    int m_g_b = g2->timeSeriesLength;
    gain = overlap * sizeof(long) + (n_g_a + n_g_b) * overlap - n_g_a * m_g_b - n_g_b * m_g_a;

    return gain;
}

static long lcm(long number1, long number2) {
    if (number1 == 0 || number2 == 0) {
        return 0;
    }
    long absNumber1 = number1 >= 0 ? number1:-number1;
    long absNumber2 = number2 >= 0 ? number2:-number2;
    long absHigherNumber = absNumber1 > absNumber2 ? absNumber1 : absNumber2;
    long absLowerNumber = absNumber1 < absNumber2 ? absNumber1 : absNumber2;
    long lcm = absHigherNumber;
    while (lcm % absLowerNumber != 0) {
        lcm += absHigherNumber;
    }
    return lcm;
}

static double NormSDist(double z) {
    // compute approximate normal distribution cdf F(z)
    if (z > 6) return 1;
    if (z < -6) return 0;
    double gamma = 0.231641900,
            a1 = 0.319381530,
            a2 = -0.356563782,
            a3 = 1.781477973,
            a4 = -1.821255978,
            a5 = 1.330274429;
    double x = z > 0 ? z : -z;
    double t = 1 / (1 + gamma * x);
    double n = 1 - (1 / (sqrt(2 * PI)) * exp(-z * z / 2))
              * (a1 * t + a2 * pow(t, 2) + a3 * pow(t, 3) + a4 * pow(t, 4) + a5 * pow(t, 5));
    if (z < 0) return 1.0 - n;
    return n;
}


static float computeProbability(Feature *f1, Feature *f2) {
    if (f1->epsilon == f2->epsilon && f1->n == f2->n && f1->sigma == f2->sigma && f1->start == f2->start) {
        return 1;
    }
    if ((f1->start > f2->start + f2->n * f2->epsilon) || (f2->start > f1->start + f1->n * f1->epsilon)) {
        return 0;
    }

    int lambda = 1;
    int tau = 10;
    float prob = 0;
    long lcmInterval = lcm(f1->epsilon, f2->epsilon);
    long scenarios;
    if (f1->epsilon > f2->epsilon) {
        scenarios = lcmInterval / f1->epsilon;
        for (int i = 0; i < scenarios; i++) {
            long delta = ((f1->epsilon - f2->epsilon) * i + f1->start - f2->start) % f2->epsilon;
            if (f1->sigma == 0 && f2->sigma == 0) {
                if (delta == 0) {
                    prob += 1;
                } else {
                    prob += 0;
                }
            } else {
                float z1 = ((lambda * tau) - delta) / sqrt(f1->sigma * f1->sigma + f2->sigma * f2->sigma);
                float z2 = (-(lambda * tau) - delta) / sqrt(f1->sigma * f1->sigma + f2->sigma * f2->sigma);
                prob += NormSDist(z1) -NormSDist(z2);
            }
        }
    } else {
        scenarios = lcmInterval / f2->epsilon;
        for (int i = 0; i < scenarios; i++) {
            long delta = ((f2->epsilon - f1->epsilon) * i + f2->start - f1->start) % f1->epsilon;
            if (f1->sigma == 0 && f2->sigma == 0) {
                if (delta == 0) {
                    prob += 1;
                } else {
                    prob += 0;
                }
            } else {
                float z1 = ((lambda * tau) - delta) / sqrt(f1->sigma * f1->sigma + f2->sigma * f2->sigma);
                float z2 = (-(lambda * tau) - delta) / sqrt(f1->sigma * f1->sigma + f2->sigma * f2->sigma);
                prob += NormSDist(z1) -NormSDist(z2);
            }
        }
    }
    prob /= (float)scenarios;
    return prob;
}


static int estimateOverlap(ApprColumnGroup *g1, ApprColumnGroup *g2) {
    float *overlapList = (float *)malloc(sizeof(float) * g1->valuesLen);
    int overlapListLen = 0;

    for (int i = 0; i < g2->columnLen; i++) {
        Feature* f2 = g2->features[i];
        long start2 = f2->start;
        long end2 = f2->start + f2->epsilon * (f2->n - 1);
        float overlapTmp = 0;

        int t = 0;

        while ((t < g1->rangeMapLen - 1) && (g1->rangeMap[t] < start2)) {
            t++;
        }

        if (t == g1->rangeMapLen - 1) {
            overlapList[overlapListLen] = 0.0;
            overlapListLen ++;
            continue;
        } else {
            if ((t > 0) && (g1->rangeMap[t-1] < start2)) {
                long n2Tmp = (start2 - g1->rangeMap[t-1]) / f2->epsilon;
                float prob = 1;
                for (int colIdx = 0; colIdx < g1->rangeMapLen; colIdx++) {
                    for (int k = 0; k < g1->columnLen; k++) {
                        if (g1->columns[k] == g1->rangeMapCols[t-1][colIdx]) {
                            Feature *f1 = g1->features[k];
                            prob *= (1 - computeProbability(f1, f2));
                            break;
                        }
                    }
                }
                prob -= 1;
                overlapTmp += prob * n2Tmp;
            }
        }
        if (g1->rangeMap[t] >= start2) {
            int j = t;
            while ((j < g1->rangeMapLen - 1) && (g1->rangeMap[j] < end2)) {
                j ++;
                long n2Tmp;
                if (g1->rangeMap[j] < end2) {
                    n2Tmp = (g1->rangeMap[j] - g1->rangeMap[j-1]) / f2->epsilon;
                } else {
                    n2Tmp = (end2 - g1->rangeMap[j-1]) / f2->epsilon;
                }
                float prob = 1;
                for (int colIdx = 0; colIdx < g1->rangeMapLen; colIdx++) {
                    for (int k = 0; k < g1->columnLen; k++) {
                        if (g1->columns[k] == g1->rangeMapCols[j-1][colIdx]) {
                            Feature *f1 = g1->features[k];
                            prob *= (1 - computeProbability(f1, f2));
                            break;
                        }
                    }
                }
                prob = 1 - prob;
                overlapTmp += prob * n2Tmp;
            }
        }
        overlapList[overlapListLen] = overlapTmp;
        overlapListLen++;
    }

    float sumOverlap = 0;
    for (int i = 0; i < overlapListLen; i++) {
        sumOverlap += overlapList[i];
    }
    int sumN = 0;
    for (int i = 0; i < g2->columnLen; i++) {
        sumN += g2->features[i]->n;
    }

    free(overlapList);
    int overlapEstimation = (int) ((sumOverlap / sumN) * g2->timeSeriesLength);
    return overlapEstimation;
}

void groupingAppr(long *timestamps, float **dataColumns, int numOfRows, int numOfCols, int *columnSchemaIdx) {
    clock_t startTime = clock();
    ApprColumnGroup **columnGroups = (ApprColumnGroup **)malloc(numOfCols * sizeof(ApprColumnGroup *));
    // init columnGroups
    int pCol;
    int groupNum = numOfCols; // current number of groups
    for (pCol = 0; pCol < groupNum; pCol++) {
        columnGroups[pCol] = (ApprColumnGroup *)malloc(sizeof(ApprColumnGroup));
        CreateApprColumnGroup(columnGroups[pCol], pCol, groupNum, dataColumns[pCol], numOfRows, timestamps);
    }

    // init gain matrix
    for (pCol = 0; pCol < groupNum; pCol++) {
        columnGroups[pCol]->gains[pCol] = 0;
        for (int rCol = pCol + 1; rCol < groupNum; rCol ++) {
            int gain = computeGain(columnGroups[pCol], columnGroups[rCol]);
            columnGroups[pCol]->gains[rCol] = gain;
            columnGroups[rCol]->gains[pCol] = gain;
            if (columnGroups[pCol]->maxGain < gain) {
                columnGroups[pCol]->maxGain = gain;
                columnGroups[pCol]->maxGainPosition = rCol;
            }
            if (columnGroups[rCol]->maxGain < gain) {
                columnGroups[rCol]->maxGain = gain;
                columnGroups[rCol]->maxGainPosition = pCol;
            }
        }
    }

    while (groupNum > 1) {
        /** merge **/
        // find the max gain
        int maxGain = -1;
        int maxGainPos = -1;
        for (pCol = 0; pCol < groupNum; pCol++) {
            if (columnGroups[pCol]->maxGain > maxGain) {
                maxGain = columnGroups[pCol]->maxGain;
                maxGainPos = pCol;
            }
        }
        // no positive gain, exit
        if (maxGain <= 0) {
            break;
        }
        // merge the group
        int source = maxGainPos;
        int target = columnGroups[source]->maxGainPosition;
        // printf("source: %d, target: %d\n", source, target);
        mergeApprColumnGroup(columnGroups[source], columnGroups[target]);

        // remove old groups
        for (pCol = target; pCol < groupNum - 1; pCol++) {
            columnGroups[pCol] = columnGroups[pCol+1];
        }
        groupNum--;

        /** update gains **/
        ApprColumnGroup *gSource = columnGroups[source];
        for (pCol = 0; pCol < groupNum; pCol++) {
            if (pCol != source) {
                ApprColumnGroup *g = columnGroups[pCol];
                // remove target gain from pCol
                for (int rCol = target; rCol < groupNum; rCol++) {
                    g->gains[rCol] = g->gains[rCol+1];
                }
                int gain = computeGain(g, gSource);
                g->gains[source] = gain;
                // update g maxgain
                if ((g->maxGainPosition != source) && (g->maxGainPosition != target)) {
                    if (gain > g->maxGain) {
                        g->maxGain = gain;
                        g->maxGainPosition = source;
                    } else {
                        if (target < g->maxGainPosition) {
                            g->maxGainPosition -= 1;
                        }
                    }
                } else {
                    g->maxGain = gain;
                    g->maxGainPosition = source;
                    for (int rCol = 0; rCol < groupNum; rCol ++) {
                        if (g->gains[rCol] > g->maxGain) {
                            g->maxGain = g->gains[rCol];
                            g->maxGainPosition = rCol;
                        }
                    }
                }
                // update the maxgain and data in source
                gSource->gains[pCol] = gain;
            }
        }
        // update maxgain of the source
        gSource->gains[source] = 0;
        gSource->maxGain = -1;
        for (pCol = 0; pCol < groupNum; pCol ++) {
            if (gSource->gains[pCol] > gSource->maxGain) {
                gSource->maxGain = gSource->gains[pCol];
                gSource->maxGainPosition = pCol;
            }
        }
    }

    clock_t endTime = clock();
    // output results to files
    outputResults(columnGroups, groupNum, columnSchemaIdx);
    outputGroupingToFile(columnGroups, APPR_GROUP_FILE, groupNum, columnSchemaIdx);
    double timeCost = (double) (endTime-startTime) / CLOCKS_PER_SEC;
    outputTimeToFile(APPR_TIME_FILE, timeCost);

    // free memory and exit
    for (pCol = 0; pCol < groupNum; pCol++){
        destroyApprColumnGroup(columnGroups[pCol]);
    }
    free(columnGroups);
}

static int cmpInt(const void* e1, const void* e2)
{
    return (int)(*(int*)e1 - *(int*)e2);
}


static void outputResults(ApprColumnGroup **columnGroups, int groupNum, int *columnSchemaIdx) {
    printf("-------Approximate algorithm results-------\n");
    for (int pCol = 0; pCol < groupNum; pCol++) {
        qsort(columnGroups[pCol]->columns, columnGroups[pCol]->columnLen, sizeof(int), cmpInt);
        for (int idx = 0; idx < columnGroups[pCol]->columnLen; idx++) {
            int tmpIdx = columnSchemaIdx[columnGroups[pCol]->columns[idx]];
            printf("c%d ", tmpIdx);
        }
        printf("\n");
    }
}

static void outputGroupingToFile(ApprColumnGroup **columnGroups, const char *filename, int groupNum, int *columnSchemaIdx) {
    FILE *fp = fopen(filename, "a+");
    for (int pCol = 0; pCol < groupNum; pCol++) {
        qsort(columnGroups[pCol]->columns, columnGroups[pCol]->columnLen, sizeof(int), cmpInt);
        for (int idx = 0; idx < columnGroups[pCol]->columnLen; idx++) {
            int tmpIdx = columnSchemaIdx[columnGroups[pCol]->columns[idx]];
            fprintf(fp, "c%d", tmpIdx);
            if (idx < columnGroups[pCol]->columnLen - 1){
                fprintf(fp, ",");
            }
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "-----------\n");
    fclose(fp);
}

static void outputTimeToFile(const char *filename, double t) {
    FILE *fp = fopen(filename, "a+");
    fprintf(fp, "%f\n", t);
    fclose(fp);
}