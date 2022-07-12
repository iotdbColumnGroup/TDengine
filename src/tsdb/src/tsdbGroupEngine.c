#include "tsdbGroupEngine.h"
#include "stdio.h"
#include "time.h"

#define EXACT_GROUP_FILE "/workspaces/TDengine/results/group_result_exact.csv"
#define EXACT_TIME_FILE "/workspaces/TDengine/results/time_exact.csv"

static void outputGroupingToFile(ColumnGroup **columnGroups, const char *filename, int groupNum, int *columnSchemaIdx);
static void outputTimeToFile(const char *filename, double t);
static void CreateColumnGroup(ColumnGroup *columnGroup, int cid, int groupNum, float *values, int valuesLen);
static void updateTimestampInfo(ColumnGroup *columnGroup);
static void destroyColumnGroup(ColumnGroup *columnGroup);
static void outputResults(ColumnGroup **columnGroups, int groupNum, int *columnSchemaIdx);
static int computeGain(ColumnGroup *g1, ColumnGroup *g2);
static int computeOverlap(ColumnGroup *g1, ColumnGroup *g2);

static int isFloatNull(float *val) {
    return *(__uint32_t *)val == FLOAT_NULL;
}

static void CreateColumnGroup(ColumnGroup *columnGroup, int cid, int groupNum, float *values, int valuesLen) {
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
    updateTimestampInfo(columnGroup);
}

static void updateTimestampInfo(ColumnGroup *columnGroup) {
    // count number of not-nan values
    columnGroup->timeSeriesLength = 0;
    for (int i = 0; i < columnGroup->valuesLen; i++) {
        if (!isFloatNull(&columnGroup->values[i])) {
            columnGroup->timeSeriesLength++;
        }
    }
}

static void destroyColumnGroup(ColumnGroup *columnGroup) {
    free(columnGroup->gains);
    free(columnGroup->columns);
    free(columnGroup->values);
    free(columnGroup);
}

static void mergeColumnGroup(ColumnGroup *g1, ColumnGroup *g2) {
    /** merge g2 values into g1 and update **/
    // merge columns
    for (int i = 0; i < g2->columnLen; i++) {
        g1->columns[i + g1->columnLen] = g2->columns[i];
    }
    g1->columnLen = g1->columnLen + g2->columnLen;
    // merge values
    for (int i = 0; i < g1->valuesLen; i++) {
        if (isFloatNull(&g1->values[i])) {
            g1->values[i] = g2->values[i];
        }
    }
    updateTimestampInfo(g1);
    // destroy
    destroyColumnGroup(g2);
}

static int computeGain(ColumnGroup *g1, ColumnGroup *g2) {
    int overlap = computeOverlap(g1, g2);
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

static int computeOverlap(ColumnGroup *g1, ColumnGroup *g2) {
    int overlap = 0;
    for (int pRow = 0; pRow < g1->valuesLen; pRow++) {
        if (!isFloatNull(&(g1->values[pRow])) && !isFloatNull(&(g2->values[pRow]))) {
            overlap++;
        }
    }
    return overlap;
}

void grouping(long *timestamps, float **dataColumns, int numOfRows, int numOfCols, int *columnSchemaIdx) {
    clock_t startTime = clock();
    ColumnGroup **columnGroups = (ColumnGroup **)malloc(numOfCols * sizeof(ColumnGroup *));
    // init columnGroups
    int pCol;
    int groupNum = numOfCols; // current number of groups
    for (pCol = 0; pCol < groupNum; pCol++) {
        columnGroups[pCol] = (ColumnGroup *)malloc(sizeof(ColumnGroup));
        CreateColumnGroup(columnGroups[pCol], pCol, groupNum, dataColumns[pCol], numOfRows);
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
        mergeColumnGroup(columnGroups[source], columnGroups[target]);

        // remove old groups
        for (pCol = target; pCol < groupNum - 1; pCol++) {
            columnGroups[pCol] = columnGroups[pCol+1];
        }
        groupNum--;

        /** update gains **/
        ColumnGroup *gSource = columnGroups[source];
        for (pCol = 0; pCol < groupNum; pCol++) {
            if (pCol != source) {
                ColumnGroup *g = columnGroups[pCol];
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
    outputGroupingToFile(columnGroups, EXACT_GROUP_FILE, groupNum, columnSchemaIdx);
    double timeCost = (double) (endTime-startTime) / CLOCKS_PER_SEC;
    outputTimeToFile(EXACT_TIME_FILE, timeCost);

    // free memory and exit
    for (pCol = 0; pCol < groupNum; pCol++){
        destroyColumnGroup(columnGroups[pCol]);
    }
    free(columnGroups);
}

static int cmpInt(const void* e1, const void* e2)
{
    return (int)(*(int*)e1 - *(int*)e2);
}


static void outputResults(ColumnGroup **columnGroups, int groupNum, int *columnSchemaIdx) {
    printf("-------Exact algorithm results-------\n");
    for (int pCol = 0; pCol < groupNum; pCol++) {
        qsort(columnGroups[pCol]->columns, columnGroups[pCol]->columnLen, sizeof(int), cmpInt);
        for (int idx = 0; idx < columnGroups[pCol]->columnLen; idx++) {
            int tmpIdx = columnSchemaIdx[columnGroups[pCol]->columns[idx]];
            printf("c%d ", tmpIdx);
        }
        printf("\n");
    }
}

static void outputGroupingToFile(ColumnGroup **columnGroups, const char *filename, int groupNum, int *columnSchemaIdx) {
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