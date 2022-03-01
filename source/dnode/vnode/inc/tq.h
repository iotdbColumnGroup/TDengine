/*
 * Copyright (c) 2019 TAOS Data, Inc. <jhtao@taosdata.com>
 *
 * This program is free software: you can use, redistribute, and/or modify
 * it under the terms of the GNU Affero General Public License, version 3
 * or later ("AGPL"), as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _TD_TQ_H_
#define _TD_TQ_H_

#include "tcommon.h"
#include "executor.h"
#include "tmallocator.h"
#include "meta.h"
#include "scheduler.h"
#include "taoserror.h"
#include "tlist.h"
#include "tmsg.h"
#include "trpc.h"
#include "ttimer.h"
#include "tutil.h"
#include "vnode.h"
#include "wal.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct STQ STQ;

// memory allocator provided by vnode
typedef struct {
  SMemAllocatorFactory* pAllocatorFactory;
  SMemAllocator*        pAllocator;
} STqMemRef;

// init once
int  tqInit();
void tqCleanUp();

// open in each vnode
STQ* tqOpen(const char* path, SWal* pWal, SMeta* pMeta, STqCfg* tqConfig, SMemAllocatorFactory* allocFac);
void tqClose(STQ*);

// required by vnode
int tqPushMsg(STQ*, void* msg, int64_t version);
int tqCommit(STQ*);

int32_t tqProcessConsumeReq(STQ* pTq, SRpcMsg* pMsg);
int32_t tqProcessSetConnReq(STQ* pTq, char* msg);
int32_t tqProcessRebReq(STQ* pTq, char* msg);

#ifdef __cplusplus
}
#endif

#endif /*_TD_TQ_H_*/
