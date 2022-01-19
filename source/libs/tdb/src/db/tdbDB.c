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

#if 0
#include "tdbDB.h"
#include "tdb.h"

TDB_EXTERN int tdbCreateDB(TDB** dbpp, tdb_db_t type) {
  TDB* dbp;
  int  ret;

  dbp = calloc(1, sizeof(*dbp));
  if (dbp == NULL) {
    return -1;
  }

  dbp->pageSize = TDB_DEFAULT_PGSIZE;
  dbp->type = type;

  switch (type) {
    case TDB_BTREE_T:
      //   ret = tdbInitBtreeDB(dbp);
      //   if (ret < 0) goto _err;
      break;
    case TDB_HASH_T:
      //   ret = tdbInitHashDB(dbp);
      //   if (ret < 0) goto _err;
      break;
    case TDB_HEAP_T:
      //   ret = tdbInitHeapDB(dbp);
      //   if (ret < 0) goto _err;
      break;
    default:
      break;
  }

  *dbpp = dbp;
  return 0;

_err:
  if (dbp) {
    free(dbp);
  }
  *dbpp = NULL;
  return 0;
}

TDB_EXTERN int tdbOpenDB(TDB* dbp, const char* fname, const char* dbname, uint32_t flags) {
  int ret = 0;

  if ((dbp->fname = strdup(fname)) == NULL) {
    ret = -1;
    return ret;
  }

  // Create the backup file if the file not exists

  // Open the file as a sub-db or a master-db
  if (dbname) {
    if ((dbp->dbname = strdup(dbname)) == NULL) {
      ret = -1;
      return ret;
    }
    // TODO: Open the DB as a SUB-DB in this file
  } else {
    // TODO: Open the DB as a MASTER-DB in this file
  }

  return ret;
}

TDB_EXTERN int tdbCloseDB(TDB* dbp, uint32_t flags) {
  // TODO
  return 0;
}
#endif