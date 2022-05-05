import taos
from util.log import *
from util.cases import *
from util.sql import *
import numpy as np


class TDTestCase:
    def init(self, conn, logSql):
        tdLog.debug("start to execute %s" % __file__)
        tdSql.init(conn.cursor())

        self.rowNum = 10
        self.ts = 1537146000000  # 2018-9-17 09:00:00.000
        
    def run(self):
        tdSql.prepare()

        intData = []        
        floatData = []

        tdSql.execute('''create table stb(ts timestamp, col1 tinyint, col2 smallint, col3 int, col4 bigint, col5 float, col6 double, 
                    col7 bool, col8 binary(20), col9 nchar(20), col11 tinyint unsigned, col12 smallint unsigned, col13 int unsigned, col14 bigint unsigned) tags(loc nchar(20))''')
        tdSql.execute("create table stb_1 using stb tags('beijing')")
        tdSql.execute('''create table ntb(ts timestamp, col1 tinyint, col2 smallint, col3 int, col4 bigint, col5 float, col6 double, 
                    col7 bool, col8 binary(20), col9 nchar(20), col11 tinyint unsigned, col12 smallint unsigned, col13 int unsigned, col14 bigint unsigned)''')
        for i in range(self.rowNum):
            tdSql.execute("insert into ntb values(%d, %d, %d, %d, %d, %f, %f, %d, 'taosdata%d', '涛思数据%d', %d, %d, %d, %d)" 
                        % (self.ts + i, i + 1, i + 1, i + 1, i + 1, i + 0.1, i + 0.1, i % 2, i + 1, i + 1, i + 1, i + 1, i + 1, i + 1))
            intData.append(i + 1)            
            floatData.append(i + 0.1)
        for i in range(self.rowNum):
            tdSql.execute("insert into stb_1 values(%d, %d, %d, %d, %d, %f, %f, %d, 'taosdata%d', '涛思数据%d', %d, %d, %d, %d)" 
                        % (self.ts + i, i + 1, i + 1, i + 1, i + 1, i + 0.1, i + 0.1, i % 2, i + 1, i + 1, i + 1, i + 1, i + 1, i + 1))
            intData.append(i + 1)            
            floatData.append(i + 0.1)  

        tdSql.query("select timetruncate(1,1d) from ntb")
        tdSql.checkRows(10)
        tdSql.query("select timetruncate(1,1u) from ntb")
        tdSql.checkRows(10)
        tdSql.query("select timetruncate(1,1a) from ntb")
        tdSql.checkRows(10)
        tdSql.query("select timetruncate(1,1m) from ntb")
        tdSql.checkRows(10)
        tdSql.query("select timetruncate(1,1h) from ntb")
        tdSql.checkRows(10)
        tdSql.query("select timetruncate(ts,1d) from ntb")
        tdSql.checkRows(10)
        tdSql.checkData(0,0,"2018-09-17 08:00:00.000")
        tdSql.query("select timetruncate(ts,1h) from ntb")
        tdSql.checkRows(10)
        tdSql.checkData(0,0,"2018-09-17 09:00:00.000")
        tdSql.query("select timetruncate(ts,1m) from ntb")
        tdSql.checkRows(10)
        tdSql.checkData(0,0,"2018-09-17 09:00:00.000")
        tdSql.query("select timetruncate(ts,1s) from ntb")
        tdSql.checkRows(10)
        tdSql.checkData(0,0,"2018-09-17 09:00:00.000")
        tdSql.query("select timetruncate(ts,1a) from ntb")
        tdSql.checkRows(10)
        tdSql.checkData(0,0,"2018-09-17 09:00:00.000")
        tdSql.checkData(1,0,"2018-09-17 09:00:00.001")
        tdSql.checkData(2,0,"2018-09-17 09:00:00.002")
        tdSql.checkData(3,0,"2018-09-17 09:00:00.003")
        tdSql.checkData(4,0,"2018-09-17 09:00:00.004")
        tdSql.checkData(5,0,"2018-09-17 09:00:00.005")
        tdSql.checkData(6,0,"2018-09-17 09:00:00.006")
        tdSql.checkData(7,0,"2018-09-17 09:00:00.007")
        tdSql.checkData(8,0,"2018-09-17 09:00:00.008")
        tdSql.checkData(9,0,"2018-09-17 09:00:00.009")
        # tdSql.query("select timetruncate(ts,1u) from ntb")
        # tdSql.checkRows(10)
        # tdSql.checkData(0,0,"2018-09-17 09:00:00.000000")
        # tdSql.checkData(1,0,"2018-09-17 09:00:00.001000")
        # tdSql.checkData(2,0,"2018-09-17 09:00:00.002000")
        # tdSql.checkData(3,0,"2018-09-17 09:00:00.003000")
        # tdSql.checkData(4,0,"2018-09-17 09:00:00.004000")
        # tdSql.checkData(5,0,"2018-09-17 09:00:00.005000")
        # tdSql.checkData(6,0,"2018-09-17 09:00:00.006000")
        # tdSql.checkData(7,0,"2018-09-17 09:00:00.007000")
        # tdSql.checkData(8,0,"2018-09-17 09:00:00.008000")
        # tdSql.checkData(9,0,"2018-09-17 09:00:00.009000")
        # tdSql.query("select timetruncate(ts,1b) from ntb")
        # tdSql.checkRows(10)
        # tdSql.checkData(0,0,"2018-09-17 09:00:00.000000000")
        # tdSql.checkData(1,0,"2018-09-17 09:00:00.001000000")
        # tdSql.checkData(2,0,"2018-09-17 09:00:00.002000000")
        # tdSql.checkData(3,0,"2018-09-17 09:00:00.003000000")
        # tdSql.checkData(4,0,"2018-09-17 09:00:00.004000000")
        # tdSql.checkData(5,0,"2018-09-17 09:00:00.005000000")
        # tdSql.checkData(6,0,"2018-09-17 09:00:00.006000000")
        # tdSql.checkData(7,0,"2018-09-17 09:00:00.007000000")
        # tdSql.checkData(8,0,"2018-09-17 09:00:00.008000000")
        # tdSql.checkData(9,0,"2018-09-17 09:00:00.009000000")


          
    def stop(self):
        tdSql.close()
        tdLog.success("%s successfully executed" % __file__)

tdCases.addWindows(__file__, TDTestCase())
tdCases.addLinux(__file__, TDTestCase())