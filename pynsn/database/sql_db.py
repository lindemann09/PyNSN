import sqlite3

class DotArraySQLDB(object):

    def __init__(self, db_file, setup=False):
        self.db_file = db_file
        if setup:
            self.setup()

    def setup(self):
        conn = sqlite3.connect(self.db_file)

        cursor = conn.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = [v[0] for v in cursor.fetchall() if v[0] != "sqlite_sequence"]

        if 'ARRAYS' not in tables:
            conn.execute('''CREATE TABLE ARRAYS
                             (HASH varchar(32) PRIMARY KEY,
                             N   INT     NOT NULL,
                             TSA REAL NOT NULL,
                             ISA  REAL NOT NULL,
                             FA  REAL NOT NULL,
                             SPAR  REAL NOT NULL,
                             logSIZE REAL NOT NULL,
                             logSPACE REAL NOT NULL,
                             COV REAL NOT NULL,
                             COLOUR TEXT);''')

        if 'DOTS' not in tables:
            conn.execute('''CREATE TABLE DOTS
                     (HASH varchar(32) NOT NULL,
                     x REAL NOT NULL,
                     y REAL NOT NULL,
                     diameter REAL NOT NULL,
                     COLOUR TEXT);''')

        cursor.close()
        conn.close()

    def add_arrays(self, dot_arrays):

        if not isinstance(dot_arrays, (list, tuple)):
            dot_arrays = [dot_arrays]

        conn = sqlite3.connect(self.db_file)
        cur = conn.cursor()

        for da in dot_arrays:
            attributes = da.as_dict()['attributes']
            # TODO if multi color, add to dots color column. multi_colour = len(attributes)!=1

            # assume single colour
            colour = attributes[0]["colour"]

            sql = "INSERT INTO ARRAYS (" + \
                  "HASH, N, TSA, ISA, FA, SPAR, logSIZE, logSPACE, COV, COLOUR" + \
                  ") \n VALUES\n" + \
                  "('{}',{},{},{},{},{},{},{},{},'{}');".format(
                          da.hash,
                          da.features.numerosity,
                          da.features.total_surface_area,
                          da.features.mean_item_surface_area,
                          da.features.field_area,
                          da.features.sparsity,
                          da.features.logSize,
                          da.features.logSpacing,
                          da.features.converage,
                          colour)
            cur.execute(sql)

            ## add dots
            sql = "INSERT INTO DOTS (HASH,x,y,diameter) \nVALUES"
            for xy, d in zip(da.xy, da.diameters):
                sql += "\n  ('{}', {}, {}, {}),".format(da.hash, xy[0], xy[1], d)
            sql = sql[:-1] + ";"
            cur.execute(sql)

        conn.commit()
        conn.close()

    def get_all_hashes(self):
        conn = sqlite3.connect(self.db_file)
        values = conn.execute("SELECT hash FROM ARRAYS;").fetchall()
        conn.close()
        return [x[0] for x in values]

    def get_array(self, hash):
        conn = sqlite3.connect(self.db_file)
        conn.row_factory = sqlite3.Row
        cur = conn.execute("SELECT * FROM ARRAYS WHERE hash ='{}';".format(hash))
        v = cur.fetchone()
        conn.close()
        return {k:v[k] for k in v.keys()}

    def get_dots(self, hash):
        # get ID
        conn = sqlite3.connect(self.db_file)
        values = conn.execute("SELECT x, y, diameter, colour FROM DOTS WHERE hash ='{}';".format(hash)).fetchall()
        conn.close()
        return values #

