import sqlite3
from ..image._colour import Colour, ImageColours


class DotArraySQLDB(object):
    """DotArray DB (two tables: _arrays, dots)"""

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

    def vacuum(self):
        conn = sqlite3.connect(self.db_file)
        conn.execute("VACUUM;")
        conn.close()

    def add_arrays(self, dot_arrays):

        if not isinstance(dot_arrays, (list, tuple)):
            dot_arrays = [dot_arrays]

        conn = sqlite3.connect(self.db_file)
        cur = conn.cursor()

        for da in dot_arrays:
            attributes = da.as_dict()['attributes']

            # assuming single colour dot _arrays
            # TODO if multi color, add to dots color column. multi_colour = len(attributes)!=1
            try:
                colour = Colour(attributes[0])
            except:
                # if attribute is no colour, use default coolour
                colour = Colour(ImageColours.COL_DEFAULT_OBJECT)

            sql = "INSERT INTO ARRAYS (" + \
                  "HASH, N, TSA, ISA, FA, SPAR, logSIZE, logSPACE, COV, COLOUR" + \
                  ") \n VALUES\n" + \
                  "('{}',{},{},{},{},{},{},{},{},'{}');".format(
                          da.hash,
                          da._properties.numerosity,
                          da._properties.total_surface_area,
                          da._properties.average_surface_area,
                          da._properties.field_area,
                          da._properties.sparsity,
                          da._properties.log_size,
                          da._properties.log_spacing,
                          da._properties.coverage,
                          colour.colour)
            cur.execute(sql)

            ## add dots
            sql = "INSERT INTO DOTS (HASH,x,y,diameter) \nVALUES"
            for xy, d in zip(da._xy, da._diameter):
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

