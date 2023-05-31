import sqlite3

class DBHandler():
    def setup_alignments_table(self):
        with self.con:
            self.con.execute(
                """
                CREATE TABLE IF NOT EXISTS alignments (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    pool TEXT NOT NULL,
                    primer_id INT NOT NULL,
                    aligned_to TEXT NOT NULL,
                    chromosome TEXT NOT NULL,
                    position INT NOT NULL,
                    sequence TEXT NOT NULL,
                    matches INT NOT NULL,
                    score REAL,
                    read_quality TEXT NO NULL,
                    mismatches_descriptor TEXT,

                    FOREIGN KEY (primer_id) REFERENCES proto_primers (id)
                );
                """
            )

            self.con.execute("""CREATE INDEX IF NOT EXISTS idx_alignments_id ON alignments(id)""")
            self.con.execute("""CREATE INDEX IF NOT EXISTS idx_alignments_pool ON alignments(pool)""")
            self.con.execute("""CREATE INDEX IF NOT EXISTS idx_alignments_primer_id ON alignments(primer_id)""")
            self.con.execute("""CREATE INDEX IF NOT EXISTS idx_alignments_aligned_to ON alignments(aligned_to)""")
            self.con.execute("""CREATE INDEX IF NOT EXISTS idx_alignments_multi ON alignments(pool, id, position, matches)""")

    def setup_proto_primers_table(self):
        with self.con:
            self.con.execute("""
                CREATE TABLE IF NOT EXISTS proto_primers(
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    pool INT NOT NULL,
                    region_name TEXT NOT NULL,
                    amplicon_name TEXT NOT NULL,
                    strand TEXT NOT NULL,
                    sequence TEXT NOT NULL,
                    length INT,
                    tm REAL,
                    gc_percent REAL,
                    hairpin_th REAL,
                    badness REAL,
                    UNIQUE(pool, region_name, amplicon_name, strand, sequence)
                )
            """)
            
            self.con.execute("CREATE INDEX IF NOT EXISTS idx_primers_id ON proto_primers(id)")
            self.con.execute("CREATE INDEX IF NOT EXISTS idx_primers_sequence ON proto_primers(sequence)")
            self.con.execute("CREATE INDEX IF NOT EXISTS idx_primers_pool ON proto_primers(pool)")
            self.con.execute("CREATE INDEX IF NOT EXISTS idx_primers_multi ON proto_primers(pool, region_name, amplicon_name, strand)")

    def __init__(self, path_to_db: str) -> None:
        con = sqlite3.connect(path_to_db)
        self.con = con

    def select(self, *args, **kwargs) -> list:
        try:
            with self.con:
                return self.con.execute(*args, **kwargs).fetchall()
        except Exception as e:
            raise e

    def execute(self, *args, **kwargs) -> sqlite3.Cursor:
        try:
            with self.con:
                return self.con.execute(*args, **kwargs)
        except Exception as e:
            raise e
    
    def executemany(self, *args, **kwargs) -> sqlite3.Cursor:
        try:
            with self.con:
                return self.con.executemany(*args, **kwargs)
        except Exception as e:
            raise e

    def __del__(self):
        self.con.close()