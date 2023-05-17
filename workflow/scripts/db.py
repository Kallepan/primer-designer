import sqlite3

class DBHandler():
    def setup_alignments_table(self):
        with self.con:
            self.con.execute(
                """
                CREATE TABLE IF NOT EXISTS primer_alignments (
                    pool INT NOT NULL,
                    primer_id TEXT NOT NULL,
                    strand TEXT NOT NULL,
                    chromosome TEXT NOT NULL,
                    position INT NOT NULL,
                    sequence TEXT NOT NULL,
                    read_quality TEXT NOT NULL,
                    matches INT NOT NULL,
                    mismatches_descriptor TEXT,

                    FOREIGN KEY (primer_id) REFERENCES proto_primers (primer_id)
                );
                """
            )

            self.con.execute(
                """
                CREATE INDEX IF NOT EXISTS primer_alignments_index ON primer_alignments (pool, primer_id, position, matches);
                """
            )

    def setup_proto_primers_table(self):
        with self.con:
            self.con.execute("""
                CREATE TABLE IF NOT EXISTS proto_primers(
                    primer_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    pool INT NOT NULL,
                    region_name TEXT NOT NULL,
                    amplicon_name TEXT NOT NULL,
                    strand TEXT NOT NULL,
                    primer_sequence TEXT NOT NULL,
                    primer_length INT,
                    tm REAL,
                    gc_percent REAL,
                    hairpin_th REAL,
                    badness REAL,
                    UNIQUE(pool, region_name, amplicon_name, strand, primer_sequence)
                )
            """)
            self.con.execute("CREATE INDEX IF NOT EXISTS idx_primers ON proto_primers(pool, region_name, amplicon_name, strand)")
            self.con.execute("CREATE INDEX IF NOT EXISTS idx_primers_badness ON proto_primers(badness)")

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
        
    def __del__(self):
        self.con.close()