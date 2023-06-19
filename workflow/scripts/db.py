import sqlite3
import pandas as pd

class DBHandler:
    def setup_regions_table(self):
        with self.con:
            self.con.execute(
                """
                CREATE TABLE IF NOT EXISTS regions (
                    name TEXT PRIMARY KEY NOT NULL,
                    start INT NOT NULL,
                    end INT NOT NULL
                );
                """
            )
            self.con.execute(
                """CREATE INDEX IF NOT EXISTS idx_regions_name ON regions(name);"""
            )

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
                    score REAL DEFAULT 0.0,
                    read_quality TEXT NOT NULL,
                    mismatches_descriptor TEXT,
                    species TEXT NOT NULL,

                    FOREIGN KEY (primer_id) REFERENCES proto_primers (id)
                );
                """
            )

            self.con.execute(
                """CREATE INDEX IF NOT EXISTS idx_alignments_id ON alignments(id)"""
            )
            self.con.execute(
                """CREATE INDEX IF NOT EXISTS idx_alignments_pool ON alignments(pool)"""
            )
            self.con.execute(
                """CREATE INDEX IF NOT EXISTS idx_alignments_primer_id ON alignments(primer_id)"""
            )
            self.con.execute(
                """CREATE INDEX IF NOT EXISTS idx_alignments_aligned_to ON alignments(aligned_to)"""
            )
            self.con.execute(
                """CREATE INDEX IF NOT EXISTS idx_alignments_position ON alignments(position)"""
            )
            self.con.execute(
                """CREATE INDEX IF NOT EXISTS idx_alignments_multi ON alignments(pool, id, position, matches)"""
            )

    def setup_proto_primers_table(self):
        with self.con:
            self.con.execute(
                """
                CREATE TABLE IF NOT EXISTS proto_primers(
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    pool INT NOT NULL,
                    region_name TEXT NOT NULL,
                    amplicon_name TEXT NOT NULL,
                    strand TEXT NOT NULL,
                    sequence TEXT NOT NULL,
                    length INT NOT NULL,
                    tm REAL NOT NULL,
                    gc_percent REAL NOT NULL,
                    hairpin_th REAL NOT NULL,
                    discarded BOOLEAN NOT NULL DEFAULT FALSE,
                    badness REAL NOT NULL DEFAULT 0.0,

                    UNIQUE(pool, region_name, amplicon_name, strand, sequence),
                    FOREIGN KEY (region_name) REFERENCES regions (name)
                )
            """
            )

            self.con.execute(
                "CREATE INDEX IF NOT EXISTS idx_primers_id ON proto_primers(id)"
            )
            self.con.execute(
                "CREATE INDEX IF NOT EXISTS idx_primers_sequence ON proto_primers(sequence)"
            )
            self.con.execute(
                "CREATE INDEX IF NOT EXISTS idx_primers_pool ON proto_primers(pool)"
            )
            self.con.execute(
                "CREATE INDEX IF NOT EXISTS idx_primers_amplicon_name ON proto_primers(amplicon_name)"
            )
            self.con.execute(
                "CREATE INDEX IF NOT EXISTS idx_primers_multi ON proto_primers(pool, region_name, amplicon_name, strand)"
            )

    def __init__(self, path_to_db: str) -> None:
        con = sqlite3.connect(path_to_db, timeout=60)
        con.execute("PRAGMA journal_mode=WAL")
        con.execute("PRAGMA busy_timeout = 60000")
        self.con = con

    def get_tables(self) -> list[str]:
        """Returns a list of tables in the database"""
        try:
            with self.con:
                return [
                    x[0]
                    for x in self.con.execute(
                        "SELECT name FROM sqlite_master WHERE type='table'"
                    ).fetchall()
                ]
        except Exception as e:
            raise e

    def select(self, *args, **kwargs) -> tuple[list, list]:
        """Returns a tuple of (rows, column_names)"""
        try:
            with self.con:
                return self.con.execute(*args, **kwargs).fetchall(), [
                    x[0] for x in self.con.execute(*args, **kwargs).description
                ]
        except Exception as e:
            raise e

    def execute(self, *args, **kwargs) -> sqlite3.Cursor:
        """Executes a query"""
        try:
            with self.con:
                return self.con.execute(*args, **kwargs)
        except Exception as e:
            raise e

    def executemany(self, *args, **kwargs) -> sqlite3.Cursor:
        """Executes a query with many values"""
        try:
            with self.con:
                return self.con.executemany(*args, **kwargs)
        except Exception as e:
            raise e

    def get_columns(self, table_name: str) -> list[str]:
        """Returns a list of column names for a given table"""
        try:
            with self.con:
                return [
                    x[1]
                    for x in self.con.execute(
                        f"PRAGMA table_info({table_name})"
                    ).fetchall()
                ]
        except Exception as e:
            raise e

    def __del__(self):
        self.con.close()
