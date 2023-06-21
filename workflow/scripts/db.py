"""Handles database connections and queries"""
import sqlite3


class DBHandler:
    def setup_table(self, path_to_sql_file: str):
        with open(path_to_sql_file, "r") as f:
            sql = f.read()
        with self.con:
            self.con.executescript(sql)

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
