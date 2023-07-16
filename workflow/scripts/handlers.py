"""Handles database connections and queries"""
import sqlite3


class DBHandler:
    def setup_table(self, path_to_sql_file: str):
        with open(path_to_sql_file, "r") as f:
            sql = f.read()
        with self.conn:
            self.conn.executescript(sql)

    def __init__(self, path_to_db: str) -> None:
        # Connect to database and increase timeout
        conn = sqlite3.connect(path_to_db, timeout=60)
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA busy_timeout = 60000")
        self.conn = conn

    def get_tables(self) -> list[str]:
        """Returns a list of tables in the database"""
        try:
            with self.conn:
                return [
                    x[0]
                    for x in self.conn.execute(
                        "SELECT name FROM sqlite_master WHERE type='table'"
                    ).fetchall()
                ]
        except Exception as e:
            raise e

    def select(self, *args, **kwargs) -> tuple[list, list]:
        """Returns a tuple of (rows, column_names)"""
        try:
            with self.conn:
                return self.conn.execute(*args, **kwargs).fetchall(), [
                    x[0] for x in self.conn.execute(*args, **kwargs).description
                ]
        except Exception as e:
            raise e

    def execute(self, *args, **kwargs) -> sqlite3.Cursor:
        """Executes a query"""
        try:
            with self.conn:
                return self.conn.execute(*args, **kwargs)
        except Exception as e:
            raise e

    def executemany(self, *args, **kwargs) -> sqlite3.Cursor:
        """Executes a query with many values"""
        try:
            with self.conn:
                return self.conn.executemany(*args, **kwargs)
        except Exception as e:
            raise e

    def get_columns(self, table_name: str) -> list[str]:
        """Returns a list of column names for a given table"""
        try:
            with self.conn:
                return [
                    x[1]
                    for x in self.conn.execute(
                        f"PRAGMA table_info({table_name})"
                    ).fetchall()
                ]
        except Exception as e:
            raise e

    def __del__(self):
        try:
            with self.conn:
                self.conn.execute("PRAGMA optimize")
        except Exception as e:
            raise e
        finally:
            self.conn.close()
