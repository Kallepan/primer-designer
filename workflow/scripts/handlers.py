""" This module contains the handlers for the workflow scripts. """
import sqlite3

from collections import defaultdict


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

    def __del__(self) -> None:
        if self.conn is None:
            return

        try:
            with self.conn:
                self.conn.execute("PRAGMA optimize")
        except Exception as e:
            # Ignore if the database is locked
            if "database is locked" not in str(e):
                raise e
        finally:
            self.conn.close()


class Graph(object):
    """Graph data structure, undirected by default."""

    def __init__(self, connections, directed=False) -> None:
        self._graph = defaultdict(set)
        self._directed = directed
        self.__add_connections(connections)

    def __add_connections(self, connections) -> None:
        """Add connections (list of tuple pairs) to graph"""

        for node1, node2 in connections:
            self.add(node1, node2)

    def add(self, node1, node2) -> None:
        """Add connection between node1 and node2"""

        self._graph[node1].add(node2)
        if not self._directed:
            self._graph[node2].add(node1)

    def remove(self, node) -> None:
        """Remove all references to node"""

        for n, cxns in self._graph.items():  # python3: items(); python2: iteritems()
            try:
                cxns.remove(node)
            except KeyError:
                pass
        try:
            del self._graph[node]
        except KeyError:
            pass

    def is_connected(self, node1, node2) -> bool:
        """Is node1 directly connected to node2"""

        return node1 in self._graph and node2 in self._graph[node1]

    def find_path(self, node1, node2, path=[]) -> list[object]:
        """Find any path between node1 and node2 (may not be shortest)"""

        path = path + [node1]
        if node1 == node2:
            return path
        if node1 not in self._graph:
            return None
        for node in self._graph[node1]:
            if node not in path:
                new_path = self.find_path(node, node2, path)
                if new_path:
                    return new_path
        return None

    def get_all_edges(self) -> list[tuple[object, object]]:
        """Get all edges of the graph"""
        edges = []
        for node1 in self._graph:
            for node2 in self._graph[node1]:
                edges.append((node1, node2))
        return edges

    def vertex_cover_approx(self) -> list[object]:
        """
        Calculate the approximate vertex cover of the graph and returns the vertices to be removed.
        """
        edges = self.get_all_edges()
        visited = defaultdict(bool)

        # Consider all edges one by one
        for u in self._graph.keys():
            # An edge is only picked when both visited[u] and visited[v] are false
            if visited[u]:
                continue

            # Go through all adjacents of u and pick the first not visited vertex
            for v in self._graph[u]:
                if not visited[v]:
                    visited[u] = True
                    visited[v] = True
                    break

        return [x for x in visited if visited[x]]

    def __str__(self) -> str:
        """String representation of the graph"""
        return "{}({})".format(self.__class__.__name__, dict(self._graph))


### Tests ###
def test_graph() -> None:
    """Test the graph class and the vertex cover algorithm"""
    ### Test 1 ###
    connections = [
        (0, 1),
        (0, 2),
        (1, 3),
        (3, 4),
        (4, 5),
        (4, 6),
        (5, 6),
        (5, 7),
        (6, 7),
    ]

    graph = Graph(connections, directed=False)
    print(graph)
    vertex = graph.vertex_cover_approx()
    print(vertex)

    ### Test 2 ###
    connections = [
        (0, 1),
        (1, 2),
        (2, 3),
    ]
    graph = Graph(connections, directed=False)
    print(graph)
    vertex = graph.vertex_cover_approx()
    print(vertex)

    ### Test 3 ###
    connections = [(1, 2), (1, 4), (1, 5), (2, 3), (4, 6)]
    graph = Graph(connections, directed=False)
    print(graph)
    vertex = graph.vertex_cover_approx()
    print(vertex)
