CREATE TABLE IF NOT EXISTS regions (
    name TEXT PRIMARY KEY,
    start INT NOT NULL,
    end INT NOT NULL,
    sequence TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_regions_name ON regions(name);
CREATE INDEX IF NOT EXISTS idx_regions_start ON regions(start);
CREATE INDEX IF NOT EXISTS idx_regions_end ON regions(end);