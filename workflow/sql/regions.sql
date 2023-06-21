CREATE TABLE IF NOT EXISTS regions (
    name TEXT PRIMARY KEY NOT NULL,
    start INT NOT NULL,
    end INT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_regions_name ON regions(name);