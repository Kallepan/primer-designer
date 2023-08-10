CREATE TABLE IF NOT EXISTS amplicons (
    region_name TEXT NOT NULL,
    name TEXT PRIMARY KEY,
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    pool TEXT NOT NULL,
    failed BOOLEAN NOT NULL DEFAULT FALSE,

    FOREIGN KEY (region_name) REFERENCES regions (name)
);

CREATE INDEX IF NOT EXISTS amplicons_region_name_idx ON amplicons(region_name);
CREATE INDEX IF NOT EXISTS amplicons_name_idx ON amplicons(name);
CREATE INDEX IF NOT EXISTS amplicons_pool_idx ON amplicons(pool);
CREATE INDEX IF NOT EXISTS amplicons_failed_idx ON amplicons(failed);