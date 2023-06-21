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
    position INT NOT NULL DEFAULT -1, -- Position in the region

    UNIQUE(pool, region_name, amplicon_name, strand, sequence),
    FOREIGN KEY (region_name) REFERENCES regions (name)
);

CREATE INDEX IF NOT EXISTS idx_proto_primers_pool_idx ON proto_primers (pool);
CREATE INDEX IF NOT EXISTS idx_proto_primers_id ON proto_primers (id);
CREATE INDEX IF NOT EXISTS idx_proto_primers_amplicon_name ON proto_primers (amplicon_name);
CREATE INDEX IF NOT EXISTS idx_proto_primers_region_name ON proto_primers (region_name);
CREATE INDEX IF NOT EXISTS idx_proto_primers_multi ON proto_primers (pool, region_name, amplicon_name, strand);