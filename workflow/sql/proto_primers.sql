CREATE TABLE IF NOT EXISTS proto_primers (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    pool INT NOT NULL,
    region_name TEXT NOT NULL,
    amplicon_name TEXT NOT NULL,
    strand TEXT NOT NULL,
    sequence TEXT NOT NULL,
    tm REAL NOT NULL,
    gc_percent REAL NOT NULL,
    hairpin_th REAL NOT NULL,
    discarded BOOLEAN NOT NULL DEFAULT FALSE,
    position INT NOT NULL,

    UNIQUE(pool, region_name, amplicon_name, strand, sequence),
    FOREIGN KEY (region_name) REFERENCES regions (name)
);

CREATE INDEX IF NOT EXISTS idx_proto_primers_pool_idx ON proto_primers (pool);
CREATE INDEX IF NOT EXISTS idx_proto_primers_id ON proto_primers (id);
CREATE INDEX IF NOT EXISTS idx_proto_primers_amplicon_name ON proto_primers (amplicon_name);
CREATE INDEX IF NOT EXISTS idx_proto_primers_region_name ON proto_primers (region_name);
CREATE INDEX IF NOT EXISTS idx_proto_primers_strand ON proto_primers (strand);
CREATE INDEX IF NOT EXISTS idx_proto_primers_position ON proto_primers (position);
CREATE INDEX IF NOT EXISTS idx_proto_primers_sequence ON proto_primers (sequence);

CREATE INDEX IF NOT EXISTS idx_proto_primers_multi_one ON proto_primers (pool, region_name, amplicon_name, strand);