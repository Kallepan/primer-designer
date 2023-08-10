CREATE TABLE IF NOT EXISTS amplicons (
    region_name TEXT NOT NULL,
    name TEXT NOT NULL,
    start INTEGER NOT NULL,x
    end INTEGER NOT NULL,
    pool TEXT NOT NULL,
    failed BOOLEAN NOT NULL DEFAULT FALSE,

    FOREIGN KEY (region_name) REFERENCES regions (name)
);

CREATE INDEX IF NOT EXISTS amplicons_region_name_idx ON amplicons(region_name);
CREATE INDEX IF NOT EXISTS amplicons_amplicon_name_idx ON amplicons(amplicon_name);
CREATE INDEX IF NOT EXISTS amplicons_amplicon_id_idx ON amplicons(amplicon_id);
