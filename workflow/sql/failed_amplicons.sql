CREATE TABLE IF NOT EXISTS failed_amplicons (
    region_name TEXT NOT NULL,
    amplicon_name TEXT NOT NULL,
    amplicon_id TEXT NOT NULL,
    amplicon_start: INTEGER NOT NULL,
    amplicon_end: INTEGER NOT NULL,
);

CREATE INDEX IF NOT EXISTS failed_amplicons_region_name_idx ON failed_amplicons(region_name);
CREATE INDEX IF NOT EXISTS failed_amplicons_amplicon_name_idx ON failed_amplicons(amplicon_name);
CREATE INDEX IF NOT EXISTS failed_amplicons_amplicon_id_idx ON failed_amplicons(amplicon_id);
