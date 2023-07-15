CREATE TABLE IF NOT EXISTS alignments (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    pool TEXT NOT NULL,
    primer_id INT NOT NULL,
    aligned_to TEXT NOT NULL,
    chromosome TEXT NOT NULL,
    position INT NOT NULL,
    sequence TEXT NOT NULL,
    matches INT NOT NULL,
    score REAL DEFAULT 0.0,
    read_quality TEXT NOT NULL,
    mismatches_descriptor TEXT,
    species TEXT NOT NULL,

    FOREIGN KEY (primer_id) REFERENCES proto_primers (id)
);

CREATE INDEX IF NOT EXISTS idx_alignments_id ON alignments(id);
CREATE INDEX IF NOT EXISTS idx_alignments_pool ON alignments(pool);
CREATE INDEX IF NOT EXISTS idx_alignments_primer_id ON alignments(primer_id);
CREATE INDEX IF NOT EXISTS idx_alignments_aligned_to ON alignments(aligned_to);
CREATE INDEX IF NOT EXISTS idx_alignments_position ON alignments(position);
CREATE INDEX IF NOT EXISTS idx_alignments_chromosome ON alignments(chromosome);
CREATE INDEX IF NOT EXISTS idx_alignments_species ON alignments(species);
CREATE INDEX IF NOT EXISTS idx_alignments_matches ON alignments(matches);
CREATE INDEX IF NOT EXISTS idx_alignments_mismatches_descriptor ON alignments(mismatches_descriptor);

CREATE INDEX IF NOT EXISTS idx_alignments_multi_one ON alignments(pool, primer_id, aligned_to, chromosome, position, matches);
CREATE INDEX IF NOT EXISTS idx_alignments_multi_two ON alignments(pool, primer_id, aligned_to, chromosome, position, matches, species);
CREATE INDEX IF NOT EXISTS idx_alignments_multi_three ON alignments(chromosome, aligned_to, position);