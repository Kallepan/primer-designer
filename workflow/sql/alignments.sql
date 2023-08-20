CREATE TABLE IF NOT EXISTS alignments (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    primer_id INT NOT NULL,
    aligned_to TEXT NOT NULL,
    chromosome TEXT NOT NULL,
    position INT NOT NULL,
    sequence TEXT NOT NULL,
    matches INT NOT NULL,
    mismatches_descriptor TEXT,
    species TEXT NOT NULL,

    FOREIGN KEY (primer_id) REFERENCES proto_primers (id)
);

CREATE INDEX IF NOT EXISTS idx_alignments_id_idx ON alignments(id);
CREATE INDEX IF NOT EXISTS idx_alignments_primer_id_idx ON alignments(primer_id);
CREATE INDEX IF NOT EXISTS idx_alignments_aligned_to_idx ON alignments(aligned_to);
CREATE INDEX IF NOT EXISTS idx_alignments_position_idx ON alignments(position);
CREATE INDEX IF NOT EXISTS idx_alignments_chromosome_idx ON alignments(chromosome);
CREATE INDEX IF NOT EXISTS idx_alignments_species_idx ON alignments(species);
CREATE INDEX IF NOT EXISTS idx_alignments_matches_idx ON alignments(matches);
CREATE INDEX IF NOT EXISTS idx_alignments_mismatches_descriptor_idx ON alignments(mismatches_descriptor);