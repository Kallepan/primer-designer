# Workflow

## 0. Preprocessing
- Create an UI in which the user can see targetable regions in the genome.
- The user should be able to select regions of interest and the target organism.

## 1. Proto-Primers Generation

- Generate all possible primers using Primer3
- Allow the user to change primer3 config settings (tm, primer length, gc content)
- Simple wrapper written in python to call primer3
- Generation:
    - Two pools
    - Define regions of interest which are divided into amplicons
    - Each amplicons is flanked by forward and reverse primers
- Store output in json file
- Generate as many proto-primers as possible

### Algorithm
- Generate optimal coordinates for amplicons:
    - Pool1: AMP1-BUFFER-AMP2-BUFFER-AMP3-BUFFER-AMP4...
    - Pool2: XXAMP1-BUFFER-AMP2-BUFFER-AMP3-BUFFER-AMP4...
    - Basically pool2 is the same as pool1, but with a X bp overlap between amplicons.
- Generate primers for each amplicon using primer3:
    - Forward primer: 5'->3'
    - Reverse primer: 3'->5' (outputed in 5'->3' direction)
    If no primers could be placed apply some pre-defined logic:
        - Increase Amplicon X bases in both directions up to a certain amount
        - Increase the are in which the primers can be placed (e.g.: 0-20 is now 0-30)
        - Skip the generation of the primer for this amplicon
- Store the primers in a json file along with meta information

## 2. Filter Proto-Primers
- Filter the primers using bowtie:
    The basic idea is to align all primers to the reference genome and filter out primer pairs! that have multiple alignments in the target organism or one alignment in a genome to be excluded.

- Pre-Processing:
    - Generate Fasta file from the primer json.
    - Determine which primers are to be included in the fasta file.
    - If an amplicon only has one primer in a direction, then do not included that one.
    - Generate a reference index file for the reference fasta file.
    - Indexes to include
        - human -> negative
        - target organism -> ensure specificity (e.g.: only allow one match or matches with a certain distance above a threshold)

- Processing:
    1. Find all primers not aligning to the original site:
        - Primers aligning multiple times
        - Primers with mismatches
        - Primers aligning to the wrong strand
    2. For each primer, find all primers aligning to the other strand within a certain distance
        - use a "hard" cutoff to remove reverse primers above a certain distance
        - use a "soft" cutoff to score primers above a certain distance
        - use another cutoff to remove problematic primers with a certain amount of mismatches
    3. Calculate a badness score for each primer considering:
        - the amount of mismatches
        - misalignments for the primer
        - the amount of adjacent primers
        - the amount of mismatches in the sequence
    4. Append the score in the filtered_proto_primers.json file and output

## 3. SADDLE
- Run the SADDLE algorithm on the filtered primers.
- SADDLE is a simulated annealing algorithm that optimizes primer sets for a given set of parameters.
### Algorithm
1. Generate a Set
2. Calculate loss for set (L = SUM(Badness(p1,p2))), where p1 and p2 are primers in the set
3. Generate a temporary set by replacing one or more primers from the current set by random replacement
4. Calculate loss for temporary set
5. Accept or reject temporary set
```python
if L(temporary_set) <= L(current_set):
    current_set  = temporary_set
else:
    current_set = current_set     #(with probability 1-p)
    current_set = temporary_set   #(with probability p)
```
6. Repeat 3-5 until convergence

## 4. Profit
- Output the final primer set to a file.


# Notes
- log primer matches which were not filtered out but have some kind of issue
- optional hard filter or just score
- Score: manuell berechnen
- Problematischer Primer zusätzlicher Filter:
    3 ' Ende mismatches sind gut -> nicht in problematische Primer
- Faktoren: 
    - Wie oft misaligned der Primer (problematic primer) (x^misalignments) wobei x < 1 und > 0
    - Wie weit entfernt liegt der Reverse Primer vom Forward Primer
        - Ab einem Cutoff betrachten wir uns die problematischen Alignments gar nicht
        - Ab X bp Abstand wird es weniger problematisch (konstant) Amplicon bp + 75 %
        - Unter X bp Abstand wird es problematisch (konstant) Amplicon bp + 25 %
    - Wie viele mismatches hat der Primer !
- Self badness wird zur SADDLE Funktion hinzugefügt