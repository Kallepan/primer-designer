# Program
- written in python (snakemake)

# Requirements
## Amplicon
    - 1000 bp TODO
## Primers
    - no appends necessary
    - reduce hairpin
    - melting temperature
    - test for intra and inter primer (pair) interaction
    - option to test against other strains (binding capability) exclude unwanted targets

    > checked in evaluation function
    > calculated for each primer pair
    > complexity: O(n!)
# Procedure
## Workflow
1. Generate DNA Targets
    - Download genome of target species
    - take points of interest as an input
    - generate desired amplicon regions +- K bp Amplicon size
2. Generate PCR Primers (see Algorigthms)
    - generate upper primers
    - generate lower primers
3. Evaluate initial set of Primers
    - Randomly select pcr targets -> Set x
    - Set composition:
        - even numbers: forward
        - odd numbers: reverse
    - Evaluation function for whole set
        - compare each primer with each other -> B(p1, p1), B(p1, p2) ...
        - L(Sx) = SUM(B(pa,pb)) | a >= b
    - set current_set = initial_set
4. Generate Temporary Set
    - replace one or more primers from the current set by random replacement
5. Calculate new set
```python
if L(temporary_set) <= L(current_set)
    current_set  = temporary_set
else
    current_set = current_set     (with probability 1-p)
    current_set = temporary_set   (with probability p)

probability = C(g, L(temporary_set) - L(current_set))
```
6. Repeat 4-5 until convergence


# Algorithms
## Primer generation
1. Proto Primers (Adatapted form Oli2go)
    - determine pivot point in genome to calculate primers in
    - determine amplicon size in upper and lower direction
    - specificity check by blasting against selected ncbi databases
    - gc between 35 und 80%
    - tm between 55 and 65
    - accepted properties in primer design:
        - melting temperature around X degress
        - dG values around -10.5 to -12.5 kcal/mol
2. Generate Primer sets from given primers
    - associate primers in pairs
3. Calculate Badness
    - sum of all badness values of all primer possible pairs
    - badness between two primers:complementary sequences
        - **(2^len * 2^numGC)/((d1+1)(d2+1))**
            - d = distance between subsequence and 3' end of primer
            - len = length of subsequence
            - numGC = number of GC bases in subsequence
        - max subsequence length of 8
    - generate all possible sub sequences of length 4 to 8 nt
    - determine if the reverse complement of the subsequence is present in the other primer
    - calculate the badness value for each subsequence and store the distance in a hashtable as **Sum(1/d+1)**
    - When looking the badness up: **2^len * 2^numGC / (1+d) * HashTable[reverse_complement(subsequence)]**
4. Generate Temp Primer set
    - replace one or more primers from the current set by random replacement
    - replacement of more than one primer set is possible but can lead to slowdown
5. Calculate new set
    - set new set to be either temp set or old set depending on badness and probability of switching
```python
if L(temporary_set) <= L(current_set):
    current_set  = temporary_set
else:
    current_set = current_set     #(with probability 1-p)
    current_set = temporary_set   #(with probability p)
```
probability = C(g, L(temporary_set) - L(current_set))
p = exp(-L(temporary_set) - L(current_set)) / C(g)
g -> number of generations
at some point stop simulated annealing and use gradient descent

6. repeat steps 4 to 5 until convergence
    - repeat calculation for 1.5 times the number of generations to reach convergence

# Questions
- Was ist die optimale Ampliconlänge? 1000 bp? user defined.

- In Saddle werden statt thermodynamische Parameter verinfachte Parameter verwendet. Kann ich das hier auch so machen?

- Problematisch bei der Methode: Target Spezifität wird nicht garantiert. Etwas vorne dranbasteln um das zu entgehen.

- Primer3/Primer BLAST einbinden? Als Test? Oder sogar als initialer Generator

- Programmierung in Python oder C++ oder Go? 