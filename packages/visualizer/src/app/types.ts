//
// Type definitions to import the results json file
//
export type Primer = {
    id: string;
    sequence: string;
    tm: number;
    gc_percent: number;
    hairpin_th: number;
    length: number;
    position: number;

    // Needed for visualization
    relativePosition?: number;
};
export type PrimerPair = {
    region_name: string;
    amplicon_name: string;
    forward_primer: Primer;
    reverse_primer: Primer;

    // Needed for visualization
    pool_id?: string;
    color?: string;
};
export type Pool = {
    primer_pairs: PrimerPair[];
    pool_id: string;
    loss: number;
};

//
// Type definitions to import the regions json file
//
export type RegionInfo = {
    name: string;
    start: number;
    end: number;
    sequence: string;
};

//
// Type definitions to import loss json file
//
export type Loss = {
    pool_id: string;
    losses: number[];
};

//
// Type definitions of the data used for visualization
//
export type Region = {
    name: string;
    start: number;
    end: number;
    sequence: string;
    primerPairsByPool: Map<string, PrimerPair[]>;
    primersByPool: Map<string, Primer[]>;
};

//
// Type definitions to import amplicons json file
//
export type Amplicon = {
    name: string;
    region: string;
    pool: string;
    discarded: boolean;

    n_forward?: number;
    n_reverse?: number;
    n_discarded_forward?: number;
    n_discarded_reverse?: number;
};

export type FormattedAmplicon = {
    name: string;
    region: string;
    pool: string;
    discarded: boolean;


    forwardPrimers: string;
    discardedForwardPrimers: string;
    reversePrimers: string;
    discardedReversePrimers: string;
}