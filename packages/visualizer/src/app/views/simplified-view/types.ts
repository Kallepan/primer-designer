/*
    Custom types for the simplified view
*/
export type SimplifiedPrimerData = {
    pool: string;
    id: string;
    x1: number;
    x2: number;

    sequence: string;
    score: number;
};

export type SimplifiedAmpliconData = {
    name: string;
    pool: string;
    x1: number;
    x2: number;
};

export type SimplifiedRegionData = {
    name: string;
    start: number;
    end: number;
    primers: SimplifiedPrimerData[];
    amplicons: SimplifiedAmpliconData[];
};