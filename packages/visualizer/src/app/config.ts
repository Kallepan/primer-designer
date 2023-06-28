//
// Stores all constants used in the application in a config object
//

export class CONFIG {
    public static PRIMER_COLOR = '#FF0000';
    public static COLORS = [
        '#FFB722',
        '#D8FF22',
        '#6AFF22',
        '#22FF48',
        '#22FFB7',
        '#22D8FF',
        '#EF22FF',
        '#2234FF',
    ];

    // Configurations for the detailed view
    public static DETAILED_VIEW = {
        SEQUENCE_LENGTH_LIMIT: 30_000,
        SEQUENCE_VISIBLE_LENGTH: 3_000,
        BATCH_SIZE: 1000,
        PRIMER_COLOR: this.PRIMER_COLOR,
    };    

    // Configurations for the simplified view
    public static SIMPLIFIED_VIEW = {
        PRIMER_OFFSET: 10,
        SCALE_FACTOR: 1000,
    };

    public static MESSAGES = {
        NO_RESULTS: 'No results found for the given region.',
        NO_REGION_INFO: 'No region info found for the given region.',
    };
};