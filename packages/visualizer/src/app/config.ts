//
// Stores all constants used in the application in a config object
//

export class CONFIG {
    public static readonly PRIMER_COLOR = '#FF0000';
    public static readonly COLORS = ["#b30000", "#7c1158", "#4421af", "#1a53ff", "#0d88e6", "#00b7c7", "#5ad45a", "#8be04e", "#ebdc78"];
    
    // Configurations for the detailed view
    public static DETAILED_VIEW = {
        SEQUENCE_MAX_LENGTH: 50_000,
        DEFAULT_COLOR: '#FFFFFF',
        PRIMER_COLOR: this.PRIMER_COLOR,
    };    

    // Configurations for the simplified view
    public static SIMPLIFIED_VIEW = {
        PRIMER_OFFSET: 50,
        PLOT_WIDTH: 1200,

        PLOT_MARGINS: {
            TOP: 50,
            RIGHT: 20,
            BOTTOM: 50,
            LEFT: 70,
        },
    };

    // Configurations for the saddle view
    public static SADDLE_VIEW = {
        PLOT_WIDTH: 1200,
        PLOT_HEIGHT: 600,

        PLOT_MARGINS: {
            TOP: 20,
            RIGHT: 30,
            BOTTOM: 30,
            LEFT: 40,
        },
    };

    public static MESSAGES = {
        NO_RESULTS: 'No results found for the given region.',
        NO_REGION_INFO: 'No region info found for the given region.',
    };
};