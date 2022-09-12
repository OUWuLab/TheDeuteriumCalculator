# Input Files (Only change the path, not the names or text after the '#' leave the r before the quotes.)
IDENTIFICATION_MZML_FILE = r"D:\Projects\ManuscriptPreparation\2019_Software_development\Final\Supplementary Files\Input Files\Antigen ID.mzML"  # must be .mzML
IDENTIFICATION_CSV_FILE = r"D:\Projects\ManuscriptPreparation\2019_Software_development\Final\Supplementary Files\Input Files\Antigen Peptides Library.csv"  # Must be .csv
PROTEIN_SEQUENCE_FILE = r"D:\Projects\ManuscriptPreparation\2019_Software_development\Final\Supplementary Files\Input Files\Antigen Sequence.txt"  # Must be .txt


# Output Files (Only change the path, not the names or text after the '#' leave the r before the quotes.)
FULL_HDX_OUTPUT = r"D:\Projects\ManuscriptPreparation\2019_Software_development\Final\Test\output"  # No file extension
RECOMMENDATION_TABLE_1 = r"D:\Projects\ManuscriptPreparation\2019_Software_development\Final\Test\t1"  # No file extension
RECOMMENDATION_TABLE_2 = r"D:\Projects\ManuscriptPreparation\2019_Software_development\Final\Test\t2"  # No file extension
SUMMARY_TABLE = r"D:\Projects\ManuscriptPreparation\2019_Software_development\Final\Test\sum"  # No file extension
WOODS_PLOT_NAME = r"D:\Projects\ManuscriptPreparation\2019_Software_development\Final\Test\f1"  # No file extension
WOODS_TABLE_NAME = r"D:\Projects\ManuscriptPreparation\2019_Software_development\Final\Test\f1"  # No file extension


# Parameters (Only change the numbers and path, not the names or text after the '#'. Defaults given in parenthesis.)
# Step 1 Parameters
NOISE_LIMIT = 10000  # All individual peaks with less than this intensity ignored (10000)
PPM_MATCH_TOLERANCE = 10  # Peaks with less difference than this value matched to sequence (10)
SLIDING_WINDOW_PPM_TOLERANCE = 1  # Peaks with less difference than this value combined within each sliding window (1)
SLIDING_WINDOW_SIZE = 30  # width of sliding window in seconds, should be integer divisible by SLIDE_FRACTION (60)
SLIDE_FRACTION = 3  # Fraction of the window that the window moves each each slide (3)
RETENTION_TOLERANCE = 30  # window of retention times to search for given peptide (+-) (30)
DEUTERIUM_RECOVERY_RATE = 1  # The experimentally determined back exchange rate. (1)
DEUTERIUM_FRACTION = 0.833  # The fraction of D2O used in the experiment. (1)
CONDITION1 = "Free"  # This will be used in single condition experiments
CONDITION2 = "Complex"
# Step 2 Parameters
WOODS_PLOT_CONFIDENCE_HIGH = 0.99  # Use to calculate confidence interval for differential woods plot (0-1)
WOODS_PLOT_CONFIDENCE_LOW = 0.98  # Use to calculate confidence interval for differential woods plot (0-1)
WOODS_PLOT_TITLE = r"Differential Uptake Woods' Plot"  # Title on the woods' plot
FRACTIONAL_WOODS_PLOT_TITLE = r"Differential Fractional Uptake Woods' Plot"  # Title on the fractional woods' plot
WOODS_PLOT_HEIGHT = 4  # Inches  (4)
WOODS_PLOT_WIDTH = 5  # Inches  (5)
RETENTION_SHIFT_SLOPE = 1
RETENTION_SHIFT_INTERCEPT = 0

# Constants (Change at your own risk; may cause incorrect results)
DEUTERIUM_MASS_DIFFERENCE = 1.00627
MASS_OF_WATER = 18.01528
MASS_OF_HYDROGEN = 1.007276
MINUTES_TO_SECONDS = 60
PEPTIDE_MASS_DICTIONARY = {
    'A': 71.0779,
    'C': 103.1429,
    'D': 115.0874,
    'E': 129.114,
    'F': 147.1739,
    'G': 57.0513,
    'H': 137.1393,
    'I': 113.1576,
    'K': 128.1723,
    'L': 113.1576,
    'M': 131.1961,
    'N': 114.1026,
    'P': 97.1152,
    'Q': 128.1292,
    'R': 156.1857,
    'S': 87.0773,
    'T': 101.1039,
    'U': 150.0379,
    'W': 186.2099,
    'Y': 163.1733,
    'V': 99.1311
}
