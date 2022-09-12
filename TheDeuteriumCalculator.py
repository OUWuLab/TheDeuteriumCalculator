import numpy as np
import pandas as pd
import csv
from scipy.optimize import curve_fit
from scipy import stats
import PARAMETERS as CON
from pyteomics import mzml
from datetime import datetime
from os import path
from matplotlib import pyplot as plt
import warnings
warnings.filterwarnings("ignore", message="Covariance of the parameters could not be estimated")
SLIDE_AMOUNT = CON.SLIDING_WINDOW_SIZE / CON.SLIDE_FRACTION


# checks whether an output directory exists, ignores actual filename
def check_directory(file_path):
    end_path = 0
    for index, char in enumerate(file_path):
        if char == '\\' or char == '/':
            end_path = index
    if end_path > 0:
        if not path.exists(file_path[0:end_path]):
            raise FileExistsError("ERROR: Output directory does not exit.")


# returns a gaussian with the given parameters
def gaussian_trend(x, a, sigma, mean):
    return a * np.exp(-(x - mean) ** 2 / (2 * sigma ** 2))


# returns the r^2 for a set of data compared to its guassian fit.
def calculate_r_squared(x_data, y_data, *p):
    y_predicted = gaussian_trend(x_data, *p)
    residuals = (y_predicted - y_data) ** 2
    ss_res = sum(residuals)
    ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)
    r_squared = round(r_squared, 3)
    return r_squared


# fits a Guassian curve to ordered list parameter, returns r^2
def fit_gaussian(y_data):
    x_data = list(range(len(y_data)))
    consecutive_non_zeroes = 0
    for y in y_data:
        if consecutive_non_zeroes == 3:
            break
        if y > 0:
            consecutive_non_zeroes += 1
        else:
            consecutive_non_zeroes = 0
    if consecutive_non_zeroes != 3:
        return "Insufficient Match"
    try:
        popt, pcov = curve_fit(gaussian_trend, x_data, y_data)
        r_squared = calculate_r_squared(x_data, y_data, *popt)
        if r_squared < 0:
            return 0
        return r_squared
    except RuntimeError:
        return "Couldn't Fit"

# returns the maximum possible deuterium for a given sequence with the consideration of Back exchange
# This is represented as the # of amino acids - 2 - the number of proline molecules * recovery rate * D2O Fraction.

def sequence_to_max_deuterium(sequence: str):
    max_deuterium = len(sequence) - 2
    for letter in sequence:
        if letter.lower() == 'p':
            max_deuterium -= 1
    return int(max_deuterium * CON.DEUTERIUM_RECOVERY_RATE * CON.DEUTERIUM_FRACTION) + 1

# returns the maximum possible deuterium for a given sequence
# This is represented as the # of amino acids - 2 - the number of proline molecules * recovery rate * D2O Fraction.
def max_deuterium_mlf(sequence: str):
    max_deuterium = len(sequence) - 2
    for letter in sequence:
        if letter.lower() == 'p':
            max_deuterium -= 1
    max_d_mlf = max_deuterium * CON.DEUTERIUM_FRACTION
    return max_d_mlf


def check_output_extension(file: str):
    for letter in file:
        if letter == '.':
            raise NameError("Output files should not have an extension, fix and restart the program.")


# Checks that user PARAMETER configuration is (relatively) correct
def check_parameters():
    output_files = [CON.RECOMMENDATION_TABLE_1, CON.RECOMMENDATION_TABLE_2, CON.SUMMARY_TABLE, CON.WOODS_PLOT_NAME,
                    CON.FULL_HDX_OUTPUT, CON.WOODS_TABLE_NAME]
    for file in output_files:
        check_output_extension(file)
        check_directory(file)
    if not check_extension(CON.IDENTIFICATION_MZML_FILE, "mzml"):
        raise NameError("IDENTIFICATION_MZML_FILE must be a .mzml, fix and restart the program.")
    if not check_extension(CON.IDENTIFICATION_CSV_FILE, "csv"):
        raise NameError("IDENTIFICATION_CSV_FILE must be a .csv, fix and restart the program.")
    if not check_extension(CON.PROTEIN_SEQUENCE_FILE, "txt"):
        raise NameError("PROTEIN_SEQUENCE_FILE must be a .txt, fix and restart the program.")
    input_files = [CON.IDENTIFICATION_MZML_FILE, CON.IDENTIFICATION_CSV_FILE, CON.PROTEIN_SEQUENCE_FILE]
    for file in input_files:
        if not path.exists(file):
            raise NameError("The following file path is incorrect: {}".format(file))
    if CON.PPM_MATCH_TOLERANCE > 100:
        warnings.warn(("PPM_MATCH_TOLERANCE is set to {}; this is high".format(CON.PPM_MATCH_TOLERANCE)))
    if CON.SLIDING_WINDOW_PPM_TOLERANCE > 10:
        warnings.warn(("SLIDING_WINDOW_PPM_TOLERANCE is set to {}, "
                       "recommendation is under 10.".format(CON.SLIDING_WINDOW_PPM_TOLERANCE)))
    if CON.SLIDING_WINDOW_SIZE % SLIDE_AMOUNT != 0:
        raise ValueError("SLIDING_WINDOW_SIZE must be divisible by SLIDE_AMOUNT")
    if CON.RETENTION_TOLERANCE < 10:
        raise ValueError("RETENTION_TOLERANCE must be a number over 10, greater than 30 recommended.")
    if CON.RETENTION_TOLERANCE < 30:
        warnings.warn("A retention tolerance of greater than 30s is recommended")
    if not 0 < CON.WOODS_PLOT_CONFIDENCE_LOW < 1:
        raise ValueError("Woods' Plot Confidence LOW must be between 0 and 1")
    if not 0 < CON.WOODS_PLOT_CONFIDENCE_HIGH < 1:
        raise ValueError("Woods' Plot Confidence HIGH must be between 0 and 1")


def check_extension(string, extension):
    # helper function for check_parameters()
    if string[-len(extension):].lower() != extension:
        return False
    else:
        return True


# Gets user input of path and changes it to usable string
def get_path_input():
    input_path = r"{}".format(input())
    compatible_path = ""
    for letter in input_path:
        if letter == '\\':
            compatible_path += '\\'
        if letter != '"':
            compatible_path += letter
    return compatible_path


# Determines the location of the sequence within the full protein
def find_start_end(peptide: str, protein: str):
    sequential_matches = 0
    start, end = 0, 0
    index = 0
    while index < len(protein):
        if peptide[sequential_matches] == protein[index]:
            sequential_matches += 1
            if sequential_matches == len(peptide):
                end = index + 1
                return start, end
            elif sequential_matches == 1:
                start = index + 1
            index += 1
        else:
            if sequential_matches != 0:
                index = start
            else:
                index += 1
            sequential_matches = 0
    return "NULL", "NULL"


# Removes any non-alpha characters from the protein sequence
def get_ppm(mz1, mz2):
    return abs(mz1 - mz2) / mz1 * 1000000


def parse_protein(file: str):
    sequence = ""
    with open(file, 'r') as f:
        file_reader = f.readlines()
        for line in file_reader:
            for character in line:
                if character.isalpha():
                    sequence += character
    return sequence
# Returns the ppm difference between two m/z values


# Takes in a list of tuples and combines each where first elements is within a ppm tolerance
# This is used to combine peaks of the exact same mz in consecutive nearby scans.
def tuple_combine(some_list):
    changed = True
    just_changed = False
    start_list = sorted(some_list, key=lambda x: x[0])
    new_list = []
    for i in range(len(start_list)):
        tup = (start_list[i][0], start_list[i][1], 1)
        start_list[i] = tup

    while changed and len(start_list) != 0:
        changed = False
        for i in range(len(start_list) - 1):
            if just_changed:
                just_changed = False
                continue
            ppm = abs(get_ppm(start_list[i][0], start_list[i+1][0]))
            if ppm < CON.SLIDING_WINDOW_PPM_TOLERANCE:
                count1 = start_list[i][2]
                count2 = start_list[i + 1][2]
                mz1 = start_list[i][0] * count1
                mz2 = start_list[i + 1][0] * count2
                count = count1 + count2
                mz = (mz1 + mz2) / count
                intensity = start_list[i][1] + start_list[i + 1][1]
                new_list.append((mz, intensity, count))
                changed = True
                just_changed = True
            else:
                new_list.append(start_list[i])
        if not just_changed:
            new_list.append(start_list[len(start_list) - 1])
        start_list = new_list
        new_list = []
        just_changed = False
    return start_list


# Modified binary search which returns the highest-intensity peak within a ppm tolerance
def compare(target, charge, array, full_array):
    midpoint = int(len(array) / 2)
    # This represents any match
    try:
        if abs(get_ppm(target, array[midpoint][0] * charge)) <= CON.PPM_MATCH_TOLERANCE:
            return_list = [(array[midpoint][0], array[midpoint][1])]
            offset = 1
            # These two while loops then check any adjacent peak
            while offset != 0 and (midpoint - offset) > 0:
                peak = array[midpoint - offset]
                ppm = abs(get_ppm(target, peak[0] * charge))
                if ppm <= CON.PPM_MATCH_TOLERANCE:
                    return_list.append((peak[0], peak[1]))
                    offset += 1
                else:
                    offset = 0
            offset = 1
            while offset != 0 and (midpoint + offset < len(array)):
                peak = array[midpoint + offset]
                ppm = abs(get_ppm(target, peak[0] * charge))
                if ppm <= CON.PPM_MATCH_TOLERANCE:
                    return_list.append((peak[0], peak[1]))
                    offset += 1
                else:
                    offset = 0
            high_intensity = (0, 0)
            for key, value in return_list:
                if value > high_intensity[1]:
                    high_intensity = key, value

            ppm_error = abs(get_ppm(target, high_intensity[0] * charge))
            if ppm_error > CON.PPM_MATCH_TOLERANCE:
                print("PPM ERROR!")
            return ppm_error, high_intensity[0],  high_intensity[1]

        elif len(array) == 1 or len(array) == 0:
            return 0, 0, 0
        elif array[midpoint][0] * charge <= target:
            return compare(target, charge, array[midpoint:], full_array)
        else:
            return compare(target, charge, array[0: midpoint], full_array)
    except IndexError:
        return 0, 0, 0


# Converts scan number to retention time using the mzml file
def set_retention_times(file: str):
    retention_scan_dictionary = {}
    with mzml.read(file) as f:
        for scan in f:
            if scan["ms level"] == 2:
                scan_time = float(scan["scanList"]["scan"][0]["scan start time"])
                scan_time = (scan_time - CON.RETENTION_SHIFT_INTERCEPT) / CON.RETENTION_SHIFT_SLOPE
                scan_time *= CON.MINUTES_TO_SECONDS
                retention_scan_dictionary[scan["index"] + 1] = scan_time
    return retention_scan_dictionary


#################################################################################################
class FullExperiment:

    def __init__(self, time_points: list, differential: bool, free_replications: int, complex_replications: int):
        self.runs = {}
        self._is_differential = differential
        self._num_free_replications = free_replications
        self._num_complex_replications = complex_replications
        self._time_points = time_points
        self._peptides = []
        self._output_files = {}
        self.difference_deviations = {}
        for time in self._time_points:
            self.runs[time] = {True: {}, False: {}}
            self._output_files[time] = {True: {}, False: {}}
            self.difference_deviations[time] = []
        self._file_names = []
        self.deviations_by_time = {}
        self.fractional_deviations_by_time = {}
        self.averages = {}
        for time in self._time_points:
            self.deviations_by_time[time] = []
            self.fractional_deviations_by_time[time] = []
            self.averages[time] = []
        self.protein = parse_protein(CON.PROTEIN_SEQUENCE_FILE)

    def get_num_replications(self, is_complex: bool):
        if is_complex:
            return self._num_complex_replications
        else:
            return self._num_free_replications

    def add_runs(self, time):
        count = 0
        for replication in range(self._num_free_replications):
            self.add_run(time, False, replication, count)
            count += 1
        if self._is_differential:
            for replication in range(self._num_complex_replications):
                self.add_run(time, True, replication, count)
                count += 1

    def read_runs(self):
        for time in self._time_points:
            is_complex = False
            for replication in range(self._num_free_replications):
                self.read_run(time, is_complex, replication)
            if self._is_differential:
                is_complex = True
                for replication in range(self._num_complex_replications):
                    self.read_run(time, is_complex, replication)

    def add_file_names(self):
        for time in self._time_points:
            is_complex = False
            for replication in range(self._num_free_replications):
                self.add_file(time, is_complex, replication)
            if self._is_differential:
                is_complex = True
                for replication in range(self._num_complex_replications):
                    self.add_file(time, is_complex, replication)

    def get_file(self, index: int):
        return self._file_names[index]

    def add_deviation(self, time, deviation, pep):
        if pep.get_mass_shift() > 0:
            self.deviations_by_time[time].append(deviation)
            self.fractional_deviations_by_time[time].append(deviation / pep.get_max_deuterium())

    def add_average(self, time, average):
        self.averages[time].append(average)

    # adds the name of each file to the list
    def add_file(self, time, is_complex: bool, replication):
        complexity = CON.CONDITION1
        if is_complex:
            complexity = CON.CONDITION2
        print("Enter path to mzML for   Time:", time, "  Replication:",
              replication + 1, "  Condition:", complexity)
        file = ""
        while not path.exists(file) or not check_extension(file, ".mzml"):
            file = get_path_input()
            if not check_extension(file, ".mzml"):
                print("Must be a '.mzML' file")
            if not path.exists(file):
                print("Path does not exist.")
        if file in self._file_names:
            warnings.warn("This file has already been added, if this is an error please restart the program.")
        self._file_names.append(file)

    # Adds a new experimental run, adds the peptides on the first run
    def add_run(self, time, complexity: bool, replication, index: int):
        run = ExperimentalRun(run=replication, complexity=complexity, time=time)
        file = self._file_names[index]
        run.read_mzml(file)
        run.hydrogen_deuterium_exchange(CON.IDENTIFICATION_MZML_FILE)
        uptakes = []
        self._peptides = []
        for pep in run.get_peptides():
            self._peptides.append(pep)
            uptakes.append(pep.get_mass_shift())
        self.runs[time][complexity][replication] = uptakes

    def read_run(self, time, complexity: bool, replication):
        run = ReadRun(time, complexity, replication)
        file = generate_output_name(time, complexity, replication)
        uptakes = []
        run.read_output(file)
        self._peptides = []
        for pep in run.get_peptides():
            self._peptides.append(pep)
            uptakes.append(pep.get_mass_shift())
        self.runs[time][complexity][replication] = uptakes

    def get_uptake(self, time, complexity, replication):
        return self.runs[time][complexity][replication]

    def update_headers(self, top_header, bottom_header, header_name):
        top_header.append(header_name)
        for time in self._time_points:
            bottom_header.append(str(time) + " s")
            top_header.append("")
        top_header.pop()
    # Contains multiple subroutines that generate each type of output

    def average_replications(self, averages, deviations, pep, index, complexity):
        for time in self._time_points:
            replications = []
            num_replications = self.get_num_replications(complexity)
            for replication in range(num_replications):
                uptake = self.runs[time][complexity][replication][index]
                if uptake != -1:
                    replications.append(uptake)
            if len(replications) != 0:
                average = sum(replications) / len(replications)
            else:
                average = -1
            averages.append(average)
            if len(replications) > 1:
                deviation = np.std(replications, ddof=1)
            else:
                deviation = -1
            self.add_deviation(time, deviation, pep)
            deviations.append(deviation)

    def generate_recommendation_table_1(self):
        top_header = ["Start", "End", "Sequence", "Peptide mass (Da)", "Retention time (min)"]
        bot_header = ["Start", "End", "Sequence", "Peptide mass (Da)", "Retention time (min)"]
        self.update_headers(top_header, bot_header, "Uptake " + CON.CONDITION1 + " (D)")
        if self._is_differential:
            self.update_headers(top_header, bot_header, "Uptake " + CON.CONDITION2 + " (D)")
        self.update_headers(top_header, bot_header, "Uptake error (SD) - " + CON.CONDITION1 + " (D)")
        if self._is_differential:
            self.update_headers(top_header, bot_header, "Uptake error (SD) - " + CON.CONDITION2 + " (D)")
        peptide_lines = []
        # Each loop generates one line
        for index, pep in enumerate(self._peptides):
            peptide = [pep.get_start(), pep.get_end(), pep.get_sequence(), pep.get_average_mass()]
            start, end = pep.get_rt_start_end()
            rt = (start + end) / 2
            peptide.append(rt)
            averages = []
            deviations = []
            self.average_replications(averages, deviations, pep, index, False)
            if self._is_differential:
                self.average_replications(averages, deviations, pep, index, True)
            peptide.extend(averages)
            peptide.extend(deviations)
            peptide_lines.append(peptide)

        with open(CON.RECOMMENDATION_TABLE_1 + ".csv", "w+", newline='') as f:
            csv_writer = csv.writer(f)
            csv_writer.writerow(top_header)
            csv_writer.writerow(bot_header)
            for line in peptide_lines:
                csv_writer.writerow(line)

    def generate_rows(self, df, time, complexity):
        for index, pep in enumerate(self._peptides):
            row = {}
            if complexity:
                row["Protein state"] = CON.CONDITION2
            else:
                row["Protein state"] = CON.CONDITION1
            row["Start"], row["End"] = pep.get_start(), pep.get_end()
            row["Sequence"] = pep.get_sequence()
            row["Peptide mass (Da)"] = pep.get_monoisotopic_mass()
            start, end = pep.get_rt_start_end()
            average_rt = (start + end) / 2
            row["Retention time (min)"] = average_rt
            row["HDX time (s)"] = time
            row["Uptake (D)"] = pep.get_mass_shift()
            uptakes = []
            for replication in range(self.get_num_replications(complexity)):
                uptakes.append(self.get_uptake(time, complexity, replication)[index])
            stdev = np.std(uptakes, ddof=1)
            row["Uptake SD (D)"] = stdev
            df = df.append(row, ignore_index=True)
        return df

    def generate_recommendation_table_2(self):
        header = ["Protein state", "Sequence", "Start", "End", "Peptide mass (Da)", "Retention time (min)",
                  "HDX time (s)", "Uptake (D)", "Uptake SD (D)"]
        df = pd.DataFrame(columns=header)
        for time in self._time_points:
            df = self.generate_rows(df, time, False)
            if self._is_differential:
                df = self.generate_rows(df, time, True)
        df.to_csv(CON.RECOMMENDATION_TABLE_2 + ".csv", index=False)

    def generate_summary_table(self):
        labels = ["HDX reaction details",
                  "HDX time course",
                  "HDX control samples",
                  "Back-exchange (mean / IQR)",
                  "# of Peptides",
                  "Sequence coverage",
                  "Average peptide length / Redundancy",
                  "Replicates (biological or technical",
                  "Repeatability",
                  "Significant differences in HDX (delta HDX > X D)"]
        time_strings = []
        for time in self._time_points:
            time_strings.append(str(time))
        info = ["",
                ",".join(time_strings),
                "",
                "",
                len(self._peptides),
                "",
                "",
                "{}: {} {}: {}".format(CON.CONDITION1, self.get_num_replications(True),
                                       CON.CONDITION2, self.get_num_replications(False)),
                "",
                ""]
        data = {"Data Set": labels, CON.CONDITION1: info, CON.CONDITION2: info}
        df = pd.DataFrame(data=data)
        df.to_csv(CON.SUMMARY_TABLE + ".csv", index=False)

    def generate_output(self):
        self.generate_recommendation_table_1()
        self.generate_recommendation_table_2()
        if self._is_differential:
            self.generate_differential_woods_plot(CON.WOODS_PLOT_TITLE)
            self.generate_differential_woods_plot(CON.WOODS_PLOT_TITLE, False)
            self.generate_summary_table()

    # confidence is between 0 and 1, df is
    def calculate_confidence_limit(self, time, fractional, confidence):
        n = self._num_complex_replications + self._num_free_replications
        df = n - 2
        if fractional:
            deviations = self.fractional_deviations_by_time[time]
        else:
            deviations = self.deviations_by_time[time]
        if not deviations:
            raise ValueError("add deviations before calculating confidence limit.")
        differences = self.difference_deviations[time]
        stdev = np.std(differences, ddof=1)
        alpha = 1 - confidence
        critical_value = stats.t.ppf(1 - (alpha / 2), df)
        standard_error = stdev / n ** 0.5
        return standard_error * critical_value

    # Saves a plot of name "file"_"time"s.png with a given title. can be fractional or absolute
    # Also generates a table of the values used in the plot
    def generate_differential_woods_plot(self, title: str, is_fractional=True):
        # Formatting
        plt.figure(figsize=(CON.WOODS_PLOT_WIDTH, CON.WOODS_PLOT_HEIGHT))
        plt.title(title)
        plt.xlabel("Sequence")
        plt.tight_layout()
        plt.xlim(0, len(self.protein))
        gray = '#D3D3D3'
        if is_fractional:
            plt.ylabel("Relative Fractional Uptake")
        else:
            plt.ylabel("Relative Uptake (Da)")
        # Generates a plot for each time point
        for time in self._time_points:
            df = pd.read_csv(CON.RECOMMENDATION_TABLE_1 + ".csv", header=[0, 1])
            time_col = str(time) + " s"
            # Plots each peptide
            for _, row in df.iterrows():
                free_deviation = row["Uptake error (SD) - " + CON.CONDITION1 + " (D)"][time_col]
                complex_deviation = row["Uptake error (SD) - " + CON.CONDITION2 + " (D)"][time_col]
                difference = (free_deviation ** 2 + complex_deviation ** 2) ** 0.5
                length = max_deuterium_mlf(row["Sequence"]["Sequence"])
                if is_fractional:
                    difference /= length
                self.difference_deviations[time].append(difference)
            high_confidence = self.calculate_confidence_limit(time, is_fractional, CON.WOODS_PLOT_CONFIDENCE_HIGH)
            low_confidence = self.calculate_confidence_limit(time, is_fractional, CON.WOODS_PLOT_CONFIDENCE_LOW)
            output_table = {
                'Sequence': [],
                'Start': [],
                'End': [],
                'Relative Uptake': [],
                'Relative Fractional Uptake': [],
                'Significant': []
            }
            # Adds each sequence to table
            for _, row in df.iterrows():
                sequence = row['Sequence']['Sequence']
                start, end = row["Start"]["Start"], row["End"]["End"]
                difference = (row["Uptake " + CON.CONDITION2 + " (D)"][time_col] -
                              row["Uptake " + CON.CONDITION1 + " (D)"][time_col])
                absolute_difference = difference
                if is_fractional:
                    difference /= max_deuterium_mlf(sequence)
                if abs(difference) > high_confidence:
                    line_color = 'r'
                else:
                    line_color = gray
                x = start, end
                y = (difference, difference)
                plt.plot(x, y, line_color)
                output_table['Sequence'].append(sequence)
                output_table['Start'].append(start)
                output_table['End'].append(end)
                output_table['Relative Uptake'].append(absolute_difference)
                output_table['Relative Fractional Uptake'].append(difference)
                if abs(difference) > high_confidence:
                    output_table['Significant'].append("Yes")
                else:
                    output_table['Significant'].append("No")
            # Plots the significance lines
            plt.plot((0, len(self.protein)), (0, 0), 'k:')
            plt.plot((0, len(self.protein)), (high_confidence, high_confidence), 'k')
            plt.plot((0, len(self.protein)), (-high_confidence, -high_confidence), 'k')
            plt.plot((0, len(self.protein)), (low_confidence, low_confidence), '--', color='k')
            plt.plot((0, len(self.protein)), (-low_confidence, -low_confidence),  '--', color='k')
            # Generates Output
            data_frame = pd.DataFrame(data=output_table)
            table_file_name = CON.WOODS_TABLE_NAME
            if is_fractional:
                table_file_name += "_fractional"
            data_frame.to_csv(table_file_name + ".csv", index=False)
            plot_file_name = CON.WOODS_PLOT_NAME + "_" + str(time) + "s"
            if is_fractional:
                plot_file_name += "_fractional"
            plot_file_name += ".png"
            plt.savefig(plot_file_name)


##########################################################################
class ReadRun:

    def __init__(self, time, complexity, replication):
        self.peptides = []
        self._replication = replication
        self._complexity = complexity
        self._time = time

    def get_peptides(self):
        return self.peptides

    def read_output(self, output_file):
        df = pd.read_csv(output_file, sep=',')
        pep = Peptide("NA", 0, 0, 0)
        is_first = True
        for index, row in df.iterrows():
            if int(row['Deuterium']) == 0:
                if is_first:
                    is_first = False
                else:
                    pep.set_fit()
                    pep.set_weighted_mass()
                    self.peptides.append(pep)
                pep = Peptide(row['Sequence'], row['SequenceMz'], row['Charge'], "N/A")
                pep.set_sequence_index(row['Start'], row['End'])
                pep.set_retention_times(row['RT'], row['RT'])
            pep.set_deuterium(row['Deuterium'], row['Mz'], row['Intensity'], row['PpmError'])


##########################################################################
class ExperimentalRun:

    def __init__(self, run, complexity, time):
        self.all_peaks = []
        self.peptides = []
        self.windows = {}
        self._replication = run
        self._complexity = complexity
        self._time = time

    # Getters
    def get_peptides(self):
        return self.peptides

    def get_run(self):
        return self._replication

    def get_complexity(self):
        return self._complexity

    def get_time(self):
        return self._time

    def get_tuple_dictionary(self):
        return self.windows

    def get_retention_time(self, index: int):
        return self.all_peaks[index]["retention time"]

    # Setters
    def set_tuple_dictionary(self, tuple_dict: dict):
        self.windows = tuple_dict

    # converts peptide scan number to retention times
    def set_pep_retention_times(self, file: str):
        conversion_dictionary = set_retention_times(file)
        for pep in self.peptides:
            scan = pep.get_scan()
            rt = conversion_dictionary[scan]
            pep.set_retention_times(rt - CON.RETENTION_TOLERANCE, rt + CON.RETENTION_TOLERANCE)

    # adds peptides to peptide_list
    def read_input(self, file: str):
        with open(file, 'r') as f:
            csv_reader = csv.DictReader(f)
            for row in csv_reader:
                self.add_peptide(row['Peptide'], float(row['Precursor']),
                                 int(row['Charge']), float(row["ScanNum"]))

    def process_scan(self, scan):
        retention_time = scan["scanList"]["scan"][0]["scan start time"] * CON.MINUTES_TO_SECONDS
        tuple_list = []
        for j in range(len(scan["m/z array"])):
            if scan["intensity array"][j] > CON.NOISE_LIMIT:
                tuple_list.append((scan["m/z array"][j], scan["intensity array"][j]))
        self.all_peaks.append({"retention time": retention_time, "tuple list": tuple_list})

    # creates a dictionary for each scan with the RT and a
    # list with tuples containing the m/z and intensity
    def read_mzml(self, file: str):
        total = 0
        print()
        print("(initializing)")
        with mzml.read(file) as f:
            for scan in f:
                if scan["ms level"] == 1:
                    total += 1
        count = 0
        retention_time = None
        with mzml.read(file) as f:
            for scan in f:
                if scan["ms level"] == 1:
                    retention_time = scan["scanList"]["scan"][0]["scan start time"]
                    retention_time *= CON.MINUTES_TO_SECONDS
                    count += 1
                    if count % 200 == 0 or count == 1 or count == total:
                        print(count, "/", total, "scans")
                    self.process_scan(scan)
        self.set_tuple_dictionary(self.sliding_window(retention_time))

    # averages intensity of similar masses within RT ranges
    def sliding_window(self, retention_time: float):
        window_count = int(int((retention_time // CON.SLIDING_WINDOW_SIZE) + 1) *
                           (CON.SLIDING_WINDOW_SIZE / SLIDE_AMOUNT)) - 1
        start_time = datetime.now()
        start = 0
        stop = CON.SLIDING_WINDOW_SIZE
        windows = []
        for i in range(window_count):
            windows.append((start, stop))
            start += SLIDE_AMOUNT
            stop += SLIDE_AMOUNT
        window_dictionary = {}
        for window in windows:
            key = window
            window_dictionary[key] = []
            for dictionary in self.all_peaks:
                if window[0] <= dictionary["retention time"] < window[1]:
                    window_dictionary[key].extend(dictionary["tuple list"])
        ret = ((retention_time // CON.SLIDING_WINDOW_SIZE) + 1) * CON.SLIDING_WINDOW_SIZE
        print("(Sliding Window)")
        print("RT {}s / {}s".format(0, int(ret)))
        for rt, tuple_list in window_dictionary.items():
            window_dictionary[rt] = tuple_combine(tuple_list)
            if int(float(rt[0])) != 0 and int(float(rt[0])) % 100 == 0:
                print("RT {}s / {}s".format(int(rt[0]), int(ret)))
        print("Sliding Window Complete")
        print("Time elapsed:", datetime.now() - start_time, "\n")
        return window_dictionary

    def add_peptide(self, sequence, mz, charge, scan):
        self.peptides.append(Peptide(sequence, mz, charge, scan))

    def iterlists(self, index):
        yield from self.all_peaks[index]["tuple list"]

    # outputs a .csv that has the data for a TIC scatter plot
    def generate_total_ion_chromatogram(self):
        mass_ratio = float(input("Enter desired m/z ratio: "))
        tolerance = float(input("Enter m/z tolerance: "))
        start_time = datetime.now()
        tic_list = []
        for i in self.all_peaks:
            retention_time = i["retention time"]
            for mz, intensity in i["tuple list"]:
                if abs(mass_ratio - mz) <= tolerance:
                    tic_list.append((retention_time, intensity))
        elapsed_time = datetime.now() - start_time
        print("\nTime to generate TIC: {}\n".format(elapsed_time))
        user_input = input("What would you like to name this file? ")
        with open(user_input + ".csv", 'w', newline='') as f:
            csv_writer = csv.writer(f)
            for i in tic_list:
                csv_writer.writerow(i)

    def match_peptides(self, pep):
        start, end = pep.get_rt_start_end()

        charge = pep.get_charge()
        pep_mass_over_charge = pep.get_mass_over_charge()
        pep_mass = pep_mass_over_charge * charge
        tuple_list = []
        # collects windows that match the rt of the sequence to user rt tolerance.
        for rt, tup_list in self.get_tuple_dictionary().items():
            if start <= rt[0] <= end or start <= rt[1] <= end:
                tuple_list.extend(tup_list)
        tuple_list.sort(key=lambda x: x[0])
        # searches for a match for each deuteration
        for det in range(pep.get_max_deuterium() + 1):
            ppm_error, mz, intensity = compare(pep_mass + det * CON.DEUTERIUM_MASS_DIFFERENCE,
                                               charge, tuple_list, tuple_list)
            if ppm_error != 0:
                pep.set_deuterium(det, mz, intensity, ppm_error)
        pep.set_fit()
        pep.set_weighted_mass()

    def hydrogen_deuterium_exchange(self, identification_file: str):
        """This function reads the identification file and converts the
        scan #'s to retention times. Then it goes through each peptide and
        matches each deuterium if possible. This ends by writing to an output
        file which contains information on each potential deuteration. This
        should be run for each experimental run."""
        self.read_input(CON.IDENTIFICATION_CSV_FILE)
        self.set_pep_retention_times(identification_file)
        count = 0
        start_time = datetime.now()
        print("Begin matching")
        for pep in self.peptides:
            self.match_peptides(pep)
            count += 1
            if count % 20 == 0 or count == 1 or count == len(self.peptides):
                print(count, "/", len(self.peptides))
        print("Time to match:", datetime.now() - start_time, '\n')
        file = generate_output_name(self._time, self._complexity, self._replication)
        self.write_hdx(file)

    # outputs tabular info of HDX results to .csv
    def write_hdx(self, file: str):
        output_exists = path.exists(CON.FULL_HDX_OUTPUT + ".csv")

        with open(file, "w+", newline='') as f:
            csv_writer = csv.writer(f)
            if not output_exists:
                header = ["Start", "End", "Sequence", "Charge", "SequenceMz", "Condition",
                          "Deuterium", "RT", "Mz", "Intensity", "PpmError", "Average", "Shift",
                          "Gaussian Fit"]
                csv_writer.writerow(header)
            lines = []
            for pep in self.peptides:
                line_list = pep.get_rows(self._complexity)
                for line in line_list:
                    lines.append(line)
            for line in lines:
                csv_writer.writerow(line)


def generate_output_name(time, is_complex, replication):
    complexity = CON.CONDITION1
    if is_complex:
        complexity = CON.CONDITION2
    file = (CON.FULL_HDX_OUTPUT + "_" + str(time) + "s_"
            + str(complexity) + "_" + str(replication + 1) + ".csv")
    return file


def read_sequence(string):
    has_period = ('.' in string)
    if not has_period:
        return string
    found_period = False
    return_string = ""
    for letter in string:
        if letter == '.' and found_period:
            return return_string
        elif letter == '.':
            found_period = True
        elif found_period:
            return_string += letter
    raise ValueError("Sequence formatted incorrectly")


######################################################################
class Peptide:
    def __init__(self, sequence, mz, charge, scan):
        self._windows = []
        self._sequence = read_sequence(sequence)
        self._charge = charge
        self._mass_over_charge = mz
        self._rt_start = 0
        self._rt_end = float("inf")
        self._scan = scan
        self._weighted_mass_to_charge = 0
        self._deuterium_dictionary = {}
        self._mass_shift = 0
        self.max_deuterium = sequence_to_max_deuterium(self._sequence)
        for det in range(self.max_deuterium + 1):
            self._deuterium_dictionary[det] = {"m/z": 0, "intensity": 0, "ppm": 0}
        self._average_mass = 0
        self.set_average_mass()
        self._protein = parse_protein(CON.PROTEIN_SEQUENCE_FILE)
        self._start, self._end = find_start_end(self._sequence, self._protein)
        self._fit = 0  # Gaussian fit

    def __str__(self):
        self.__repr__()
        return ""

    def __repr__(self):
        print("Sequence:", self._sequence)
        print("Average Mass:", self._average_mass)
        print("Max Deuterium:", self.get_max_deuterium())
        print("Fit:", self._fit)
        print("Mass Shift:", self._mass_shift)

    # Getters
    def get_fit(self):
        return self._fit

    def get_mass_shift(self):
        return self._mass_shift

    def get_weighted_mass(self):
        return self._weighted_mass_to_charge

    def get_average_mass(self):
        return self._average_mass

    def get_start(self):
        return self._start

    def get_end(self):
        return self._end

    def get_sequence(self):
        return self._sequence

    def get_charge(self):
        return self._charge

    def get_mass_over_charge(self):
        return self._mass_over_charge

    def get_max_deuterium(self):
        return self.max_deuterium

    def get_rt_start_end(self):
        return self._rt_start, self._rt_end

    # returns peptide data for the detailed output
    def get_rows(self, complexity: str):
        if complexity:
            condition = CON.CONDITION2
        else:
            condition = CON.CONDITION1
        line_list = []
        self.set_fit()
        for i in range(self.get_max_deuterium() + 1):

            mz, intensity, ppm = self.get_deuterium(i)
            line = ["", "", "", 0, 0, "", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            line[0] = self._start
            line[1] = self._end
            line[2] = self._sequence
            line[3] = self._charge
            line[4] = self._mass_over_charge
            line[5] = condition
            line[6] = i
            line[7] = (self._rt_end + self._rt_start) / 2 / CON.MINUTES_TO_SECONDS
            line[8] = mz
            line[9] = intensity
            line[10] = ppm
            line[11] = self._weighted_mass_to_charge
            line[12] = self._mass_shift
            line[13] = self._fit
            line_list.append(line)
        return line_list

    def get_scan(self):
        return self._scan

    # returns mz, intensity, ppm_error of of one match
    def get_deuterium(self, deuterium):
        mz = self._deuterium_dictionary[deuterium]["m/z"]
        intensity = self._deuterium_dictionary[deuterium]["intensity"]
        ppm = self._deuterium_dictionary[deuterium]["ppm"]
        return mz, intensity, ppm

    def get_monoisotopic_mass(self):
        return self._mass_over_charge * self._charge + CON.MASS_OF_WATER

    # Setters
    def set_sequence_index(self, start, end):
        self._start = start
        self._end = end

    def set_fit(self):
        intensities = []
        for key in self._deuterium_dictionary.keys():
            intensities.append(self._deuterium_dictionary[key]["intensity"])
        self._fit = fit_gaussian(intensities)

    def set_windows(self, window):
        self._windows.append(window)

    def set_deuterium(self, det, mz, intensity, ppm_error):
        self._deuterium_dictionary[det]["m/z"] = mz
        self._deuterium_dictionary[det]["intensity"] = intensity
        self._deuterium_dictionary[det]["ppm"] = ppm_error

    def set_retention_times(self, start: float, end: float):
        self._rt_start = start
        self._rt_end = end

    def set_weighted_mass(self):
        total_intensity = 0
        mass = 0
        for det in self._deuterium_dictionary.keys():
            mz, intensity, ppm_error = self.get_deuterium(det)
            total_intensity += intensity
            mass += intensity * mz
        if total_intensity == 0:
            mass = 0
        else:
            mass /= total_intensity
        self._weighted_mass_to_charge = mass
        weighted_mass = self._weighted_mass_to_charge - CON.MASS_OF_HYDROGEN
        self._mass_shift = weighted_mass * self._charge - self._average_mass
        if total_intensity == 0 or self._fit == "Insufficient Match":
            self._weighted_mass_to_charge = -1
            self._mass_shift = -1

    # calculates the average mass from the sequence (Uses values in the PARAMETERS.py file)
    def set_average_mass(self):
        mass = 0
        for amino in self._sequence:
            mass += CON.PEPTIDE_MASS_DICTIONARY[amino]
        mass += CON.MASS_OF_WATER
        self._average_mass = mass


def get_time_points():
    num_time_points = 0
    while num_time_points < 1:
        try:
            num_time_points = int(input("How many time points? "))
        except ValueError:
            pass
        if num_time_points < 1:
            print("Please enter a positive integer.")
    print("Please enter each time point")
    time_points = []
    for i in range(num_time_points):
        usr_input = -1
        while usr_input < 0:
            try:
                usr_input = int(input("Time point #{} (s): ".format(i + 1)))
            except ValueError:
                pass
            if usr_input < 0:
                print("Please enter 0 or a positive integer.")
            if usr_input in time_points:
                print("That time point has already been entered.")
                usr_input = -1
        time_points.append(usr_input)
    return time_points


def get_is_differential():
    is_differential = -1
    while is_differential == -1:
        differential_input = input("Was this a differential experiment (Y/N)? ").lower().strip()
        if differential_input[0] == 'y':
            is_differential = True
        elif differential_input[0] == 'n':
            is_differential = False
        else:
            print("Please enter 'Yes' or 'No'")
    return is_differential


def get_num_replications_input(is_complex: bool):
    num_replications = 0
    if is_complex:
        condition = CON.CONDITION2
    else:
        condition = CON.CONDITION1
    while num_replications < 1:
        try:
            num_replications = int(input("How many {} replications? ".format(condition)))
        except ValueError:
            pass
        if num_replications < 1:
            print("Please enter a positive integer.")
    return num_replications


def show_menu():
    print("Please enter the number of one of the following selections:")
    print("(1) Detailed output from mzML (Only used once per experiment)")
    print("(2) Final summary and figures from detailed outputs")
    print("(3) Quit")


##############################################################################
def main():
    try:
        check_parameters()
        # Get User Input
        time_points = get_time_points()
        is_differential = get_is_differential()
        num_free_replications = get_num_replications_input(False)
        if is_differential:
            num_complex_replications = get_num_replications_input(True)
        else:
            num_complex_replications = 0
        print()
        menu_input = None
        while menu_input != 'q':
            show_menu()
            try:
                menu_input = input().strip()[0]
            except IndexError:
                pass
            if menu_input == '1':
                start_time = datetime.now()
                # Perform peak matching and generate output
                experiment = FullExperiment(time_points, is_differential, num_free_replications, num_complex_replications)
                experiment.add_file_names()
                for time in time_points:
                    experiment.add_runs(time)
                print("Total Time Elapsed:", datetime.now() - start_time)
            if menu_input == '2':
                print("Generating Output files")
                experiment = FullExperiment(time_points, is_differential, num_free_replications, num_complex_replications)
                experiment.read_runs()
                experiment.generate_output()
                print("\nSuccess!\n")
            if menu_input == '3':
                quit()
    except AttributeError:
        print("ERROR: The parameter file seems to be damaged. Please re-download from the source. "
              "Ensure you do not alter the names of variables")


if __name__ == '__main__':
    main()
