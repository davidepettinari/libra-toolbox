import pandas as pd


class LSCFileReader:
    def __init__(self, file_path, vial_labels=None):
        self.file_path = file_path
        self.vial_labels = vial_labels
        self.data = None
        self.header_content = None

    def read_file(self):
        # first read the file without dataframe to find the line starting with S#
        header_lines = []
        with open(self.file_path, "r") as file:
            lines = file.readlines()
            for i, line in enumerate(lines):
                if line.startswith("S#"):
                    start = i
                    break
                header_lines.append(line)
        self.header_content = "".join(header_lines)

        # read the file with dataframe starting from the line with S#
        self.data = pd.read_csv(self.file_path, skiprows=start)

    def get_bq1_values(self):
        return self.data["Bq:1"].tolist()

    def get_bq1_values_with_labels(self):
        if self.vial_labels is None:
            raise ValueError("Vial labels must be provided")

        assert len(self.vial_labels) == len(
            self.get_bq1_values()
        ), "Vial labels and Bq:1 values are not equal in length, remember to give None as a label for missing vials"

        values = self.get_bq1_values()
        labelled_values = {label: val for label, val in zip(self.vial_labels, values)}

        return labelled_values

    def get_count_times(self):
        return self.data["Count Time"].tolist()

    def get_lum(self):
        return self.data["LUM"].tolist()
