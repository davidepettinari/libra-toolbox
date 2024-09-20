import os
import numpy as np


class DataProcessor:
    def __init__(self) -> None:
        self.files = []

        self.time_values = []
        self.energy_values = []

    def add_file(
        self,
        filename: str,
        time_column: int,
        energy_column: int,
        scale_time: bool = True,
        **kwargs,
    ):
        self.files.append(filename)

        # Should we store the data for each file separately too?
        data = np.genfromtxt(filename, **kwargs)

        time_values = data[:, time_column]
        if scale_time:
            # convert times from ps to s
            time_values *= 1 / 1e12
        energy_values = data[:, energy_column]

        # Append time and energy values to the list
        self.time_values = np.concatenate((self.time_values, time_values))
        self.energy_values = np.concatenate((self.energy_values, energy_values))

        # sort time and energy values
        inds = np.argsort(self.time_values)
        self.time_values = np.array(self.time_values)[inds]
        self.energy_values = np.array(self.energy_values)[inds]

        print(f"Added file: {filename} containing {len(time_values)} events")

    def get_count_rate(self, bin_time: float, energy_window: tuple = None):
        time_values = self.time_values.copy()
        energy_values = self.energy_values.copy()

        time_bins = np.arange(0, time_values[-2], bin_time)

        if energy_window is not None:
            peak_mask = (energy_values > energy_window[0]) & (
                energy_values < energy_window[1]
            )
            time_values = time_values[peak_mask]

        count_rates, count_rate_bins = np.histogram(time_values, bins=time_bins)
        count_rates = count_rates / bin_time

        return count_rates, count_rate_bins

    def get_avg_rate(self, t_min: float, t_max: float, energy_window: tuple = None):
        time_values = self.time_values.copy()
        energy_values = self.energy_values.copy()

        if energy_window is not None:
            peak_mask = (energy_values > energy_window[0]) & (
                energy_values < energy_window[1]
            )
            time_values = time_values[peak_mask]

        # Create mask to only count pulses of any energy in section time window
        idx = np.logical_and(self.time_values > t_min, self.time_values < t_max)

        counts = len(self.time_values[idx])
        error = np.sqrt(len(self.time_values[idx]))
        delta_t = t_max - t_min
        count_rate = counts / delta_t
        count_rate_err = error / delta_t

        return count_rate, count_rate_err


if __name__ == "__main__":
    pass
