from libra_toolbox.neutron_detection.diamond.process_data import *
import pytest


def test_get_avg_neutron_rate():

    time_values = np.arange(1, 10)
    print(time_values)
    t_min = 2.5
    t_max = 9.0

    avg_neutron_rate = get_avg_neutron_rate(time_values, t_min, t_max)

    expected_number_of_counts = 6  # there are 6 counts between 2.5 and 9.0

    expected_avg_neutron_rate = {
        "counts": expected_number_of_counts,
        "err": np.sqrt(expected_number_of_counts),
        "count rate": expected_number_of_counts / (t_max - t_min),
        "count rate err": np.sqrt(expected_number_of_counts) / (t_max - t_min),
    }

    assert avg_neutron_rate == expected_avg_neutron_rate


@pytest.mark.parametrize("delimiter", [",", ";", "\t"])
@pytest.mark.parametrize("extention", ["csv", "CSV"])
def test_get_time_energy_values(tmpdir, delimiter, extention):

    # make data
    time_values_out = np.random.rand(10)
    energy_values_out = np.random.rand(10)
    extra_column = np.random.rand(10)

    # dump 2 columns to the same csv file
    filename = tmpdir.join(f"test.{extention}")
    np.savetxt(
        filename,
        np.column_stack((time_values_out, energy_values_out, extra_column)),
        delimiter=delimiter,
    )

    # read the data back
    time_values_in, energy_values_in = get_time_energy_values(
        tmpdir, time_column=0, energy_column=1, delimiter=delimiter
    )

    assert np.allclose(time_values_in, time_values_out)
    assert np.allclose(energy_values_in, energy_values_out)


def test_main(tmpdir):
    total_time_s = 100  # s
    total_time_ps = total_time_s * 1e12  # s to ps

    # peak 1
    nb_counts_peak1 = int(2e4)
    size_peak1 = nb_counts_peak1
    time_values_out_peak1 = np.random.rand(size_peak1) * total_time_ps
    mean_energy_peak1 = 4e6
    std_energy_peak1 = 0.1e6
    energy_values_peak1 = np.random.normal(
        mean_energy_peak1, std_energy_peak1, size_peak1
    )

    # make data
    nb_counts_peak2 = int(7e4)
    size_peak2 = nb_counts_peak2
    time_values_out_peak2 = np.random.rand(size_peak2) * total_time_ps
    mean_energy_peak2 = 14e6
    std_energy_peak2 = 1e6

    energy_values_peak2 = np.random.normal(
        mean_energy_peak2, std_energy_peak2, size_peak2
    )

    # import matplotlib.pyplot as plt

    # plt.hist(energy_values_peak1)
    # plt.hist(energy_values_peak2)
    # plt.show()

    filename1 = tmpdir.join(f"peak1.csv")
    np.savetxt(
        filename1,
        np.column_stack((time_values_out_peak1, energy_values_peak1)),
        delimiter=",",
    )
    filename2 = tmpdir.join(f"peak2.csv")
    np.savetxt(
        filename2,
        np.column_stack((time_values_out_peak2, energy_values_peak2)),
        delimiter=",",
    )

    # run
    bin_time = 20  # s
    res1 = main(
        tmpdir,
        bin_time=bin_time,
        energy_peak_min=mean_energy_peak1 - std_energy_peak1 * 2,
        energy_peak_max=mean_energy_peak1 + std_energy_peak1 * 2,
        delimiter=",",
        time_column=0,
        energy_column=1,
    )

    res2 = main(
        tmpdir,
        bin_time=bin_time,
        energy_peak_min=mean_energy_peak2 - std_energy_peak2 * 2,
        energy_peak_max=mean_energy_peak2 + std_energy_peak2 * 2,
        delimiter=",",
        time_column=0,
        energy_column=1,
    )

    # test
    expected_count_rate_total = (nb_counts_peak1 + nb_counts_peak2) / total_time_s
    expected_count_rate_peak1 = nb_counts_peak1 / total_time_s
    expected_count_rate_peak2 = nb_counts_peak2 / total_time_s

    assert "all_count_rates" in res1
    assert "all_count_rate_bins" in res1
    assert "time_values" in res1
    assert "energy_values" in res1
    assert "peak_count_rates" in res1
    assert "peak_count_rate_bins" in res1
    assert "peak_time_values" in res1

    # check that the count rates are as expected
    assert np.allclose(res1["all_count_rates"], expected_count_rate_total, rtol=0.1)
    assert np.allclose(res1["peak_count_rates"], expected_count_rate_peak1, rtol=0.1)

    assert np.allclose(res2["all_count_rates"], expected_count_rate_total, rtol=0.1)
    assert np.allclose(res2["peak_count_rates"], expected_count_rate_peak2, rtol=0.1)

    assert len(res1["all_count_rate_bins"]) == total_time_s / bin_time
    assert len(res2["all_count_rate_bins"]) == total_time_s / bin_time
