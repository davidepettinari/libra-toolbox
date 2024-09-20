import pint
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.integrate import solve_ivp

ureg = pint.UnitRegistry()
ureg.setup_matplotlib()
ureg.define("neutron = 1 * particle = n")

SPECIFIC_ACT = 3.57e14 * ureg.Bq * ureg.g**-1
MOLAR_MASS = 6.032 / 2 * ureg.g * ureg.mol**-1


class Model:
    def __init__(
        self,
        radius: pint.Quantity,
        height: pint.Quantity,
        TBR: pint.Quantity,
        neutron_rate: pint.Quantity,
        k_wall: pint.Quantity,
        k_top: pint.Quantity,
        irradiations: list,
    ) -> None:
        self.radius = radius
        self.height = height

        self.L_wall = 0.06 * ureg.inches

        self.neutron_rate = neutron_rate

        self.TBR = TBR
        self.irradiations = irradiations

        self.k_wall = k_wall
        self.k_top = k_top

        self.concentrations = []
        self.times = []

    @property
    def volume(self):
        return self.A_top * self.height

    @property
    def A_top(self):
        return np.pi * self.radius**2

    @property
    def A_wall(self):
        perimeter_wall = 2 * np.pi * (self.radius + self.L_wall)
        return perimeter_wall * self.height + self.A_top

    def source(self, t):
        for irradiation in self.irradiations:
            irradiation_start = irradiation[0]
            irradiation_stop = irradiation[1]
            if irradiation_start < t < irradiation_stop:
                return self.TBR * self.neutron_rate
        return 0 * self.TBR * self.neutron_rate

    def Q_wall(self, c_salt):
        return self.A_wall * self.k_wall * c_salt

    def Q_top(self, c_salt):
        return self.A_top * self.k_top * c_salt

    def rhs(self, t, c):
        t *= ureg.s
        c *= ureg.particle * ureg.m**-3

        return self.volume.to(ureg.m**3) ** -1 * (
            self.source(t).to(ureg.particle * ureg.s**-1)
            - self.Q_wall(c).to(ureg.particle * ureg.s**-1)
            - self.Q_top(c).to(ureg.particle * ureg.s**-1)
        )

    def _generate_time_intervals(self, t_final):
        time_intervals = []
        previous_tf = None

        for irr in self.irradiations:
            t0 = irr[0]
            tf = irr[1]

            if previous_tf is not None:
                time_intervals.append((previous_tf, t0))
            time_intervals.append((t0, tf))
            previous_tf = tf

        # Add the final interval from the last tf to t_final
        if previous_tf is not None and previous_tf < t_final:
            time_intervals.append((previous_tf, t_final))

        return time_intervals

    def run(self, t_final):
        concentration_units = ureg.particle * ureg.m**-3
        time_units = ureg.s
        initial_concentration = 0
        time_intervals = self._generate_time_intervals(t_final)

        for interval in time_intervals:
            t0 = interval[0].to(time_units).magnitude
            tf = interval[1].to(time_units).magnitude

            res = solve_ivp(
                fun=self.rhs,
                t_span=(t0, tf),
                y0=[initial_concentration],
                t_eval=np.linspace(t0, tf, 1000),
                # method="RK45",  # RK45 doesn't catch the end of irradiations properly... unless constraining the max_step
                # max_step=(0.5 * ureg.h).to(time_units).magnitude,
                # method="Radau",
                method="BDF",
            )
            self.times.append(res.t)
            self.concentrations.append(res.y[0])
            initial_concentration = res.y[0][-1]

        self.times = np.concatenate(self.times) * time_units
        self.concentrations = np.concatenate(self.concentrations) * concentration_units

    def reset(self):
        self.concentrations = []
        self.times = []

    def integrated_release_top(self):
        top_release = self.Q_top(self.concentrations)
        integrated_top = cumulative_trapezoid(
            top_release.to(ureg.particle * ureg.h**-1).magnitude,
            self.times.to(ureg.h).magnitude,
            initial=0,
        )
        integrated_top *= ureg.particle  # attach units
        return integrated_top

    def integrated_release_wall(self):
        wall_release = self.Q_wall(self.concentrations)
        integrated_wall = cumulative_trapezoid(
            wall_release.to(ureg.particle * ureg.h**-1).magnitude,
            self.times.to(ureg.h).magnitude,
            initial=0,
        )
        integrated_wall *= ureg.particle  # attach units
        return integrated_wall


def quantity_to_activity(Q):
    return Q * SPECIFIC_ACT * MOLAR_MASS


def activity_to_quantity(A):
    return A / (SPECIFIC_ACT * MOLAR_MASS)
