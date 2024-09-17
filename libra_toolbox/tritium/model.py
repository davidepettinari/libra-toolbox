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
        radius: pint.Quantity = None,
        height: pint.Quantity = None,
        TBR: pint.Quantity = None,
    ) -> None:
        self.radius = radius
        self.height = height

        self.L_wall = 0.06 * ureg.inches

        self.neutron_rate = 3e8 * ureg.neutron * ureg.s**-1

        self.TBR = TBR
        self.irradiations = []

        self.k_wall = 1.9e-8 * ureg.m * ureg.s**-1  # from Kumagai
        self.k_top = 4.9e-7 * ureg.m * ureg.s**-1  # from Kumagai

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

    def run(self, t_final):
        concentration_units = ureg.particle * ureg.m**-3
        time_units = ureg.s
        initial_concentration = 0
        res = solve_ivp(
            fun=self.rhs,
            t_span=(0, t_final.to(time_units).magnitude),
            y0=[initial_concentration],
            t_eval=np.sort(
                np.concatenate(
                    [
                        np.linspace(0, t_final.to(time_units).magnitude, 10000),
                        [irr[1].to(time_units).magnitude for irr in self.irradiations],
                    ]
                )
            ),
            # method="RK45",  # RK45 doesn't catch the end of irradiations properly... unless constraining the max_step
            # max_step=(0.5 * ureg.h).to(time_units).magnitude,
            # method="Radau",
            method="BDF",
        )
        self.times = res.t * time_units
        self.concentrations = res.y[0] * concentration_units

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
