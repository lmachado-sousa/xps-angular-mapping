"""Headless tests for EasyPizza's data processing (no Tk or display needed).

Run with: python -m unittest -v
"""
import os
import tempfile
import unittest

import numpy as np

import angular_map as am
import gaussian
import savitzky
from shirley import shirley

HERE = os.path.dirname(os.path.abspath(__file__))
MODES = (am.RAW, am.SYM60, am.SYM120, am.SYM90)


def repo_file(name):
    return os.path.join(HERE, name)


def write_xpd(directory, intensity, thetas=range(9, 46, 3), phis=np.arange(0, 241, 3.0), skip=()):
    """Write a synthetic .xpd file; `intensity(theta, phi)` takes degrees."""
    path = os.path.join(directory, "synthetic.xpd")
    with open(path, "w") as f:
        f.write("Intensity    Theta    Phi\n")
        for t in thetas:
            for p in phis:
                if (t, p) not in skip:
                    f.write("%f    %f    %f\n" % (intensity(t, p), t, p))
    return path


class RepoDataTest(unittest.TestCase):
    def test_xpd_raw_points_follow_file_angles(self):
        data = np.loadtxt(repo_file("1ML1.xpd"), skiprows=1)
        theta, phi, ki = am.load_map(repo_file("1ML1.xpd"))
        self.assertEqual(len(ki), 1943)
        np.testing.assert_array_equal(theta, data[:, 1])           # theta stays in degrees
        np.testing.assert_allclose(phi, np.radians(data[:, 2]))    # phi comes back in radians
        self.assertEqual(sorted(set(theta)), list(range(9, 79, 3)))
        for t in set(theta):                                         # chi = (I - <I>) / <I> per ring
            ring = theta == t
            intensity = data[ring, 0]
            np.testing.assert_allclose(ki[ring], intensity / intensity.mean() - 1, atol=1e-12)

    def test_hol_uses_columns_0_1_4(self):
        data = np.loadtxt(repo_file("Ni111_3p.hol"), skiprows=1)
        theta, phi, ki = am.load_map(repo_file("Ni111_3p.hol"))
        np.testing.assert_array_equal(theta, data[:, 0])
        np.testing.assert_allclose(phi, np.radians(data[:, 1]))    # 0 to 360 degrees, all kept
        np.testing.assert_array_equal(ki, data[:, 4])

    def test_dat_uses_columns_0_1_2_and_starts_phi_at_zero(self):
        data = np.loadtxt(repo_file("1test.dat"), skiprows=1)
        theta, phi, ki = am.load_map(repo_file("1test.dat"))
        np.testing.assert_array_equal(theta, data[:, 0])
        np.testing.assert_allclose(phi, np.radians(data[:, 1] - 47))  # file starts at 47 degrees
        np.testing.assert_array_equal(ki, data[:, 2])

    def test_interpolated_maps_have_expected_grid(self):
        for name in ("1ML1.xpd", "1test.dat", "Ni111_3p.hol"):
            for mode in MODES:
                with self.subTest(file=name, mode=mode):
                    theta, phi, ki = am.load_map(repo_file(name), mode)
                    self.assertTrue(np.all((phi >= 0) & (phi <= 2 * np.pi)))
                    gtheta, gphi, z = am.interpolate_map(theta, phi, ki, mode)
                    for grid in (gtheta, gphi, z):
                        self.assertEqual(grid.shape, (600, 600))
                    self.assertTrue(np.all(np.isfinite(z)))
                    if mode == am.RAW:
                        self.assertAlmostEqual(gtheta.min(), theta.min())
                        self.assertAlmostEqual(gtheta.max(), theta.max())
                        self.assertAlmostEqual(gphi.min(), 0)
                        self.assertAlmostEqual(gphi.max(), phi.max())
                    else:
                        self.assertTrue(np.all((gphi >= 0) & (gphi < 2 * np.pi)))
                        np.testing.assert_array_equal(z[gtheta < theta.min()], ki.min())

    def test_unsupported_extension(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "map.txt")
            with open(path, "w") as f:
                f.write("Intensity Theta Phi\n1 9 0\n")
            with self.assertRaisesRegex(ValueError, "Unsupported file type"):
                am.load_map(path)


class SymmetryTest(unittest.TestCase):
    """Synthetic maps whose symmetrised result is known in advance."""

    def test_60_degree_pattern_folds_onto_itself(self):
        def intensity(t, p):
            return 1000 * (1 + t / 100) * (1 + 0.3 * np.cos(np.radians(6 * p)))

        with tempfile.TemporaryDirectory() as tmp:
            theta, phi, ki = am.load_map(write_xpd(tmp, intensity), am.SYM60)
        phi_deg = np.degrees(phi)
        self.assertAlmostEqual(phi_deg.max(), 360)
        for t in set(theta):
            window = intensity(t, np.arange(0, 61, 3.0))  # the folded sector, 0 to 60 degrees
            expected = intensity(t, phi_deg[theta == t]) / window.mean() - 1
            np.testing.assert_allclose(ki[theta == t], expected, atol=1e-6)

    def test_fold_is_the_mean_of_equivalent_azimuths(self):
        rng = np.random.default_rng(1)
        for mode, period, scan in ((am.SYM60, 60, 240), (am.SYM120, 120, 240), (am.SYM90, 90, 270)):
            with self.subTest(mode=mode):
                phis = np.arange(0, scan + 1, 3.0)
                values = {p: 500 + 1000 * rng.random() for p in phis}
                theta, phi, ki = am._xpd_points(
                    np.array([values[p] for p in phis]), np.full(len(phis), 30.0), phis, mode)
                sector = np.arange(0, period + 1, 3.0)
                folded = np.array([np.mean([values[p + k] for k in range(0, scan - period + 1, period)])
                                   for p in sector])
                np.testing.assert_array_equal(phi, sector)
                np.testing.assert_allclose(ki, folded / folded.mean() - 1, atol=1e-12)

    def test_ring_without_phi_zero_keeps_its_angles(self):
        # 1ML1.xpd has this case: the theta = 21 ring starts at phi = 3 degrees
        with tempfile.TemporaryDirectory() as tmp:
            path = write_xpd(tmp, lambda t, p: 1000 + p, skip={(21, 0.0)})
            for mode, last in ((am.RAW, 240), (am.SYM60, 60), (am.SYM120, 120)):
                with self.subTest(mode=mode):
                    theta, phi, ki = am._xpd_points(*np.loadtxt(path, skiprows=1).T, mode)
                    ring = np.array(phi)[np.array(theta) == 21]
                    np.testing.assert_array_equal(ring, np.arange(3, last + 1, 3.0))

    def test_symmetrised_map_keeps_azimuth_orientation(self):
        def intensity(t, p):  # bumps at 30, 150 and 270 degrees
            return 1000 * (1 + sum(np.exp(-((p - c + 180) % 360 - 180) ** 2 / 200) for c in (30, 150, 270)))

        with tempfile.TemporaryDirectory() as tmp:
            path = write_xpd(tmp, intensity)
            for mode, centres in ((am.RAW, [30, 150]), (am.SYM120, [30, 150, 270])):
                with self.subTest(mode=mode):
                    theta, phi, ki = am.load_map(path, mode)
                    gtheta, gphi, z = am.interpolate_map(theta, phi, ki, mode)
                    ring = np.abs(gtheta - 30) < 0.5
                    for centre in centres:
                        near = ring & (np.abs(np.degrees(gphi) - centre) < 45)
                        peak = np.degrees(gphi[near][np.argmax(z[near])])
                        self.assertAlmostEqual(peak, centre, delta=3)


class SpectrumHelperTest(unittest.TestCase):
    def test_savitzky_golay_keeps_a_cubic(self):
        x = np.linspace(-3, 3, 101)
        y = 0.5 * x ** 3 - 2 * x ** 2 + x + 3
        np.testing.assert_allclose(savitzky.savitzky_golay(y, 7, 3)[3:-3], y[3:-3], atol=1e-9)

    def test_gaussian_fit_recovers_parameters(self):
        x = np.linspace(-3, 3, 121)
        A, mu, sigma = gaussian.fit_iterative(x, gaussian.gauss(x, A=50, mu=0.3, sigma=0.7))
        np.testing.assert_allclose([A, mu, sigma], [50, 0.3, 0.7], rtol=1e-6)

    def test_shirley_removes_a_shirley_background(self):
        x = np.linspace(60, 80, 200)
        peak = 1000 * np.exp(-(x - 70) ** 2 / (2 * 0.8 ** 2))
        right_area = np.cumsum(peak[::-1])[::-1] - peak  # area at higher indices
        background = 200 + 300 * right_area / right_area[0]
        y = peak + background
        np.testing.assert_allclose(shirley(y, 0, len(y) - 1)[1:], peak[1:], atol=5)


if __name__ == "__main__":
    unittest.main()
