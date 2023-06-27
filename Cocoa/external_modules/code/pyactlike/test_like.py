import pytest
import numpy as np
import pyactlike


def get_example_spectra():
    like = pyactlike.ACTPowerSpectrumData()
    filename = like.data_dir + "bf_ACTPol_WMAP_lcdm.minimum.theory_cl"
    tt_lmax = 5000
    ell, dell_tt, dell_te, dell_ee = np.genfromtxt(
        filename,
        delimiter=None,
        unpack=True,
        max_rows=tt_lmax - 1,
        usecols=(0, 1, 2, 3),
    )
    return ell, dell_tt, dell_te, dell_ee


def test_TTTEEE():
    """This function tests out the basic functionality of this likelihood code."""

    ell, dell_tt, dell_te, dell_ee = get_example_spectra()
    like = pyactlike.ACTPowerSpectrumData()
    chi2 = -2 * like.loglike(dell_tt, dell_te, dell_ee, 1.003)
    print("ACTPol chi2 = " + "{0:.12f}".format(chi2))
    print("Expected:     288.252869629064")
    assert np.isclose(chi2, 288.252869629064)


def test_bmin():
    # nonzero bmin
    ell, dell_tt, dell_te, dell_ee = get_example_spectra()
    like = pyactlike.ACTPowerSpectrumData(bmin=24)
    chi2 = -2 * like.loglike(dell_tt, dell_te, dell_ee, 1.003)
    print("ACTPol chi2 = " + "{0:.12f}".format(chi2))
    print("Expected:     235.146031846935")
    assert np.isclose(chi2, 235.146031846935)


def test_single_channel():
    """This function tests out the single channels functionality of this likelihood code."""

    # TT only
    like = pyactlike.ACTPowerSpectrumData(use_tt=True, use_te=False, use_ee=False)
    ell, dell_tt, dell_te, dell_ee = get_example_spectra()
    chi2 = -2 * like.loglike(dell_tt, dell_te, dell_ee, 1.003)
    assert np.isclose(chi2, 97.4331220842641)

    # TE only
    like = pyactlike.ACTPowerSpectrumData(use_tt=False, use_te=True, use_ee=False)
    chi2 = -2 * like.loglike(dell_tt, dell_te, dell_ee, 1.003)
    assert np.isclose(chi2, 81.6194890026420)

    # EE only
    like = pyactlike.ACTPowerSpectrumData(use_tt=False, use_te=False, use_ee=True)
    chi2 = -2 * like.loglike(dell_tt, dell_te, dell_ee, 1.003)
    assert np.isclose(chi2, 98.5427508626497)


@pytest.mark.skip(reason="cobaya optional")
def test_cobaya():
    """Test the Cobaya interface to the ACT likelihood."""
    from cobaya.yaml import yaml_load
    from cobaya.model import get_model

    info_yaml = r"""
        likelihood:
            pyactlike.ACTPol_lite_DR4:
                components: 
                    - tt
                    - te
                    - ee
                lmax: 6000

        theory:
            camb:
                extra_args:
                    lens_potential_accuracy: 1

        params:
            ns:
                prior:
                  min: 0.8
                  max: 1.2
            H0:
                prior:
                  min: 40
                  max: 100       
            yp2:
                prior:
                    min: 0.5
                    max: 1.5       
        """
    info = yaml_load(info_yaml)
    model = get_model(info)
    assert np.isfinite(model.loglike({"ns": 1.0, "H0": 70, "yp2": 1.0})[0])
