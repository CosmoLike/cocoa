# Foregrounds for SPT
import logging
import os

import numpy as np
from cobaya.log import HasLogger, LoggedError


def cosmo_scale_ksz(H0, sigma8, omegb, omegam, ns, tau):
    return (
        (H0 / 71.0) ** 1.7
        * (sigma8 / 0.8) ** 4.7
        * (omegab / 0.044) ** 2.1
        * (omegam / 0.264) ** (-0.44)
        * (ns / 0.96) ** (-0.19)
    )


def pkSZ(dZ):
    return 1.445 * (dZ / 1.05) ** 0.51


def cosmo_scale_tsz(H0, sigma8, omegab):
    return (H0 / 71.0) ** 1.73 * (sigma8 / 0.8) ** 8.34 * (omegab / 0.044) ** 2.81


def read_dl_template(filename, lmax, rescale_3000=True):
    """
    Read SPT template
    WARNING: need to check file before reading...
    """
    template = np.zeros(lmax + 1)
    with open(filename, "r") as f:
        for line in f:
            ll, realtmp = map(float, line.split())
            if ll >= 2 and ll <= lmax:
                template[int(ll)] = realtmp

    if rescale_3000:
        template = template / template[3000]

    return template


# proportional to the Planck function normalized to 1 at nu0
def Bnu(nu, nu0, T):
    # h/k
    # 4.799237 . 10-11 s K
    # expect GHz
    # so 4.799237e-2 K/GHz
    hk = 4.799237e-2

    Bnu = (nu / nu0) ** 3
    return Bnu * (np.exp(hk * nu0 / T) - 1.0) / (np.exp(hk * nu / T) - 1.0)


def dBdT(nu, nu0):
    """
    nu,nu0 in GHz
    dBdT is proportional to derivative of planck function
    but is normalized so its equal to 1 at nu0
    """
    d_plck = lambda x: x ** 4 * np.exp(x) / (np.exp(x) - 1) ** 2

    dBdT0 = d_plck(nu0 / 56.78)
    dBdT = d_plck(nu / 56.78)
    return dBdT / dBdT0


def tszFreqDep(nu, nu0):
    tsz = lambda x: x * (np.exp(x) + 1) / (np.exp(x) - 1) - 4

    tszFreqDep0 = tsz(nu0 / 56.78)
    tszFreqDep = tsz(nu / 56.78)

    return tszFreqDep / tszFreqDep0


class SPTforegounds(HasLogger):
    lmax = 13500

    def __init__(
        self,
        data_folder="",
        spt_dataset_clustered="",
        spt_dataset_clustered2="",
        spt_dataset_kSZ="",
        spt_dataset_kSZ2="",
        spt_dataset_tSZ="",
        shang_model_correlation_shape=True,
        only_1halo_tsz_cib=False,
        radio_ampl_mean=1.3,
        radio_ampl_unc=-1,
        single_clustered_freq_scaling=True,
        relative_alpha_cluster=False,
        apply_prior_cirrus_90ghz=True,
        apply_prior_cirrus_150ghz=True,
        apply_prior_cirrus_220ghz=True,
        spt_prior_clusterpoisson=-1.0,
    ):
        self.set_logger()

        self.relative_alpha_cluster = relative_alpha_cluster
        self.radio_ampl_mean = radio_ampl_mean
        self.radio_ampl_unc = radio_ampl_unc
        self.only_1halo_tsz_cib = only_1halo_tsz_cib

        self.single_clustered_freq_scaling = single_clustered_freq_scaling
        self.shang_model_correlation_shape = shang_model_correlation_shape

        self.apply_prior_cirrus_90ghz = apply_prior_cirrus_90ghz
        self.apply_prior_cirrus_150ghz = apply_prior_cirrus_150ghz
        self.apply_prior_cirrus_220ghz = apply_prior_cirrus_220ghz

        self.spt_prior_clusterpoisson = spt_prior_clusterpoisson

        #        InitForegroundData(SPTClusTemplate, SPTClus2Template, SPTkSZTemplate,SPTkSZ2Template,SPTtSZTemplate,self.spt_prior_clusterpoisson, relative_alpha_cluster)
        # clustered template
        self.clust_dg_templ = read_dl_template(
            os.path.join(data_folder, spt_dataset_clustered), self.lmax
        )
        self.clust2_dg_templ = read_dl_template(
            os.path.join(data_folder, spt_dataset_clustered2), self.lmax
        )

        # KSZ template
        self.ksz_templ = read_dl_template(os.path.join(data_folder, spt_dataset_kSZ), self.lmax)

        # 2nd KSZ template
        self.ksz2_templ = read_dl_template(os.path.join(data_folder, spt_dataset_kSZ2), self.lmax)

        # TSZ template
        self.tsz_templ = read_dl_template(os.path.join(data_folder, spt_dataset_tSZ), self.lmax)

        self.log.debug("1halo (or all clustered) DG template:{}".format(spt_dataset_clustered))
        self.log.debug("2halo (or 2nd) DG template:{}".format(spt_dataset_clustered2))
        self.log.debug("kSZ template:{}".format(spt_dataset_kSZ))
        self.log.debug("ksz template2:{}".format(spt_dataset_kSZ2))
        self.log.debug("tSZ template:{}".format(spt_dataset_tSZ))

        # prior on separation between clsutered and Poisson terms
        # ignored if <= 0
        self.AddAlphaPoisson = 1 if relative_alpha_cluster else 0

        self.l_divide_3000 = np.array([i / 3000.0 for i in range(self.lmax + 1)])

        # Cirrus template
        with np.errstate(divide="ignore"):
            self.cirrus_templ = self.l_divide_3000 ** (-1.2)
            self.cirrus_templ[:200] = 0.0

        # Cluster RG template
        self.clust_rg_templ = self.l_divide_3000 ** (0.6)
        self.clust_rg_templ[:50] = 0.0

    def ReportFGLmax(self):
        return self.lmax

    def dl_foregrounds(
        self, params, ifr, jfr, nfr, eff_fr, norm_fr, req_lmin, req_lmax, components=False
    ):
        """
        Compute dl foregrounds for SPT likelihood

        Input:
            params: dict of input parameters ()
            ifr: freq 1
            jfr: freq 2
            nfr: number of frequencies
            eff_fr: effective frequency
            nomr_fr:
            req_lmin: requested lmin
            req_lmax: requested lmax

        Output:
            dl_foreground: vector of fg
        """

        if req_lmax > self.lmax:
            raise ValueError("asked for too large ell in dl_foreground")
        if req_lmin < 2:
            raise ValueError("asked for too small ell in dl_foreground")

        dl_dg_cl = self._dl_dusty_clustered(
            params, eff_fr[ifr, 0], eff_fr[jfr, 0], norm_fr[0], ifr, jfr
        )
        dl_dg_po = self._dl_dusty_poisson(
            params, eff_fr[ifr, 1], eff_fr[jfr, 1], norm_fr[1], ifr, jfr
        )
        dl_rg_po = self._dl_radio(params, eff_fr[ifr, 2], eff_fr[jfr, 2], norm_fr[2])
        dl_k_sz = self._dl_ksz(params)
        dl_t_sz = self._dl_tsz(params, eff_fr[ifr, 4], eff_fr[jfr, 4], norm_fr[4])

        # TSZ-Dusty correlation
        if self.only_1halo_tsz_cib:
            dltmpi = self._dl_dusty_clustered(
                params, eff_fr[ifr, 0], eff_fr[ifr, 0], norm_fr[0], ifr, jfr, only1halo=True
            )
            dltmpj = self._dl_dusty_clustered(
                params, eff_fr[jfr, 0], eff_fr[jfr, 0], norm_fr[0], ifr, jfr, only1halo=True
            )
        else:
            dltmpi = self._dl_dusty_clustered(
                params, eff_fr[ifr, 0], eff_fr[ifr, 0], norm_fr[0], ifr, jfr
            ) + self._dl_dusty_poisson(params, eff_fr[ifr, 1], eff_fr[ifr, 1], norm_fr[1], ifr, jfr)
            dltmpj = self._dl_dusty_clustered(
                params, eff_fr[jfr, 0], eff_fr[jfr, 0], norm_fr[0], ifr, jfr
            ) + self._dl_dusty_poisson(params, eff_fr[jfr, 1], eff_fr[jfr, 1], norm_fr[1], ifr, jfr)

        dl_tsz_dgcl_cor = -self._tsz_dgcl_cor(params) * (
            self._tsz_cib(params)
            * np.sqrt(self._dl_tsz(params, eff_fr[ifr, 4], eff_fr[ifr, 4], norm_fr[4]) * dltmpj)
            + self._tsz_cib(params)
            * np.sqrt(self._dl_tsz(params, eff_fr[jfr, 4], eff_fr[jfr, 4], norm_fr[4]) * dltmpi)
        )

        # TSZ-Radio correlation
        if params["tsz_rg_cor"] == 0:
            dl_tsz_rg_cor = np.zeros(self.lmax + 1)
        else:
            dl_tsz_rg_cor = (
                -params["tsz_rg_cor"]
                * self._tsz_rg_cor()
                * (
                    sqrt(
                        self._dl_tsz(params, eff_fr[ifr, 4], eff_fr[ifr, 4], norm_fr[4])
                        * self._dl_radio(params, eff_fr[jfr, 0], eff_fr[jfr, 0], norm_fr[0])
                    )
                    + sqrt(
                        self._dl_tsz(params, eff_fr[jfr, 4], eff_fr[jfr, 4], norm_fr[4])
                        * self._dl_radio(params, eff_fr[ifr, 0], eff_fr[ifr, 0], norm_fr[0])
                    )
                )
            )

        # Cirrus
        dl_cirrus = self._dl_galcirrus(params, eff_fr[ifr, 0], eff_fr[jfr, 0])

        dl_dg = dl_dg_po + dl_dg_cl
        dl_rg = dl_rg_po

        if components:
            return (
                dl_dg_po,
                dl_dg_cl,
                dl_rg,
                dl_t_sz,
                dl_k_sz,
                dl_cirrus,
                dl_tsz_dgcl_cor,
                dl_tsz_rg_cor,
            )
        else:
            return dl_dg + dl_rg + dl_t_sz + dl_k_sz + dl_cirrus + dl_tsz_dgcl_cor + dl_tsz_rg_cor

    def _dl_dusty_clustered(self, params, fri, frj, fr0, ifr, jfr, only1halo=False):
        effalpha_cl = params["T_dg_cl"] + self.AddAlphaPoisson * params["T_dg_po"]
        effsigmasq_cl = params["sigmasq_dg_cl"] + self.AddAlphaPoisson * params["sigmasq_dg_po"]
        effalpha_cl_2 = params["beta_dg_cl"] + self.AddAlphaPoisson * params["beta_dg_po"]

        frqdep = ((fri * frj) / (fr0 * fr0)) ** (effalpha_cl_2)
        frqdep = frqdep * Bnu(fri, fr0, effalpha_cl) * Bnu(frj, fr0, effalpha_cl)
        frqdep = frqdep / dBdT(fri, fr0) / dBdT(frj, fr0)

        # not <1 like decor in the other case...this boosts the autospectra
        decor = ((fri * frj) / (fr0 * fr0)) ** (np.log(fri / fr0 * frj / fr0) / 2 * effsigmasq_cl)
        frqdep = frqdep * decor

        dl_dusty_clustered = params["czero_dg_cl"] * frqdep * self.clust_dg_templ

        #       if params["dg_cl_ell_power"] != 0:
        #           dl_dusty_clustered = dl_dusty_clustered * (self.l_divide_3000)**params["dg_cl_ell_power"]

        if not only1halo and params["czero_dg_cl2"] != 0:
            if self.single_clustered_freq_scaling:
                dl_dusty_clustered = (
                    dl_dusty_clustered + params["czero_dg_cl2"] * frqdep * self.clust2_dg_templ
                )
            else:
                effalpha_cl = params["T_dg_cl2"] + self.AddAlphaPoisson * params["T_dg_po"]
                effsigmasq_cl = (
                    params["sigmasq_dg_cl2"] + self.AddAlphaPoisson * params["sigmasq_dg_po"]
                )
                effalpha_cl_2 = params["beta_dg_cl2"] + self.AddAlphaPoisson * params["beta_dg_po"]
                frqdep = ((fri * frj) / (fr0 * fr0)) ** (effalpha_cl_2)
                frqdep = frqdep * Bnu(fri, fr0, effalpha_cl) * Bnu(frj, fr0, effalpha_cl)
                frqdep = frqdep / dBdT(fri, fr0) / dBdT(frj, fr0)
                decor = ((fri * frj) / (fr0 * fr0)) ** (
                    np.log(fri / fr0 * frj / fr0) / 2 * effsigmasq_cl
                )
                frqdep = frqdep * decor
                dl_dusty_clustered = dl_dusty_clustered + self.clust2_dg_templ * (
                    params["czero_dg_cl2"] * frqdep
                )

        return dl_dusty_clustered

    def _dl_dusty_poisson(self, params, fri, frj, fr0, ifr, jfr):
        # basic powerlaw
        frqdep = ((fri * frj) / (fr0 * fr0)) ** (params["beta_dg_po"])

        # not <1 like decor in the other case...this boosts the autospectra
        decor = ((fri * frj) / (fr0 * fr0)) ** (
            np.log(fri / fr0 * frj / fr0) / 2 * params["sigmasq_dg_po"]
        )
        frqdep = frqdep * decor

        frqdep = frqdep * Bnu(fri, fr0, params["T_dg_po"]) * Bnu(frj, fr0, params["T_dg_po"])

        dl_dusty_poisson = (
            params["czero_dg_po"] / dBdT(fri, fr0) / dBdT(frj, fr0) * frqdep
        ) * self.l_divide_3000 ** 2

        return dl_dusty_poisson

    def _dl_radio(self, params, fri, frj, fr0):
        dl_radio = (
            params["czero_rg_po"]
            / dBdT(fri, fr0)
            / dBdT(frj, fr0)
            * (fri / fr0 * frj / fr0)
            ** (params["alpha_rg"] + np.log(fri / fr0 * frj / fr0) / 2 * params["sigmasq_rg"] ** 2)
        ) * self.l_divide_3000 ** 2

        if params["czero_rg_cl"] > 0:
            dl_radio = (
                dl_radio
                + (
                    params["czero_rg_cl"]
                    / dBdT(fri, fr0)
                    / dBdT(frj, fr0)
                    * (fri / fr0 * frj / fr0) ** (params["alpha_rg"])
                )
                * self.clust_rg_templ
            )

        return dl_radio

    def _dl_ksz(self, params):
        dl_ksz = params["czero_ksz"] * self.ksz_templ
        if params["czero_ksz2"] != 0:
            dl_ksz = dl_ksz + params["czero_ksz2"] * self.ksz2_templ
        return dl_ksz

    def _dl_tsz(self, params, fri, frj, fr0):
        return (params["czero_tsz"] * tszFreqDep(fri, fr0) * tszFreqDep(frj, fr0)) * self.tsz_templ

    def _tsz_dgcl_cor(self, params):
        # motivated by shaw analysis of Sehgal sims
        #        tsz_dgcl_cor = max(0., (0.3 - 0.2 * exp(-(l-500)/1000))/0.3)

        # motivated by simplicity
        if self.shang_model_correlation_shape:
            tsz_dgcl_cor = -0.0703 * (self.l_divide_3000 ** 2) + 0.612 * self.l_divide_3000 + 0.458
        else:
            tsz_dgcl_cor = 1.0
            if params["tsz_cib_slope"] != 0:
                tsz_dgcl_cor = tsz_dgcl_cor + (self.l_divide_3000 - 1) * params["tsz_cib_slope"]

        return tsz_dgcl_cor

    def _tsz_cib(self, params):
        return params["tsz_dg_cor"]

    def _tsz_rg_cor(self):
        return 1.0

    def _dl_galcirrus(self, params, fri, frj):
        power = self._cirrus_power3000(params, fri, frj)
        return power * self.cirrus_templ

    def _cirrus_power3000(self, params, fri, frj):
        fr0 = 220.0
        frqdep = ((fri * frj) / (fr0 * fr0)) ** (params["beta_cirrus"])
        frqdep = frqdep * Bnu(fri, fr0, params["T_cirrus"]) * Bnu(frj, fr0, params["T_cirrus"])

        frqdep = frqdep / (dBdT(fri, fr0) * dBdT(frj, fr0))

        return params["czero_cirrus"] * frqdep

    def getForegroundPriorLnL(self, foregrounds):
        """
        foregrounds = dict of params
        output: chi2/2
        """

        PriorLnL = 0

        # dusty index prior
        del_alpha = self.spt_prior_clusterpoisson
        if del_alpha > 0:
            PriorLnL += (foregrounds["T_dg_po"] - foregrounds["T_dg_cl"]) ** 2 / (
                2 * (del_alpha) ** 2
            )

        # cirrus prior
        cirrus_factor = 1
        if foregrounds["czero_cirrus"] != 0:
            if self.apply_prior_cirrus_90ghz:
                freq = 97.9
                cirrus90 = self._cirrus_power3000(foregrounds, freq, freq)
                prior90 = 0.16 * cirrus_factor
                CirrusPriorLnL = (cirrus90 - prior90) ** 2 / (2 * (0.06) ** 2)
                self.log.debug("post cirrus90 prior:{}".format(CirrusPriorLnL))
                PriorLnL = PriorLnL + CirrusPriorLnL

            if self.apply_prior_cirrus_150ghz:
                freq = 153.44
                cirrus150 = self._cirrus_power3000(foregrounds, freq, freq)
                prior150 = 0.21 * cirrus_factor
                CirrusPriorLnL = (cirrus150 - prior150) ** 2 / (2 * (0.06) ** 2)
                self.log.debug("post cirrus150 prior:{}".format(CirrusPriorLnL))
                PriorLnL = PriorLnL + CirrusPriorLnL

            if self.apply_prior_cirrus_220ghz:
                freq = 219.67
                cirrus220 = self._cirrus_power3000(foregrounds, freq, freq)
                prior220 = 2.19 * cirrus_factor
                CirrusPriorLnL = (cirrus220 - prior220) ** 2 / (2 * (0.7) ** 2)
                self.log.debug("post cirrus220 prior:{}".format(CirrusPriorLnL))
                PriorLnL = PriorLnL + CirrusPriorLnL

        # Radio-Galaxies Poisson prior
        radio_amp = self.radio_ampl_mean
        radio_unc = self.radio_ampl_unc
        if radio_amp >= 0 and radio_unc > 0:
            fradio = (foregrounds["czero_rg_po"] - radio_amp) / radio_unc
            RGprior = fradio ** 2 / 2
            if foregrounds["czero_rg_po"] < 0:
                RGprior = RGprior + 1e6
            self.log.debug("RG poisson prior:{}".format(RGprior))
            PriorLnL = PriorLnL + RGprior

        # Radio-Galaxies Cluster prior
        #        if radio_cl_amp >= 0 and radio_cl_unc > 0:
        #            fradio = (foregrounds["czero_rg_cl"] - radio_cl_amp) / radio_cl_unc
        #            RGprior = fradio**2 / 2
        #            if foregrounds["czero_rg_cl"] < 0: RGprior = RGprior + 1e6
        #            self.log.debug("RG cluster prior:{}".format(RGprior))
        #            PriorLnL = PriorLnL + RGprior

        self.log.debug("Foreground prior:{}".format(PriorLnL))
        return PriorLnL
