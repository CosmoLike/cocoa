import os
import warnings
import numpy as np
from cobaya.likelihood import Likelihood


# ---- default data files ----

def get_hd_filenames(delensed=True, baryonic_feedback=False):
    """Returns the file names of the CMB-HD lensed or delensed data."""
    # default data directory, relative to this file:
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data/')
    data_path = lambda fname: os.path.join(data_dir, fname)
    # file names
    cmb_type = 'delensed' if delensed else 'lensed'
    ell_info = 'lmin30_lmax20k'
    Linfo = 'Lmin30_Lmax20k'
    bin_file = data_path('bin_edges.txt')
    if baryonic_feedback:
        data_file = data_path(f'hd_binnedTheorySpectra_{ell_info}_{Linfo}_{cmb_type}_feedback.txt')
    else:
        data_file = data_path(f'hd_binnedTheorySpectra_{ell_info}_{Linfo}_{cmb_type}.txt')
    covmat_file = data_path(f'hd_covmat_{ell_info}_{Linfo}_{cmb_type}_withFG.txt')
    recon_noise_file = data_path(f'hd_lensingNoise_{Linfo}_withFG.txt')
    return bin_file, data_file, covmat_file, recon_noise_file


def get_desi_filenames():
    """Returns the filenames of the mock DESI BAO data."""
    # default data directory, relative to this file:
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data/')
    data_path = lambda fname: os.path.join(data_dir, fname)
    data_file = data_path('mock_desi_bao_rs_over_DV_data.txt')
    cov_file = data_path('mock_desi_bao_rs_over_DV_cov.txt')
    return data_file, cov_file


# ---- binning ----

def get_bin_info(bin_edges, lmax, lmin=2):
    """Given a 1D array of the bin edges, returns three arrays containing
    the lower and upper edges of each bin and the corresponding bin centers,
    for bins in the given ell-range.

    Parameters
    ----------
    bin_edges : array_like of int
        A 1D array containing the edges of the bins.
    lmax, lmin : int
        Sets the ell-range used.

    Returns
    -------
    lower, upper : array_like of int
        1D Arrays containing the lower or upper edge of each bin.
    center : array_like of float
        A 1D array containing the center of each bin.
    """
    edges = bin_edges.astype(int)
    lower = edges[:-1].copy() # all but the last bin edge
    upper = edges[1:].copy() # all but the first bin edge
    # currently, upper[i] = lower[i+1]; want the i^th bin to include both 
    # edges (upper[i] and lower[i] are both in bin i), so add one to all 
    # lower edges (except the first)
    lower[1:] += 1
    # trim the bin edges to be in the correct ell-range
    lloc = np.where(lower >= int(lmin))
    lower = lower[lloc]
    upper = upper[lloc]
    uloc = np.where(upper <= int(lmax))
    lower = lower[uloc]
    upper = upper[uloc]
    # get the bin centers
    center = (lower + upper) / 2.
    return lower, upper, center


def get_binning_matrix(bin_edges, lmin=2, lmax=None, start_at_ell=2):
    lmax = int(lmax) if (lmax is not None) else int(bin_edges[-1])
    lmin = int(lmin)
    # get ells from starting value to lmax
    ells = np.arange(int(start_at_ell), lmax+1)
    num_ells = len(ells)
    # get lower and upper bin edges
    lower, upper, _ = get_bin_info(bin_edges, lmax, lmin=lmin)
    nbin = len(lower)
    # make the binninng matrix
    binmat = np.zeros((nbin, num_ells))
    for i, (bmin, bmax) in enumerate(zip(lower, upper)):
        loc = np.where((ells >= bmin) & (ells <= bmax))
        n = bmax - bmin + 1 # number of ells in this bin
        binmat[i][loc] = 1. / n
    return binmat


# functions to decompose full covmat into blocks, and put it back together
# (used during initialization)

def nbins_per_spectrum(ell_ranges, bin_edges):
    """Given the minimum and maximum multipoles used for each spectrum, along
    with the bin edges, returns a dictionary holding the number of bins used
    for each spectrum.

    Parameters
    ----------
    ell_ranges : dict of list or tuple of int
        A dictionary with keys given by the elements of `spectra`, holding
        a tuple or list `[lmin, lmax]` giving the minimum and maximum 
        multipoles used for each spectra. 
    bin_edges : array_like of int
        A 1D array holding the multipole values for the lower edge of the 
        first bin, and upper edges of each bin. 

    Returns
    -------
    nbins : dict of int
        A dictionary with the same keys as `ell_ranges`, holding the number
        of bins used for each spectrum.
    """
    nbins = {}
    for s, (lmin, lmax) in ell_ranges.items():
        _, _, lbin = get_bin_info(bin_edges, lmax, lmin=lmin)
        nbins[s] = len(lbin)
    return nbins


def cov_to_blocks(cov, spectra=['tt', 'te', 'ee', 'bb', 'kk'], ell_ranges=None, bin_edges=None):
    """Given a covariance matrix `cov` containing blocks for the covariance
    between the different `spectra` (e.g. a TT x TT block for the temperature power
    spectrum auto-covariance, a TT x kappakappa block for the temperature and
    lensing potential power spectra cross-covariance, etc.), returns a nested dict
    of the different blocks.

    Parameters
    ----------
    cov : array_like
        The full, two-dimensional covariane matrix for the given `spectra`.
    spectra : list of str, default=['tt', 'te', 'ee', 'bb', 'kk']
        A list of the spectra in the covariance matrix, in the correct order.
    ell_ranges : dict of list or tuple of int, default=None
        A dictionary with keys given by the elements of `spectra`, holding
        a tuple or list `[lmin, lmax]` giving the minimum and maximum
        multipoles used for each spectra. If `None`, it's assumed that
        each block of the covariance matrix has the same multipole range.
        Otherwise the `bin_edges` must be provided.
    bin_edges : array_like of int, default=None
        A 1D array holding the multipole values for the lower edge of the
        first bin, and upper edges of each bin. Must be passed with
        `ell_ranges`, if the blocks in the covariance matrix have different
        multipole ranges.

    Returns
    -------
    blocks : dict of dict of array_like
        A nested dict, with keys given by the elements in `spectra`, holding the
        two-dimensional covariance matrix for that block. The first key gives
        the row of the block, and the second gives the column. See the note below.

    Raises
    ------
    ValueError
        If the `ell_ranges` are provided but no `bin_edges` were given.

    Note
    ----
    For a number `nspec` of different spectra, there are `nspec * nspec` blocks
    in the covariance matrix. Using the default `spectra=['tt', 'te', 'ee', 'bb', 'kk']`,
    `nspec = 5`, so there are 25 blocks. The first row of blocks will be TT x TT,
    TT x TE, TT x EE, TT x BB, TT x kappakappa; the second will be TE x TT, TE x TE, etc.
    These are stored in the dict as `blocks['tt']['tt']` for TT x TT,
    `blocks['tt']['te']` for TT x TE, etc.

    See Also
    --------
    cov_from_blocks :
        The inverse function which takes in `blocks` and returns the `cov`.
    """
    nspec = len(spectra)
    # get number of bins for each spectrum
    if ell_ranges is None: # each block has same number of bins
        nbin = cov.shape[0] // nspec # per block
        nbins = {s: nbin for s in spectra}
    else: # each block may have different number of bins
        if bin_edges is None:
            raise ValueError('You must also provide the `bin_edges` with the `ell_ranges`.')
        else:
            nbins = nbins_per_spectrum(ell_ranges, bin_edges)
    blocks = {}
    for i1, s1 in enumerate(spectra):
        blocks[s1] = {}
        for i2, s2 in enumerate(spectra):
            # indices where this block starts and ends
            imin1 = sum([nbins[s] for s in spectra[:i1]])
            imin2 = sum([nbins[s] for s in spectra[:i2]])
            imax1 = imin1 + nbins[s1]
            imax2 = imin2 + nbins[s2]
            blocks[s1][s2] = cov[imin1:imax1, imin2:imax2].copy()
    return blocks


def cov_from_blocks(blocks, spectra=['tt', 'te', 'ee', 'bb', 'kk'], ell_ranges=None, bin_edges=None):
    """Given a nested dict `blocks` containing covariance matrices for the
    power spectra types in `spectra` (e.g., `blocks['tt']['ee']` is the
    TT x EE cross-covariance matrix for the temperature and polarization
    power spectra), return a single covariance matrix containing the blocks.

    Parameters
    ----------
    blocks : nested dict of array_like
        A nested dict, with keys given by the elements in `spectra`, holding the
        two-dimensional covariance matrix for that block. The first key gives
        the row of the block, and the second gives the column. See the note below.
    spectra : list of str, default=['tt', 'te', 'ee', 'bb', 'kk']
        A list of the spectra in the covariance matrix, in the correct order.
    ell_ranges : dict of list or tuple of int, default=None
        A dictionary with keys given by the elements of `spectra`, holding
        a tuple or list `[lmin, lmax]` giving the minimum and maximum
        multipoles used for each spectra. If `None`, it's assumed that
        each block of the covariance matrix has the same multipole range.
        Otherwise the `bin_edges` must be provided.
    bin_edges : array_like of int, default=None
        A 1D array holding the multipole values for the lower edge of the
        first bin, and upper edges of each bin. Must be passed with
        `ell_ranges`, if the blocks in the covariance matrix have different
        multipole ranges.

    Returns
    -------
    cov : array_like
        The full, two-dimensional covariance matrix for the given `spectra`.

    Raises
    ------
    ValueError
        If the `ell_ranges` are provided but no `bin_edges` were given.

    Note
    ----
    For a number `nspec` of different spectra, there are `nspec * nspec` blocks
    in the covariance matrix. Using the default `spectra=['tt', 'te', 'ee', 'bb', 'kk']`,
    `nspec = 5`, so there are 25 blocks. The first row of blocks will be TT x TT,
    TT x TE, TT x EE, TT x BB, TT x kappakappa; the second will be TE x TT, TE x TE, etc.
    These should be stored in `blocks` as `blocks['tt']['tt']` for TT x TT,
    `blocks['tt']['te']` for TT x TE, etc.

    See Also
    --------
    cov_to_blocks :
        The inverse function which takes in the `cov` and returns `blocks`.
    """
    nspec = len(spectra)
    # get number of bins for each spectrum
    if ell_ranges is None: # each block has same number of bins
        nbin = blocks[spectra[0]][spectra[0]].shape[0] # per block
        nbins = {s: nbin for s in spectra}
        nbin_tot = nbin * nspec
    else: # each block may have different number of bins
        if bin_edges is None:
            raise ValueError('You must also provide the `bin_edges` with the `ell_ranges`.')
        else:
            nbins = nbins_per_spectrum(ell_ranges, bin_edges)
            nbin_tot = sum([nbins[s] for s in spectra])
    cov = np.zeros((nbin_tot, nbin_tot))
    for i1, s1 in enumerate(spectra):
        for i2, s2 in enumerate(spectra):
            # indices where this block starts and ends
            imin1 = sum([nbins[s] for s in spectra[:i1]])
            imin2 = sum([nbins[s] for s in spectra[:i2]])
            imax1 = imin1 + nbins[s1]
            imax2 = imin2 + nbins[s2]
            cov[imin1:imax1, imin2:imax2] = blocks[s1][s2].copy()
    return cov


class HDData:
    def __init__(self, lmin=30, lmax=20100, Lmax=20100, delensed=True,
                 baryonic_feedback=False, data_file=None, covmat_file=None, 
                 bin_file=None, recon_noise_file=None, 
                 has_cmb_power_spectra=True, 
                 has_cmb_lensing_spectrum=True, 
                 use_cmb_power_spectra=True, 
                 use_cmb_lensing_spectrum=True,
                 use_desi_bao=False): #NOTE: set `use_desi_bao=False` when using Cobaya 
        """Initialize the CMB-HD likelihood with the binned lensed or delensed 
        data spectra and covariance matrix.

        Parameters
        ----------
        lmin : int, default=30
            The minimum multipole used in the data files, which sets the lowest
            bin edge. All CMB power spectra and CMB lensing spectra, and the 
            corresponding covariance matrix, are assumed to start at the bin edge 
            corresponding to this value.
        lmax, Lmax : int, default=20100
            The maximum multipole used in the CMB power spectra (`lmax`) and CMB 
            lensing power spectrum (`Lmax`). These can be set to a lower value than 
            the default, which will cut the mock spectra and covariance matrix. 
        delensed : bool, default=True
            Whether to use delensed data (binned CMB power spectra and covariance
            matrix for the delensed case), as opposed to lensed data.
        baryonic_feedback: bool, default=False
            Whether to use binned power spectra that were calculated using the
            HMCode2020 non-linear model that includes the effect of baryonic
            feedback, as opposed to the HMCode2016 CDM-only model.
        data_file : str, default=None
            The path to and name of the file containing the binned power 
            spectra as a single one-dimensional array. 
            If not specified, the default file is used.
            The CMB TT, TE, EE, and BB spectra are expected to be in units 
            of uK^2, without any multipole factor of ell * (ell + 1) / 2pi 
            applied. The CMB lensing spectrum is expected to be in the form 
            C_L^kk = [L(L+1)]^2 C_L^phiphi / 4, where L is the lensing 
            multipole and C_L^phiphi is the power spectrum of the projected 
            lensing potential. The order of the spectra is expected to be
            TT, TE, EE, BB, kk.
        cov_file : str, default=None
            The path to and name of the file containing the binned covariance
            matrix. If not specified, the default file is used.
            The covariance matrix should have blocks for the covariance between
            different spectra, i.e. TT x TT for cov(C_l1^TT, C_l2^TT). These
            blocks should be in the same order and binned in the same way as
            the `data_file`.
        bin_file : str, default=None
            The path to and name of the file containing the bin edges, used to
            bin the power spectra and covariance matrix. The file should contain 
            a one-dimensional array of lower bin edges, except the last entry, 
            which should be the upper edge of the last bin. If not specified, 
            the default is used.
        recon_noise_file : str, default=None
            The path  to and name of the file containing the unbinned lensing 
            reconstruction noise spectrum, in the same convention as the 
            lensing power spectrum. The first column in the file should hold 
            the lensing multipoles, and the second should hold the value of the
            noise at that multipole. This should be the iteratively delensed 
            residual lensing noise spectrum and is only used to obtain delensed
            theory spectra from CAMB. If not specified, the default file is 
            used.
        has_cmb_power_spectra : bool, default=True
            Whether the binned spectra and covariance matrix have blocks for 
            the TT, TE, EE, and BB CMB power spectra.
        has_cmb_lensing_spectrum : bool, default=True
            Whether the binned spectra and covariance matrix have (a) block(s)
            for the lensing 'kk' power spectrum.
        use_cmb_power_spectra : bool, default=True
            Whether to include the TT, TE, EE, and BB CMB power spectra in
            the likelihood calculation, if applicable.
        use_cmb_lensing_spectrum : bool, default=True
            Whether to include the CMB lensing power spectrum in the likelihood
            calculation, if applicable.
        use_desi_bao : bool, default=False
            Whether to load in the mock DESI BAO data. Note that the likelihood 
            calculation for the mock BAO data is separate from the calculation for
            the mock CMB data.

        Raises
        ------
        ValueError
            If the settings for `use_cmb_power_spectra` and 
            `use_cmb_lensing_spectrum` are both `False`, or if either is
            inconsistent with the settings `has_cmb_power_spectra` and 
            `has_cmb_lensing_spectrum`.


        Note
        ----
        You do not have to pass file names if you are using the mock data and
        covariance matrices provided with `hdlike`; they are found automatically 
        in the functions `hdlike.get_hd_filenames` and `hdlike.get_desi_filenames`. 
        The option to use alternative mock data is included for flexibility, but 
        you must ensure that they are binned consistently with the `hdlike` 
        covariance matrix, have the correct ordering, and follow all other conventions 
        stated above.
        """
        self.lmin = lmin
        self.lmax = lmax
        self.Lmax = Lmax
        # --- default multipole ranges for full data set ---
        self.hd_lmin = 30
        self.hd_lmax = 20100
        self.hd_Lmax = 20100
        # --- check input ---
        if (not use_cmb_power_spectra) and (not use_cmb_lensing_spectrum):
            errmsg = "You set both `use_cmb_power_spectra` and `use_cmb_lensing_spectrum` to `False`, so there is nothing to calculate!"
            raise ValueError(errmsg)
        if use_cmb_power_spectra and (not has_cmb_power_spectra):
            errmsg = "You set `use_cmb_power_spectra: True` and `has_cmb_power_spectra: False`. To use the CMB TT/TE/EE/BB data, you must also set `has_cmb_power_spectra: True."
            raise ValueError(errmsg)
        if use_cmb_lensing_spectrum and (not has_cmb_lensing_spectrum):
            errmsg = "You set `use_cmb_lensing_spectrum: True` and `has_cmb_lensing_spectrum: False`. To use CMB lensing data, you must also set `has_cmb_lensing_spectrum: True."
            raise ValueError(errmsg)
        # default file names
        default_bin_file, default_data_file, default_covmat_file, default_recon_noise_file = get_hd_filenames(delensed=delensed, baryonic_feedback=baryonic_feedback)
        # if `delensed = True`, and the user provides either a new `data_file`
        # or a new `recon_noise_file` (used for delensed theory) but not both,
        # warn the user that their theory calculation may not match the data
        if delensed:
            warn = False
            if (data_file is None) and (recon_noise_file is not None):
                warn = True
                msg = "You provided a `recon_noise_file` to calculate the delensed theory, but you're using the default `data_file`, so your delensed data and theory may be inconsistent."
            elif (data_file is not None) and (recon_noise_file is None):
                warn = True
                msg = "You provided a `data_file` for the delensed data, but you're using the default `recon_noise_file` to calculate the delensed theory, so your delensed data and theory may be inconsistent. (You may ignore this message if you're using an automatically-generated YAML file)."
            if warn and (os.path.basename(data_file) != os.path.basename(default_data_file)):
                warnings.warn(msg)
        # also warn the user if `delensed=True` but `use_cmb_power_spectra=False`
        if delensed and (not use_cmb_power_spectra):
            warnings.warn("You set `delensed = True` but `use_cmb_power_spectra = False`, so there is nothing to delens.")
        self.has_cmb_power_spectra = has_cmb_power_spectra
        self.use_cmb_power_spectra = use_cmb_power_spectra
        self.has_cmb_lensing_spectrum = has_cmb_lensing_spectrum
        self.use_cmb_lensing_spectrum = use_cmb_lensing_spectrum
        self.delensed = delensed
        # --- load the data ---
        bin_fname = default_bin_file if (bin_file is None) else bin_file
        data_fname = default_data_file if (data_file is None) else data_file
        covmat_fname  = default_covmat_file if (covmat_file is None) else covmat_file
        recon_noise_fname = default_recon_noise_file if (recon_noise_file is None) else recon_noise_file
        self.desi = use_desi_bao
        if self.desi:
            desi_data_file, desi_cov_file = get_desi_filenames()
            desi_cov = np.loadtxt(desi_cov_file)
            self.desi_invcov = np.linalg.inv(desi_cov)
            self.z, self.rs_dv = np.loadtxt(desi_data_file, unpack=True, usecols=(0,1))
        # load the data, covmat, and bin edges
        data = np.loadtxt(data_fname)
        covmat = np.loadtxt(covmat_fname)
        self.bin_edges = np.loadtxt(bin_fname)
        # get the number of bins and the binning matrix for cmb data
        cmb_lower, cmb_upper, _ = get_bin_info(self.bin_edges, self.lmax, lmin=self.lmin)
        self.cmb_nbin = len(cmb_lower)
        self.cmb_binmat = get_binning_matrix(self.bin_edges, lmin=self.lmin, lmax=self.lmax)
        # get the number of bins and the binning matrix for lensing data
        lens_lower, lens_upper, _ = get_bin_info(self.bin_edges, self.Lmax, lmin=self.lmin)
        self.lens_nbin = len(lens_lower)
        self.lens_binmat = get_binning_matrix(self.bin_edges, lmin=self.lmin, lmax=self.Lmax)
        # trim the data and covmat to keep bins below `lmax` or `Lmax`
        # NOTE that we assume data and covmat begin at `lmin`
        data, covmat = self.trim_data_lmax(data, covmat)
        # if we have unnecessary 'blocks' (e.g., data has clkk, but only 
        # want cmb), remove them
        self.data, covmat = self.trim_data_blocks(data, covmat)
        self.invcov = np.linalg.inv(covmat)
        # load the lensing reconstruction noise, if we need it
        if self.delensed:
            self.L, self.recon_noise = np.loadtxt(recon_noise_fname, unpack=True)
            # will want recon noise up to CAMB's internal lmax, set to inf
            # outside range [lmin, lmax]; we don't yet know CAMB's lmax,
            # so set it to None for now, then on first iteration get the noise
            # in the correct format and store it here:
            self.nlkk = None


    def trim_data_lmax(self, data, covmat):
        """Given a binned data vector and covariance matrix, each consisting
        of different 'blocks', trim each block to include only bins below
        the desired lmax. Note that no trimming is done for any lmin."""
        # determine how many 'blocks' are in the data (one per spectrum),
        # and the total number of bins we want to keep
        num_bins = 0
        spectra = []
        if self.has_cmb_power_spectra:
            spectra = ['tt', 'te', 'ee', 'bb']
            num_bins += 4 * self.cmb_nbin
        if self.has_cmb_lensing_spectrum:
            spectra.append('kk')
            num_bins += self.lens_nbin
        num_blocks = len(spectra)
        # determine how many bins per block in data and covmat currently
        tot_data_nbin = len(data) 
        tot_cov_nbin = covmat.shape[0] 
        # trim each block and put them back together, if necessary
        if (tot_data_nbin == num_bins) and (tot_cov_nbin == num_bins):
            return data, covmat
        else:
            # multipole ranges for full data and covmat
            hd_ell_ranges = {s: [self.hd_lmin, self.hd_lmax] for s in ['tt', 'te', 'ee', 'bb']}
            hd_ell_ranges['kk'] = [self.hd_lmin, self.hd_Lmax]
            hd_nbins = nbins_per_spectrum(hd_ell_ranges, self.bin_edges)
            # break data vector into individual spectra and trim each
            data_blocks = {}
            imin = 0 # tmp starting index of the data block
            for s in spectra:
                if s in ['tt', 'te', 'ee', 'bb']:
                    imax = imin + self.cmb_nbin
                else:
                    imax = imin + self.lens_nbin
                data_blocks[s] = data[imin:imax].copy()
                imin += hd_nbins[s]
            trimmed_data = np.concatenate([data_blocks[s] for s in spectra])
            # multipole ranges for trimmed covmat
            ell_ranges = {s: [self.lmin, self.lmax] for s in ['tt', 'te', 'ee', 'bb']}
            ell_ranges['kk'] = [self.lmin, self.Lmax]
            nbins = nbins_per_spectrum(ell_ranges, self.bin_edges)
            # break the covmat into blocks and trim each
            cov_blocks = cov_to_blocks(covmat, spectra=spectra, ell_ranges=hd_ell_ranges, bin_edges=self.bin_edges)
            for i, s1 in enumerate(spectra):
                for s2 in spectra[i:]:
                    cov_blocks[s1][s2] = cov_blocks[s1][s2][:nbins[s1],:nbins[s2]]
                    if s1 != s2:
                        cov_blocks[s2][s1] = np.transpose(cov_blocks[s1][s2].copy())
            trimmed_covmat = cov_from_blocks(cov_blocks, spectra=spectra, ell_ranges=ell_ranges, bin_edges=self.bin_edges)
            return trimmed_data, trimmed_covmat


    def trim_data_blocks(self, data, covmat):
        """Given a full data vector and covariance matrix containing blocks for
        cltt, clte, clee, clbb, and clkk, removes the unneeded blocks (e.g., to
        only include CMB without lensing).
        """
        if self.has_cmb_lensing_spectrum and (not self.use_cmb_lensing_spectrum):
            return data[:4*self.cmb_nbin], covmat[:4*self.cmb_nbin, :4*self.cmb_nbin]
        elif self.has_cmb_power_spectra and (not self.use_cmb_power_spectra):
            return data[-self.lens_nbin:], covmat[-self.lens_nbin:, -self.lens_nbin:]
        else:
            return data, covmat


    def get_desi_redshifts(self):
        """Returns an array of redshifts for the DESI BAO mock data (e.g., to
        pass to CAMB).
        
        Returns
        -------
        z: array_like of float
            An array of redshifts.

        Raises
        ------
        ValueError
            If `use_desi_bao` was `False` in the initialization of `HDData`.
        """ 
        if self.desi:
            return self.z
        else:
            raise ValueError("You must set `use_desi_bao=True` when you initialize `HDData` to use the mock DESI BAO data.")


    def get_clkk_res(self, camb_results):
        """Calculate the residual lensing power, given the lensing reconstruction noise"""
        lmax = camb_results.Params.max_l # need CAMBs internal lmax
        if self.nlkk is None:
            # need recon noise to start at ell = 0, and set to inf outside range [lmin, lmax]
            nlkk = np.ones(lmax+1) * np.inf
            loc = np.where((self.L >= self.lmin) & (self.L <= self.Lmax))
            nlkk[self.lmin:self.Lmax+1] = self.recon_noise[loc].copy()
            self.nlkk = nlkk.copy() * 4. /(2. * np.pi) # need it in CAMB convention
        # get clkk up to CAMB's lmax, in CAMB convention
        clkk = camb_results.get_lens_potential_cls(lmax=lmax)[:,0]
        # wiener filter with the noise
        filt = clkk / (clkk + self.nlkk)
        clkk_filt = clkk * filt
        # set the wiener-filtered clkk to zero outside of the range `lmin`, `Lmax`
        clkk_filt[:self.lmin] = 0.
        clkk_filt[self.Lmax+1:] = 0.
        # return the residual lensing power, in CAMB convention
        clkk_res = clkk - clkk_filt
        return clkk_res


    def get_delensed(self, camb_results):
        """Get the delensed CMB power spectra from the residual lensing power."""
        clkk_res = self.get_clkk_res(camb_results)
        delensed_cls = camb_results.get_lensed_cls_with_spectrum(clkk_res, lmax=self.lmax, CMB_unit='muK', raw_cl=True)
        theo = {}
        for i, s in enumerate(['tt', 'ee', 'bb', 'te']):
            theo[s] = delensed_cls[:,i]
            theo[s][:2] = 0
        return theo


    def log_likelihood_desi(self, theo_rs_dv):
        """Calculate the DESI BAO log(likelihood) = -chi2 / 2, given the theory.

        Parameters
        ----------
        theo_rs_dv : array_like of float
            An array of the theoretical BAO measurement r_s/d_V(z), at the same
            redshifts as the DESI BAO data.

        Returns
        -------
        loglike : float
            The log-likelihood value.

        Raises
        ------
        ValueError
            If `use_desi_bao` was `False` in the initialization of `HDData`.

        Note
        ----
        You can access the redshifts for the DESI data by calling the function
        `hdlike.HDData.get_desi_redshifts`.
        """
        if self.desi:
            diff = self.rs_dv -  theo_rs_dv
            tmp = self.desi_invcov @ diff
            chi2 = np.transpose(diff) @ tmp
            loglike = -0.5 * chi2
            return loglike
        else:
            raise ValueError("You must set `use_desi_bao=True` when you initialize `HDData` to use the mock DESI BAO data.")



    def log_likelihood(self, theo_cls, camb_results=None):
        """Calculate the CMB-HD log(likelihood) = -chi2 / 2, given 
        the theory spectra.

        Parameters
        ----------
        theo_cls : dict of array_like of float
            A dictionary of the lensed CMB power spectra in units of uK^2, 
            with keys `'tt'`, `'te'`, `'ee'`, and `'bb'`; and/or the lensing
            potential power spectrum with key `'pp'`. None of the spectra 
            should have any ell-factors applied.
        camb_results : camb.results.CAMBdata, default=None
            A `camb.results.CAMBdata` instance, used to calculate the 
            delensed theory. This may be `None` if not using delensing.

        Returns
        -------
        loglike : float
            The log-likelihood value of the theory, given the data.
        """
        # dict of cmb spectra
        if self.delensed and self.use_cmb_power_spectra:
            delensed_theo_cls = self.get_delensed(camb_results)
            for s in ['tt', 'te', 'ee', 'bb']:
                theo_cls[s] = delensed_theo_cls[s] # replace lensed cmb with delensed
        if self.use_cmb_lensing_spectrum:
            # get clkk from clpp
            ells = np.arange(len(theo_cls['pp']))
            Lfact = (ells * (ells + 1.))**2 / 4.
            theo_cls['kk'] = theo_cls['pp'] * Lfact
            theo_cls['kk'][:2] = 0
        else:
            ells = np.arange(len(theo_cls['tt']))
        # bin the theory
        binned_theo = {}
        spectra = []
        if self.use_cmb_power_spectra:
            for s in ['tt', 'te', 'ee', 'bb']:
                spectra.append(s)
                binned_theo[s] = self.cmb_binmat @ theo_cls[s][2:self.lmax+1]
        if self.use_cmb_lensing_spectrum:
            spectra.append('kk')
            binned_theo['kk'] = self.lens_binmat @ theo_cls['kk'][2:self.Lmax+1]
        theo = np.concatenate([binned_theo[s] for s in spectra])
        # calculate loglike = -chi2 / 2
        diff = self.data - theo
        tmp = self.invcov @ diff 
        chi2 = np.transpose(diff) @ tmp
        loglike = -0.5 * chi2
        return loglike


class HDLike(Likelihood):
    def initialize(self):
        """Load the CMB-HD data and covariance matrix, and determine what's 
        in it (e.g., both CMB and lensing potential). Also set lmin/lmax 
        and load the bin edges to bin the theory in the same way as the data.

        Raises
        -----
        ValueError 
            If the settings to use CMB and/or lensing data are inconsistent.
        """
        self.hd_data = HDData(lmin=self.lmin, lmax=self.lmax, Lmax=self.Lmax, 
                              delensed=self.delensed, 
                              baryonic_feedback=self.baryonic_feedback,
                              data_file=self.data_file, 
                              covmat_file=self.covmat_file, 
                              bin_file=self.bin_file, 
                              recon_noise_file=self.recon_noise_file, 
                              has_cmb_power_spectra=self.has_cmb_power_spectra, 
                              has_cmb_lensing_spectrum=self.has_cmb_lensing_spectrum,
                              use_cmb_power_spectra=self.use_cmb_power_spectra, 
                              use_cmb_lensing_spectrum=self.use_cmb_lensing_spectrum)
    

    def get_requirements(self):
        """Returns a dict with the quantities needed by the theory code."""
        reqs = {'Cl': {}}
        if self.use_cmb_power_spectra:
            for s in ['tt', 'te', 'ee', 'bb']:
                reqs['Cl'][s] = self.lmax
        if self.use_cmb_lensing_spectrum:
            reqs['Cl']['pp'] = self.Lmax
        if self.delensed: # need the CAMB results
            reqs['CAMBdata'] = None
        return reqs


    def logp(self, **params_values):
        """Calculate the log(likelihood) = -chi2 / 2, given the theory"""
        # dict of lensed cmb spectra and the lensing potential spectrum
        theo_cls = self.provider.get_Cl()
        if self.delensed and self.use_cmb_power_spectra: # for delensed theory
            camb_results = self.provider.get_CAMBdata()
        else:
            camb_results = None
        loglike = self.hd_data.log_likelihood(theo_cls, camb_results=camb_results)
        return loglike


