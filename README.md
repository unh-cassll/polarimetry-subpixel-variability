# polarimetry-subpixel-variability

**A project in which we investigate the effects of subpixel variability on the polarimetric slope sensing of ocean waves.**

Manuscript in preparation for *IEEE Transactions on Geoscience and Remote Sensing* (*TGRS*):

"The Effects of Subpixel Variability on Polarimetric Sensing of Ocean Waves"
Nathan J. M. Laxague, Z. Göksu Duvarci, Lindsay Hogan, Steven P. Anderson, and Christopher J. Zappa

## Contents

**Figure-generating codes**

Start with *aa_polarimetry_subpixel_variability_figure_gen.m* ("aa" in filename pushes it to top of alphabetically-sorted lists)

... rest of codes in root directory are MATLAB functions called by the *figure_gen* script.

**Codes used in preparation of data**

These live within *_codes*. My implementation of the Elfouhaily *et al.* [1997] model wavenumber directional spectrum is provided in the *ElfouhailyEtAl1997* subdirectory.

**Data used to produce graphics**

These files live within *_data*. This folder is meant to be the smallest collection of example data which could be feasibly used for demonstration/testing purposes. The original raw datasets from which the spectra/statistics were produced are substantial in size (of order 300 TB).

**Graphics**

These images live within *_outputs/manuscript_files*. The *figure_gen* script in the root directory handles the printing of these to file, with user-set option to print as vector (.svg) or raster (.png) images.
