# TS-lookup

- This is the code from Goes et al. (2018) to infer salinity for temperature profiles.
- This matlab package is to calculate the salinity from temperature profiles.

## REFERENCE

Goes, M., J. Christophersen, S. Dong, G. Goni, and M.O. Baringer, 2018: An Updated Estimate of Salinity for the Atlantic Ocean Sector Using Temperature–Salinity Relationships. J. Atmos. Oceanic Technol., 35, 1771–1784, <https://doi.org/10.1175/JTECH-D-18-0029.1>


## INSTALLATION

run: `csh install.sh`

## TEST

run in matlab: `example.m`

## Notes 

- Now it works for the whole globe between `70S-65N`.
- **NEEDS FOLDER WOA13 TO WORK**

## Contents of this package

- `Calc_sal_Thacker_Goes_EmDr_Stom_svd_globe.m` : Calculate the salinity from the T profiles.
- inputs:

```
        TT - temperature
        PP - pressure
        lat - latitude
        lon - longitude
        timet - time
        pad - use padding or not (1 and 0)
        Pout - vector os interpolated depths if wanted
        mehod - choose between Thacker [default], goes, annual, svd
```

- outputs:

```
         y2 - salinity (Now it is NaN free - linearly interpolated)
         y3 - smoothed version of y2
         TT - temperature profiles (can be interpolated is Pout provided)
         PP - new depth profile
```

- `example.m` : runs one example using `data_example.mat` and produces the following figures.

```
   TS_orig_reconstruct_map_%g_scatter_SA.eps
   TS_orig_reconstruct_SA.eps
   TS_orig_reconstruct_scatter_SA.eps
```

## Fitted parameters for each region in the globe

```
fit_sigma_cora_argo_extended_fitlm_MD_thackergoes_noyear.mat
fit_sigma_cora_argo_extended_fitlm_NA_thackergoes_noyear.mat
fit_sigma_cora_argo_extended_fitlm_NOP_thackergoes_noyear.mat
fit_sigma_cora_argo_extended_fitlm_NO_thackergoes_noyear.mat
fit_sigma_cora_argo_extended_fitlm_NP_thackergoes_noyear.mat
fit_sigma_cora_argo_extended_fitlm_SA_thackergoes_noyear.mat
fit_sigma_cora_argo_extended_fitlm_SOI_thackergoes_noyear.mat
fit_sigma_cora_argo_extended_fitlm_SOP_thackergoes_noyear.mat
fit_sigma_cora_argo_extended_fitlm_SO_thackergoes_noyear.mat
fit_sigma_cora_argo_extended_fitlm_SP_thackergoes_noyear.mat
fit_sigma_cora_argo_extended_fitlm_TA_thackergoes_noyear.mat
fit_sigma_cora_argo_extended_fitlm_TI_thackergoes_noyear.mat
fit_sigma_cora_argo_extended_fitlm_TP_thackergoes_noyear.mat
fit_sigma_cora_argo_extended_fitlm_FC_thackergoes_noyear_ryanonly.mat
```

## And using the Svd method

```
fit_sigma_cora_argo_extended_svd_MD.mat
fit_sigma_cora_argo_extended_svd_NA.mat
fit_sigma_cora_argo_extended_svd_NO.mat
fit_sigma_cora_argo_extended_svd_NOP.mat
fit_sigma_cora_argo_extended_svd_NP.mat
fit_sigma_cora_argo_extended_svd_SA.mat
fit_sigma_cora_argo_extended_svd_SOI.mat
fit_sigma_cora_argo_extended_svd_SO.mat
fit_sigma_cora_argo_extended_svd_SOP.mat
fit_sigma_cora_argo_extended_svd_SP.mat
fit_sigma_cora_argo_extended_svd_TA.mat
fit_sigma_cora_argo_extended_svd_TI.mat
fit_sigma_cora_argo_extended_svd_FC_ryanonly.mat
```
