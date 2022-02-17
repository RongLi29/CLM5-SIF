# CLM5-SIF
PhotosynthesisMod.F90: leaf level fluorescence and scaling up to canopy level (Kn calibrated with cotton dataset, from van der Tol et al, 2014)

PhotosynthesisMod_droughtfluorescencemodel.F90: Kn calibtrated with drought dataset (from van der Tol et al, 2014)

PhotosynthesisMod_sustainedNPQ.F90: consideres sustained NPQ according to Raczka et al (2019)

SurfaceAlbedoMod.F90: modifications of radiative transfer

SurfaceAlbedoType.F90: add variables and outputs related to radiative transfer

SurfaceRadiationMod.F90: add variables for PAR absorbed by leaves only (not affected by stem/snow); PAR/SR set to be 0.435 instead of 0.5. If coupled with CAM, the default-PAR version should be used.

SurfaceRadiationMod_oriPAR_SR_ratio.F90: The default PAR/SR=0.5 is used when not coupled with CAM (PAR is based on radiative transfer if coupled with CAM).

UrbanAlbedoMod.F90: Account for urban landuints for added variables (ftnn, ftin, refd, refi, refd_gr, refi_gr)

CanopyFluxesMod.F90: add input variables when calling PhotosynthesisTotal (for calculation of canopy scattering)

pftconMod.F90: add Clumping index for each pft as an input

Folder 'input' contains input parameter files with different modifications incorporated. All files have an added variable 'CI_pft' providing clumping index for each PFT. 

clm5_params.c171117_CI_Vcmax_Majrt.nc: with all modifications incorporated (added clumping index (CI), changed flnr (for TRY Vcmax), leaf reflecatnce and leaf transmittance)

clm5_params.c171117_CI1.nc: no modifications made (added variable 'CI_pft', with all clumping indices set to 1, i.e., no clumping)

clm5_params.c171117_CI.nc: considered clumping

clm5_params.c171117_CI_Majrt.nc: considered clumping, modified leaf reflectance and transmittance (based on Majasalmi et al, 2019)

clm5_params.c171117_CI_Vcmax.nc: considered clumping, modified flnr (for TRY Vcmax)

clm5_params.c171117_CI1_Vcmax_Majrt.nc: modified flnr (for TRY Vcmax) and leaf refleactance and transmittance, clumping not considered




