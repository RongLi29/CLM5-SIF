# CLM5-SIF
PhotosynthesisMod.F90: leaf level fluorescence and scaling up to canopy level 

PhotosynthesisMod_droughtfluorescencemodel.F90: Kn calibtrated with drought dataset

PhotosynthesisMod_sustainedNPQ.F90: consideres sustained NPQ according to Raczka et al (2019)

SurfaceAlbedoMod.F90: modifications of radiative transfer

SurfaceAlbedoMod_sigma.F90: modifications of radiative transfer (the singularity of sigma considered)

SurfaceAlbedoType.F90: add variables and outputs related to radiative transfer

SurfaceRadiationMod.F90: add variables for PAR absorbed by leaves only (not affected by stem/snow); PAR/SR set to be 0.435 instead of 0.5. If coupled with CAM, the default-PAR version should be used.

UrbanAlbedoMod.F90: Account for urban landuints for added variables (ftnn, ftin, refd, refi, refd_gr, refi_gr)

(not necessary) SolarAbsorbedType.F90: add outputs(PARSUNZ PARSHAZ) 

CanopyFluxesMod.F90: add input variables when calling PhotosynthesisTotal (for calculation of canopy scattering)

pftconMod.F90: add Clumping index for each pft as an input

clm5_params.c171117_CI_Vcmax_rt.nc: add clumping index (CI), change flnr (for TRY Vcmax), leaf reflecatnce & transmittance
