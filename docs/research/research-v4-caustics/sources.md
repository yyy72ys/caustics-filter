# Key sources read 2026-11-10 (index + opened pages)

- Guardado & Sanchez-Crespo, GPU Gems Ch.2 — opened live. Sun disc 0.53°, noon-equator assumption, backward single-ray method, useful rays 1-5%, water transparency claim, texture-space vs screen-space. https://developer.nvidia.com/gpugems/gpugems/part-i-natural-effects/chapter-2-rendering-water-caustics
- Terrill, Melville & Stramski 2001 JGR 106:16815 — opened live. Bubbles significant at void fraction >1e-6, b_bub >1e-2 m-1 in clouds, 5 decades variability in ~10 min at U10=15 m/s Hs=3.2m, radii O(10-1000) um, instrument void range 1e-8 to 1e-4, coated bubbles backscatter x4 (Zhang 1998 cited).
- Hieronymi & Macke — ResearchGate full text opened + Ocean Sci abstract. Peaks >10x mean at low wind <=5 m/s, >7x to 3m, 11x with ultra-gravity, 50% enhancements to 30-35m, +/-10% at 100m, Gernez 13x at 0.86m/532nm, diffuse light kills flashes (need >=30% direct), CV nearly independent of sun altitude (Gernez), Preisendorfer-Mobley surface reflectance depends on sun altitude.
- Veal et al 2010 / Stramski & Dera 1988 (Oceanologia PDF in search): Schenck 1957 up to ~6x mean, flashes 150% definition, 10-300ms, >350/min at 1m, short steep waves focus near surface, long flat deeper.
- Snell window: 97.5° full, half-angle 48.75° (MDPI 2023); 48.6° critical (Mobley-citing polarimetry paper). n_water 1.333.
- Fresnel table (illuminationpros, secondary but standard values): R 2.0% @0°, 2.1 @30, 2.8 @45, 6.0 @60, 21.2 @75.
- Wave-current: Gallet & Young 2014 Omega = sqrt(gk)+u.k; ECMWF/SWAN action N=E/sigma, action balance, opposing gradient -> shorter+steeper, following -> elongated+less steep, oblique -> refraction; USACE CETN: blocking U=-C/4 deep water, group velocity zero.
- Bubbles general: Koike/Flatau arXiv physics/0006060: b=2f/r_eff, Qsca->2, f=f0 W^4.4, f=1e-8 @6.2 m/s, 1e-7 @10.5 m/s at 2m, 1e-6 @0.3m to 1e-7 @2.7m, clouds to ~4Hs (some 6Hs), persist minutes. Stramski & Tegowski 2001 (EarthRef abstract): reflectance >2x in minutes at 10 m/s. Stramski 2024 slides: void 2.2e-5 @1.8m -> b_bub=1 m-1; concentration swings 3-4 orders.
- Deane & Stokes 2002 via UCSD: n(a) ~ a^-3/2 below Hinze ~0.5-1.5mm radius, a^-10/3 above.
- Foam: Koepke 1984 effective reflectance ~22% visible (14.6% in another averaging), Frouin 1996 foam reflectance spectral, ~4% coverage at ~14 m/s; Evans & Stride: foam albedo 0.4-0.6 fresh vs ocean ~0.05.
- Rain: Laxague 2020 GRL + Bliven 1993: rain damps gravity waves >few cm, generates ring waves 1-2cm wavelength, visible patch 25-50cm diameter, craters/crowns/stalks, rain rates 3-30 mm/h in tank studies.
- Cox-Munk: sigma^2 = 0.003+5.12e-3 U (total), upwind 0.00316U+..., crosswind 0.00192U (IOCCG slides).
- Photon mapping: Jensen 1996/2001 separate caustic map LS+D, density estimation; Yuksel & Keyser heightfield two-pass; Shah et al caustics mapping iterates vs light-space depth -> works over bathymetry (vs plane-only methods); Guardado single-point = noon only, highly inaccurate (Yuksel critique).
- Bottom albedo: Mobley shallow-water report: sand material reflectance 0.2-0.7, multiple scattering between ripples cannot be ignored; Mishra 2005 retrieved albedos seagrass 12%, coral 12-24%, sand 24%; Fogarty 2018 system albedo up to 0.20 vs 0.05 param over bright sand / turbid >=50 g/m3.
- Tessendorf 2001: Phillips Ph(k)=A |khat.what|^2 /k^4 exp(-1/(kL)^2), L=V^2/g, dispersion sqrt(gk).
- Premoze & Ashikhmin 2000 Rendering Natural Waters: full Fresnel, beam tracing, c=a+bb, caustic image extension.
