GOODS-S
-------------

RELEASE: V4.1

DATE: MARCH 15TH, 2014

SOURCE: http://3dhst.research.yale.edu

Downloads
^^^^^^^^^^^^^^

CATALOG DOWNLOAD: goodss_3dhst.v4.1.cats.tar

This file opens into a directory: goodss_3dhst.v4.1.cats. The following are the contents of this directory:

* Catalog/

	* goodss_3dhst.v4.1.cat: photometric catalog, ASCII format
	* goodss_3dhst.v4.1.cat.FITS: photometric catalog, FITS format
	
* Eazy/

	* goodss_3dhst.v4.1.coeff: 
Coefficients of fits for all objects
	* goodss_3dhst.v4.1.param: 
Parameter file containing all of the parameter values used to compute goodss_3dhst.v4.1.zout, as well as information on the individual filters and templates used.

	* goodss_3dhst.v4.1.pz: Binary output file containing redshift z, chi squared at z, prior, probability distribution pz

	* goodss_3dhst.v4.1.readbin.pro: IDL code demonstrating how to read in the binary output from EAZY

	* goodss_3dhst.v4.1.temp_sed
	* goodss_3dhst.v4.1.tempfilt
	* goodss_3dhst.v4.1.translate: Eazy translate file listing the flux and error columns and corresponding filter number in filter transmission file FILTERS.res.latest
	* goodss_3dhst.v4.1.zbin: Binary output file containing observed and best-fit template for each object

	* goodss_3dhst.v4.1.zout: EAZY output file, ASCII format
	* goodss_3dhst.v4.1.zout.FITS: EAZY output file, FITS format


* Fast/

	* goodss_3dhst.v4.1.fout: FAST output file, ASCII format
	* goodss_3dhst.v4.1.fout.FITS: FAST output file, FITS format
	* goodss_3dhst.v4.1.param: input parameter file for FAST
	* goodss_3dhst.v4.1.translate: filter translate file for FAST

* Restframe/
	
	* goodss_3dhst.v4.1.master.RF: rest frame colors, ASCII format
	* goodss_3dhst.v4.1.master.RF.FITS: rest frame colors, fits format



The contents of each file are explained below.

Photometric Catalog
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The detailed methods for creating the catalogs are described in a companion paper - Skelton et al. (2014). 
The photometric catalog for the goodss field contains the following datasets:


+----------------------------+------------+------------------------------------------+
|                            |            |                                          |
| Bands                      |  Survey    |  Reference                               |
|                            |            |                                          |
+============================+============+==========================================+
| U, R                       | ESO GOODS  |  Nonino et al. 2009                      |
+----------------------------+------------+------------------------------------------+
| U38, B, V, Rc, I           | GaBoDs     | Hildebrandt et al. 2006, Erben 2005      |
+----------------------------+------------+------------------------------------------+
| 14 medium band filters     | MUSYC      | Cardamone et al. 2010                    |
+----------------------------+------------+------------------------------------------+
|F435W,  F606W, F775W, F850LP| GOODS      |Giavalisco et al. 2004                    |
+----------------------------+------------+------------------------------------------+
| F606W,F814W                | CANDELS    | Koekemoer et al. 2011                    |
+----------------------------+------------+------------------------------------------+
| J, H, Ks                   | ESO/GOODS  | Retzlaff et al. 2010, Wuyts et al. 2008  |
+----------------------------+------------+------------------------------------------+
| J, Ks                      | TENIS      | Hsieh et al. 2012                        |
+----------------------------+------------+------------------------------------------+
| F140W                      | 3D-HST     | Brammer et al. 2012                      |
+----------------------------+------------+------------------------------------------+
| F125W,F160W                | CANDELS    | Grogin et al. 2011, Koekemoer et al. 2011|
+----------------------------+------------+------------------------------------------+
| 3.6,4.5um                  | SEDS       | Ashby et al. (2013)                      |
+----------------------------+------------+------------------------------------------+
| 5.8,8.0um                  | GOODS      | Dickinson et al. 2003                    |
+----------------------------+------------+------------------------------------------+

The catalog has a single line header with all column names as shown here. The ASCII and the 
FITS versions of the catalog contain the same columns.

::

	#id x y ra dec faper_F160W eaper_F160W faper_F140W eaper_F140W f_F160W e_F160W w_F160W 
	f_U38 e_U38 w_U38 f_U e_U w_U f_F435W e_F435W w_F435W f_B e_B w_B f_V e_V w_V f_F606Wcand 
	e_F606Wcand w_F606Wcand f_F606W e_F606W w_F606W f_R e_R w_R f_Rc e_Rc w_Rc f_F775W e_F775W 
	w_F775W f_I e_I w_I f_F814Wcand e_F814Wcand w_F814Wcand f_F850LP e_F850LP w_F850LP f_F850LPcand 
	e_F850LPcand w_F850LPcand f_F125W e_F125W w_F125W f_J e_J w_J f_tenisJ e_tenisJ w_tenisJ 
	f_F140W e_F140W w_F140W f_H e_H w_H f_tenisK e_tenisK w_tenisK f_Ks e_Ks w_Ks f_IRAC1 
	e_IRAC1 w_IRAC1 f_IRAC2 e_IRAC2 w_IRAC2 f_IRAC3 e_IRAC3 w_IRAC3 f_IRAC4 e_IRAC4 w_IRAC4 
	f_IA427 e_IA427 f_IA445 e_IA445 f_IA505 e_IA505 f_IA527 e_IA527 f_IA550 e_IA550 f_IA574 
	e_IA574 f_IA598 e_IA598 f_IA624 e_IA624 f_IA651 e_IA651 f_IA679 e_IA679 f_IA738 e_IA738 
	f_IA767 e_IA767 f_IA797 e_IA797 f_IA856 e_IA856 tot_cor wmin_ground wmin_hst wmin_irac 
	wmin_wfc3 z_spec star_flag  kron_radius a_image b_image theta_J2000 class_star flux_radius 
	fwhm_image flags IRAC1_contam IRAC2_contam IRAC3_contam IRAC4_contam contam_flag f140w_flag 
	use_phot near_star nexp_f125w nexp_f140w nexp_f160w
	# goodss_3dhst.v4.1.cat:
	# by R. E. Skelton (3/11/2014)

All fluxes are normalized to an AB zeropoint of 25, such that: magAB = 25.0-2.5*log10(flux).

=====================  =======================================================================================
Column Name               Column Content
=====================  =======================================================================================
1    # id				Unique identifier within a given field  	
2    # x				X centroid in image coordinates
3    # y				Y centroid in image coordinates  
4    # ra				RA J2000 degrees 
5    # dec				Dec J2000 degrees   
6    # faper_F160W		Flux within a 0.7 arcsecond aperture in F160W zeropoint 25.0 	  
7    # eaper_F160W		1 sigma error within a 0.7 arcsecond aperture in F160W zeropoint 25.0 
8    # faper_F140W		Flux within a 0.7 arcsecond aperture in F140W zeropoint 25.0 		  
9    # eaper_F140W		1 sigma error within a 0.7 arcsecond aperture in F140W zeropoint 25.0 		  
10   # f_F160W			Total flux in F160W zeropoint 25 	  
11   # e_F160W			1 sigma error in F160W total flux zeropoint 25 	  
12   # w_F160W			Weight relative to 95th percentile within F160W weight map
13   # f_U38            Total flux in U38  zeropoint 25 
14   # e_U38			1 sigma error in U38 total flux zeropoint 25 
15   # w_U38            Weight relative to 95th percentile within U38 weight map						  
16   # f_U				Total flux in U  zeropoint 25 
17   # e_U				1 sigma error in U total flux zeropoint 25 
18   # w_U				Weight relative to 95th percentile within U weight map	
19   # f_F435W			Total flux in ACS F435W  zeropoint 25 
20   # e_F435W			1 sigma error in ACS F435W total flux zeropoint 25 
21   # w_F435W			Weight relative to 95th percentile within F435W weight map	
22   # f_B              Total flux in B zeropoint 25 						  
23   # e_B				1 sigma error in B total flux zeropoint 25 				
24   # w_B              Weight relative to 95th percentile within B weight map						  
25   # f_V              Total flux in V zeropoint 25 						  
26   # e_V				1 sigma error in V total flux zeropoint 25 						  
27   # w_V              Weight relative to 95th percentile within V weight map						  
28   # f_F606Wcand      Total flux in F606W CANDELS zeropoint 25 						  
29   # e_F606Wcand		1 sigma error in F606W CANDELS total flux zeropoint 25 						  
30   # w_F606Wcand     	Weight relative to 95th percentile within F606W CANDELS weight map					
31   # f_F606W         	Total flux in F606W GOODS zeropoint 25 						  
32   # e_F606W			1 sigma error in F606W GOODS total flux zeropoint 25 						  
33   # w_F606W          Weight relative to 95th percentile within F606W GOODS weight map						  
34   # f_R              Total flux in R zeropoint 25 						  
35   # e_R				1 sigma error in R total flux zeropoint 25 						  
36   # w_R              Weight relative to 95th percentile within R weight map						  
37   # f_Rc             Total flux in Rc zeropoint 25 						  
38   # e_Rc				1 sigma error in Rc total flux zeropoint 25 						  
39   # w_Rc             Weight relative to 95th percentile within Rc weight map						  
40   # f_F775W          Total flux in F775W GOODS zeropoint 25 						  
41   # e_F775W			1 sigma error in F775W GOODS total flux zeropoint 25 						  
42   # w_F775W          Weight relative to 95th percentile within F775W GOODS weight map						  
43   # f_I              Total flux in I zeropoint 25 						  
44   # e_I				1 sigma error in I total flux zeropoint 25 						  
45   # w_I              Weight relative to 95th percentile within I weight map						  
46   # f_F814Wcand      Total flux in F814W CANDELS zeropoint 25 
47   # e_F814Wcand		1 sigma error in F814W CANDELS total flux zeropoint 25 
48   # w_F814Wcand      Weight relative to 95th percentile within F814W CANDELS weight map 
49   # f_F850LP         Total flux in F850LP GOODS zeropoint 25 
50   # e_F850LP			1 sigma error in F850LP GOODS total flux zeropoint 25 
51   # w_F850LP         Weight relative to 95th percentile within F850LP GOODS weight map
52   # f_F850LPcand     Total flux in F850LP CANDELS zeropoint 25 
53   # e_F850LPcand		1 sigma error in F850LP CANDELS total flux zeropoint 25 
54   # w_F850LPcand     Weight relative to 95th percentile within F850LP CANDELS weight map
55   # f_F125W          Total flux in F125W  zeropoint 25 
56   # e_F125W			1 sigma error in F125W  total flux zeropoint 25 
57   # w_F125W          Weight relative to 95th percentile within F125W weight map
58   # f_J              Total flux in Jzeropoint 25 
59   # e_J				1 sigma error in J total flux zeropoint 25 
60   # w_J              Weight relative to 95th percentile within J weight map
61   # f_tenisJ         Total flux in TENIS J zeropoint 25 
62   # e_tenisJ			1 sigma error in TENIS J total flux zeropoint 25 
63   # w_tenisJ         Weight relative to 95th percentile within TENIS J weight map
64   # f_F140W          Total flux in zeropoint 25 
65   # e_F140W			1 sigma error in  total flux zeropoint 25 
66   # w_F140W          Weight relative to 95th percentile within weight map
67   # f_H              Total flux in H zeropoint 25 
68   # e_H				1 sigma error in H total flux zeropoint 25 
69   # w_H              Weight relative to 95th percentile within H weight map
70   # f_tenisK         Total flux in TENIS K zeropoint 25 
71   # e_tenisK			1 sigma error in TENIS K total flux zeropoint 25 
72   # w_tenisK         Weight relative to 95th percentile within TENIS K weight map
73   # f_Ks             Total flux in K zeropoint 25 
74   # e_Ks				1 sigma error in K total flux zeropoint 25 
75   # w_Ks             Weight relative to 95th percentile within K weight map
76   # f_IRAC1          Total flux in IRAC1 zeropoint 25 
77   # e_IRAC1			1 sigma error in IRAC1 total flux zeropoint 25 
78   # w_IRAC1          Weight relative to 95th percentile within IRAC1 weight map
79   # f_IRAC2          Total flux in IRAC2 zeropoint 25 
80   # e_IRAC2			1 sigma error in IRAC2 total flux zeropoint 25 
81   # w_IRAC2          Weight relative to 95th percentile within IRAC2 weight map
82   # f_IRAC3          Total flux in IRAC3 zeropoint 25 
83   # e_IRAC3			1 sigma error in IRAC3 total flux zeropoint 25 
84   # w_IRAC3          Weight relative to 95th percentile within IRAC3 weight map
85   # f_IRAC4          Total flux in IRAC4 zeropoint 25 
86   # e_IRAC4			1 sigma error in IRAC4 total flux zeropoint 25 
87   # w_IRAC4          Weight relative to 95th percentile within IRAC4 weight map
88   # f_IA427 			Total flux in IA427  zeropoint 25 
89   # e_IA427 			1 sigma error in IA427 total flux zeropoint 25 
90   # f_IA445 			Total flux in IA445  zeropoint 25 
91   # e_IA445 			1 sigma error in IA445 total flux zeropoint 25 
92   # f_IA505 			Total flux in IA505  zeropoint 25 
93   # e_IA505 			1 sigma error in IA505 total flux zeropoint 25 
94   # f_IA527 			Total flux in IA527  zeropoint 25 
95   # e_IA527 			1 sigma error in IA527 total flux zeropoint 25 
96   # f_IA550 			Total flux in IA550  zeropoint 25 
97   # e_IA550 			1 sigma error in IA550 total flux zeropoint 25 
98   # f_IA574 			Total flux in IA574  zeropoint 25 
99   # e_IA574 			1 sigma error in IA574 total flux zeropoint 25 
100  # f_IA598 			Total flux in IA598  zeropoint 25 
101  # e_IA598 			1 sigma error in IA598 total flux zeropoint 25 
102  # f_IA624 			Total flux in IA624  zeropoint 25 
103  # e_IA624 			1 sigma error in IA624 total flux zeropoint 25 
104  # f_IA651 			Total flux in IA651  zeropoint 25 
105  # e_IA651 			1 sigma error in IA651 total flux zeropoint 25 
106  # f_IA679 			Total flux in IA679  zeropoint 25 
107  # e_IA679 			1 sigma error in IA679 total flux zeropoint 25 
108  # f_IA738 			Total flux in IA738  zeropoint 25 
109  # e_IA738 			1 sigma error in IA738 total flux zeropoint 25 
110  # f_IA767 			Total flux in IA767  zeropoint 25 
111  # e_IA767 			1 sigma error in IA767 total flux zeropoint 25 
112  # f_IA797 			Total flux in IA797  zeropoint 25 
113  # e_IA797 			1 sigma error in IA797 total flux zeropoint 25 
114  # f_IA856 			Total flux in IA856  zeropoint 25 
115  # e_IA856			1 sigma error in IA856 total flux zeropoint 25 
116  # tot_cor          Correction from AUTO to total flux based on F160W F140W 
117  # wmin_ground      Minimum weight for all ground-based photometry
118  # wmin_hst         Minimum weight for ACS/WFC3 bands excluding no coverage 
119  # wmin_wfc3        Minimum weight for F160W, F125W and F140W excluding no coverage 
120  # wmin_irac		Minimum weight for IRAC bands excluding no coverage 
121  # z_spec           Spectroscopic redshift, when available see Skelton et al., 2014 for sources 
122  # star_flag        For F160W<25, star=1 and galaxy=0; for F160W>25, flag=2
123  # kron_radius      KRON_RADIUS 
124  # a_image			A_IMAGE semi-major axis, pixels 
125  # b_image			B_IMAGE semi-minor axis, pixels 
126  # theta_J2000		Position angle of the major axis counter-clockwise, 0.0 = X world axis 
127  # class_star       Sextractor stellarity-index CLASS_STAR 
128  # flux_radius      Circular aperture radius enclosing half the total flux
129  # fwhm_image       FWHM pixels  from a gaussian fit to the core
130  # flags			SExtractor extraction flags measured
131  # IRAC1_contam		Ratio of contaminating flux from neighbors to the object’s flux IRAC CH1
132  # IRAC2_contam		Ratio of contaminating flux from neighbors to the object’s flux IRAC CH2
133  # IRAC3_contam		Ratio of contaminating flux from neighbors to the object’s flux IRAC CH3
134  # IRAC4_contam		Ratio of contaminating flux from neighbors to the object’s flux IRAC CH4
135  # contam_flag		A flag for IRAC phot contam ratio >50% all bands . 0=OK, 1=bad
136  # f140w_flag		Set if F140W is used for the corrections to total i.e., no F160W coverage 
137  # use_phot			Use flag: 1=USE, 0=DON'T USE; see Skelton et al. 2014  for definition
138  # near_star		1=Close to a star
139  # nexp_f125w		Number of individual exposures in F125W, based on NEXP maps
140  # nexp_f140w		Number of individual exposures in F140W, based on NEXP maps
141  # nexp_f160w		Number of individual exposures in F160W, based on NEXP maps
=====================  =======================================================================================

Photometric Redshifts
-----------------------

See the EAZY documentation (http://www.astro.yale.edu/eazy and Brammer, van Dokkum & Coppi 2008) for more information on the EAZY output.

The output catalogs are in the *_3dhst.v4.1.zout.

The first line of the catalog is a header with all the column names as shown here.

::

	# id z_spec z_a z_m1 chi_a z_p chi_p z_m2 odds l68 u68  l95 u95  l99 u99  nfilt q_z 
		z_peak peak_prob z_mc
	# EAZY $Rev: 34 $

=================  ==================================================================================================
Column Name        Column Content
=================  ==================================================================================================
1    # id			Unique identifier within a given field
2    # z_spec 		Spectroscopic redshift from photometric catalog, if it exists; If no zspec = -1
3    # z_a 			Redshift where chi2 is minimized for all template linear combination modes before applying prior
4    # z_m1 			Marginalized Redshift without prior
5    # chi_a 		Minimum chi squared for all template linear combinations
6    # z_p 			Redshift where likelihood is maximized after applying prior
7    # chi_p 		Original chi squared at z_p
8    # z_m2 			Marginalized redshift including prior
9    # odds 			Redshift quality parameter from Benitez 2000   - fraction of probability within 0.2 of zphot
10   # l68 			lower limit of 1 sigma confidence interval
11   # u68  		upper limit of 1 sigma confidence interval
12   # l95 			lower limit of 2 sigma confidence interval
13   # u95  		upper limit of 2 sigma confidence interval
14   # l99 			lower limit of 3 sigma confidence interval
15   # u99  		upper limit of 3 sigma confidence interval
16   # nfilt 		Number of filters used in fit
17   # q_z 			Redshift quality parameter; see Brammer et al. 2008  
18   # z_peak 		Photometric redshift at peak probability distribution
19   # peak_prob 	Value of peak probability
20   # z_mc			Monte Carlo photometric redshift drawn randomly from probability distribution
=================  ==================================================================================================

Stellar Population Fits
------------------------

Output from FAST (see Kriek et al. 2009 for details)

The header of the .fout file is as follows:

::
	
	#    id         z      ltau     metal      lage        Av     lmass      lsfr     
		lssfr      la2t      chi2
	# FAST version: 0.9b
	# Photometric catalog file: uds_3dhst.v4.1.cat
	# Photometric redshift file: uds_3dhst.v4.1.zout
	# Template error function: TEMPLATE_ERROR.fast.v0.2
	# AB ZP:        25.00
	# Library:      Bruzual & Charlot (2003)
	# SFH:          Exponentially declining SFH: SFR ~ exp(-t/tau)
	# Stellar IMF:  Chabrier
	# metallicity:  0.020
	# log(tau/yr):  7.0    - 10.0, in steps of 0.20
	# log(age/yr):  7.6    - 10.1, in steps of 0.10
	# A_V:          0.0    -  4.0, in steps of 0.10
	# z:            0.0100 -  6.0000, in steps of 0.0100
	# Filters:     205  88 122  79 236 123 124 239 125 203 263 204 264 265  18  19  20  21
	# ltau: log[tau/yr], lage: log[age/yr], lmass: log[mass/Msol], lsfr: log[sfr/(Msol/yr)], 
		lssfr: log[ssfr*yr], la2t: log[age/tau]
	# For sfr=0. lsfr is set to -99	
	#    id         z      ltau     metal      lage        Av     lmass      lsfr     lssfr  
			la2t      chi2

================  ==========================================================================================
Column Name       Column Content
================  ==========================================================================================
1   # id 			Unique identifier within a given field
2   # z				Redshift used in fit z_spec if exists in input catalog, z_peak from eazy otherwise  
3   # ltau			logtau/yr  
4   # metal			metallicity 
5   # lage			logage/yr  
6   # Av     		Dust attenuation in the V-band
7   # lmass      	log Mstar/Msun  
8   # lsfr     		log SFR/Msun/yr    
9   # lssfr      	log SSFR * yr  
10   # la2t      	log age/tau  
11   # chi2			chi squared of fit
================  ==========================================================================================


Rest Frame Colors
--------------------

The first line of the catalog *_3dhst.v4.1.master.RF is a header with all column names as shown here.

::
	
	# id z DM L153 n_153 L154 n_154 L155 n_155 L161 n_161 L162 n_162 L163 n_163 L156 n_156 L157 n_157 L158 n_158 L159 n_159 L160 n_160 L135 n_135 L136 n_136 L137 n_137 L138 n_1	
	#
	# 153: REST_FRAME/maiz-apellaniz_Johnson_U.res,   3.59854e+03
	# 154: REST_FRAME/maiz-apellaniz_Johnson_B.res,   4.38592e+03
	# 155: REST_FRAME/maiz-apellaniz_Johnson_V.res,   5.49056e+03
	# 161: 2MASS/J.res,   1.23751e+04
	# 162: 2MASS/H.res,   1.64763e+04
	# 163: 2MASS/K.res,   2.16203e+04
	# 156: SDSS/u.dat,   3.56179e+03
	# 157: SDSS/g.dat,   4.71887e+03
	# 158: SDSS/r.dat,   6.18519e+03
	# 159: SDSS/i.dat,   7.49966e+03
	# 160: SDSS/z.dat,   8.96122e+03
	# 135: REST_FRAME/Bessel_UX.dat,   3.59291e+03
	# 136: REST_FRAME/Bessel_B.dat,   4.38477e+03
	# 137: REST_FRAME/Bessel_V.dat,   5.48882e+03
	# 138: REST_FRAME/Bessel_R.dat,   6.48893e+03
	# 139: REST_FRAME/Bessel_I.dat,   8.03337e+03
	# 270: RestUV/Tophat_1400_200.dat,   1.39971e+03
	# 271: RestUV/Tophat_1700_200.dat,   1.69989e+03
	# 272: RestUV/Tophat_2200_200.dat,   2.20011e+03
	# 273: RestUV/Tophat_2700_200.dat,   2.70000e+03
	# 274: RestUV/Tophat_2800_200.dat,   2.80016e+03
	# 275: RestUV/Tophat_5500_200.dat,   5.50055e+03	
	#	
	# z = z_peak / z_spec	
	#

The filters are numbered in the EAZY format and listed in *_3dhst.v4.1.master.RF. 

All fluxes are normalized to an AB zeropoint of 25, such that: magAB = 25.0-2.5*log10(flux). 

Rest-frame colors are calculated as color = -2.5*log(Lfilter1/Lfilter2)

================  ====================================================================================================
Column Name       Column Content
================  ====================================================================================================
1   # id			Unique identifier within a given field
2   # z				Redshift z_spec if it exists, z_peak from EAZY otherwise   
3   # DM 			Distance modulus for a W_M/W_L/H0 = 0.3/0.7/70 cosmology
4   # L153			Fnu flux density for filter REST_FRAME/maiz-apellaniz_Johnson_U.res with an AB zeropoint of 25 
5   # n_153 			Number of filters considered in RF fit for L153 
6   # L154 				Fnu flux density for filter REST_FRAME/maiz-apellaniz_Johnson_B.res with an AB zeropoint of 25 
7   # n_154 			Number of filters considered in RF fit for L154
8   # L155 				Fnu flux density for filter RREST_FRAME/maiz-apellaniz_Johnson_V.res with an AB zeropoint of 25 
9   # n_155 			Number of filters considered in RF fit for L155
10   # L161 				Fnu flux density for filter 2MASS/J.res with an AB zeropoint of 25 
11   # n_161 			Number of filters considered in RF fit for L161
12   # L162 				Fnu flux density for filter 2MASS/H.res with an AB zeropoint of 25 
13   # n_162 			Number of filters considered in RF fit for L162
14   # L163 				Fnu flux density for filter 2MASS/K.res with an AB zeropoint of 25 
15   # n_163 			Number of filters considered in RF fit for L163
16   # L156 				Fnu flux density for filter SDSS/u.dat with an AB zeropoint of 25 
17   # n_156 			Number of filters considered in RF fit for L156
18   # L157 				Fnu flux density for filter SDSS/g.da with an AB zeropoint of 25 
19   # n_157 			Number of filters considered in RF fit for L157
20   # L158 				Fnu flux density for filter SDSS/r.dat with an AB zeropoint of 25 
21   # n_158 			Number of filters considered in RF fit for L158
22   # L159 				Fnu flux density for filter SDSS/i.dat with an AB zeropoint of 25 
23   # n_159 			Number of filters considered in RF fit for L159
24   # L160 				Fnu flux density for filter SDSS/z.dat with an AB zeropoint of 25 
25   # n_160 			Number of filters considered in RF fit for L160
26   # L135 				Fnu flux density for filter REST_FRAME/Bessel_UX.dat with an AB zeropoint of 25 
27   # n_135 			Number of filters considered in RF fit for L135
28   # L136 				Fnu flux density for filter REST_FRAME/Bessel_B.dat with an AB zeropoint of 25 
29   # n_136 			Number of filters considered in RF fit for L136
30   # L137 				Fnu flux density for filter REST_FRAME/Bessel_V.da with an AB zeropoint of 25 
31   # n_137 			Number of filters considered in RF fit for L137
32   # L138 				Fnu flux density for filter REST_FRAME/Bessel_R.dat with an AB zeropoint of 25 
33   # n_138 			Number of filters considered in RF fit for L138
34   # L139 				Fnu flux density for filter REST_FRAME/Bessel_I.dat with an AB zeropoint of 25 
35   # n_139 			Number of filters considered in RF fit for L139
36   # L270 				Fnu flux density for filter RestUV/Tophat_1400_200.dat with an AB zeropoint of 25 
37   # n_270 			Number of filters considered in RF fit for L270
38   # L271 				Fnu flux density for filter RestUV/Tophat_1700_200.dat with an AB zeropoint of 25 
39   # n_271 			Number of filters considered in RF fit for L271
40   # L272 				Fnu flux density for filter RestUV/Tophat_2200_200.dat with an AB zeropoint of 25 
41   # n_272 			Number of filters considered in RF fit for L272
42   # L273 				Fnu flux density for filter RestUV/Tophat_2700_200.dat with an AB zeropoint of 25 
43   # n_273 			Number of filters considered in RF fit for L273
44   # L274 				Fnu flux density for filter RestUV/Tophat_2800_200.dat with an AB zeropoint of 25 
45   # n_274 			Number of filters considered in RF fit for L274
46   # L275 				Fnu flux density for filter RestUV/Tophat_5500_200.dat with an AB zeropoint of 25 
47   # n_275  			Number of filters considered in RF fit for L275
================  ====================================================================================================
