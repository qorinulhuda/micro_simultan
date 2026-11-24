********************************************************************************
* PROJECT : SIMULASI DATA MIKRO-MAKRO & SISTEM PERSAMAAN SIMULTAN (SEM)
*           Eksperimen Monte Carlo Data Sintetis Individu ke Agregat Kab/Kota
* METODE  : Weighted Collapse, 2SLS, 3SLS, IV-GMM, SUR, Quantile, & Survey-IV
* AUTHOR  : Qorinul Huda
********************************************************************************

clear all
set more off
set seed 20261124

** =============================================================================
** TAHAP 1: PEMBANGKITAN DATASET MIKRO BESAR (N = 50.000 INDIVIDU, 100 KAB/KOTA)
** =============================================================================
	* 1. Membuat kerangka 100 Kabupaten/Kota sintetis dengan variabel makro & instrumen
	set obs 100
	gen int id_kab = 1100 + _n
	gen byte pulau = ceil(_n / 20)
	la def pulau 1 "Wilayah 1" 2 "Wilayah 2" 3 "Wilayah 3" 4 "Wilayah 4" 5 "Wilayah 5"
	la val pulau pulau
	
	* Variabel Eksogen & Instrumen Tingkat Kabupaten/Kota:
	* - infra_kab : Indeks Infrastruktur Daerah (instrumen sisi penawaran/upah)
	* - rasio_sek : Rasio Ketersediaan Sekolah Daerah (instrumen pendidikan/konsumsi)
	* - umk_riil  : Upah Minimum Kabupaten/Kota Riil (instrumen kebijakan upah)
	gen double infra_kab    = rnormal(50, 15) + 5*(pulau==2)
	gen double rasio_sek    = rnormal(65, 12) + 4*(pulau==2)
	gen double umk_riil     = 2.2 + 0.02*infra_kab + rnormal(0, 0.3)
	
	* Heterogenitas tak teramati (unobserved heterogeneity) tingkat kabupaten/kota
	gen double u_kab_earn   = rnormal(0, 0.25)
	gen double u_kab_educ   = rnormal(0, 0.40)
	
	* 2. Ekspansi menjadi 500 observasi individu per kabupaten/kota (Total N = 50.000)
	expand 500
	bysort id_kab: gen long id_ind = _n
	order id_kab id_ind pulau
	
	* Penimbang survei (final_weight): merepresentasikan jumlah populasi per responden
	gen double final_weight = runiform(80, 450)
	
	* Karakteristik Demografi Individu
	gen byte female   = (runiform() < 0.49)
	gen int age       = floor(runiform(18, 65))
	gen double exp    = max(0, age - 15)
	gen double exp_sq = exp^2
	
	* Kemampuan bawaan tak teramati (unobserved ability) yang memicu bias endogenitas
	* karena berkorelasi dengan lama sekolah (educ_yr) sekaligus upah (ln_wage)
	gen double ability = rnormal(0, 1)
	
	* 3. Data Generating Process (DGP) Level Individu
	* Lama sekolah individu dipengaruhi instrumen daerah (rasio_sek), demografi, & ability
	gen double educ_yr = 7.5 + 0.06*rasio_sek - 0.35*female - 0.03*(age-35) + 1.2*ability + u_kab_educ + rnormal(0, 1.8)
	replace educ_yr = round(max(0, min(20, educ_yr)), 1)
	
	* Log Upah Individu (Mincer Earnings Function + Infrastruktur Daerah + Ability)
	gen double ln_wage = 0.45 + 0.085*educ_yr + 0.032*exp - 0.0005*exp_sq - 0.18*female + 0.012*infra_kab + 0.25*ability + u_kab_earn + rnormal(0, 0.35)
	
	* Status Pekerja Formal (Variabel Biner 0/1 untuk model non-linear Probit/Logit)
	gen double latent_formal = -2.2 + 0.18*educ_yr + 0.015*infra_kab - 0.25*female + 0.30*ability + rnormal(0, 1)
	gen byte formal = (latent_formal > 0)
	
	* Log Pengeluaran Konsumsi Per Kapita Individu
	gen double ln_cons = 0.30 + 0.72*ln_wage + 0.015*educ_yr - 0.004*age + rnormal(0, 0.22)
	
	* Hapus variabel laten tak teramati agar menyerupai struktur data survei
	drop latent_formal ability u_kab_earn u_kab_educ
	
	* Pelabelan Variabel Dataset Mikro Sintetis
	la var id_kab       "Kode Kabupaten/Kota (Sintetis)"
	la var id_ind       "ID Unik Individu"
	la var final_weight "Penimbang Survei (Survey Weight)"
	la var female       "Perempuan (1=Ya, 0=Tidak)"
	la var age          "Usia Individu (Tahun)"
	la var educ_yr      "Lama Sekolah Individu (Tahun)"
	la var ln_wage      "Log Pendapatan/Upah Bulanan"
	la var ln_cons      "Log Konsumsi Per Kapita"
	la var formal       "Status Pekerja Formal (1=Formal, 0=Informal)"
	la var infra_kab    "Indeks Infrastruktur Kabupaten/Kota (Instrumen Z1)"
	la var rasio_sek    "Rasio Ketersediaan Sekolah Daerah (Instrumen Z2)"
	la var umk_riil     "Upah Minimum Kabupaten Riil (Instrumen Z3)"
	
	compress
	save "data_mikro_simulasi_50k.dta", replace
	describe


** =============================================================================
** TAHAP 2: AGREGASI TINGKAT KABUPATEN/KOTA (WEIGHTED COLLAPSE & INFERENSI)
** =============================================================================
	use "data_mikro_simulasi_50k.dta", clear
	
	* CATATAN METODOLOGIS (WEIGHTED AGGREGATION):
	* Untuk menghasilkan indikator makro kabupaten/kota yang representatif terhadap
	* populasi, proses collapse menggunakan penimbang survei [aw = final_weight].
	* Total penimbang (est_populasi) disimpan sebagai bobot inferensi pada tahap makro.
	
	preserve
		collapse (count) n_sampel = id_ind (sum) est_populasi = final_weight, by(id_kab)
		tempfile pop_kab
		save `pop_kab'
	restore
	
	collapse (mean) avg_wage = ln_wage avg_cons = ln_cons rls_kab = educ_yr share_formal = formal share_female = female avg_age = age (first) pulau infra_kab rasio_sek umk_riil [aw = final_weight], by(id_kab)
			 
	merge 1:1 id_kab using `pop_kab', nogen
	
	* Simulasi Umpan Balik Simultan (Two-Way Causality / Simultaneity Feedback):
	* Di tingkat regional, kenaikan konsumsi agregat (avg_cons) mendorong permintaan
	* lokal & menaikkan upah (avg_wage), sementara kenaikan upah (avg_wage) meningkatkan
	* daya beli dan konsumsi (avg_cons). Ini menciptakan bias simultanitas pada OLS.
	gen double shock_macro = rnormal(0, 0.08)
	replace avg_wage = avg_wage + 0.35*avg_cons + shock_macro
	replace avg_cons = avg_cons + 0.45*avg_wage - 0.5*shock_macro
	
	* Konversi proporsi ke satuan persentase (0 - 100%)
	replace share_formal = share_formal * 100
	replace share_female = share_female * 100
	
	la var est_populasi "Estimasi Total Populasi Kab/Kota (Bobot Inferensi Makro)"
	la var avg_wage     "Rata-rata Log Upah Kab/Kota (Endogen Y1)"
	la var avg_cons     "Rata-rata Log Konsumsi Kab/Kota (Endogen Y2)"
	la var rls_kab      "Rata-rata Lama Sekolah Kab/Kota (Tahun)"
	la var share_formal "Pangsa Pekerja Formal Kab/Kota (%)"
	la var infra_kab    "Indeks Infrastruktur Kab/Kota (Instrumen Z1)"
	la var rasio_sek    "Rasio Ketersediaan Sekolah Kab/Kota (Instrumen Z2)"
	la var umk_riil     "UMK Riil Kab/Kota (Instrumen Z3)"
	
	order id_kab pulau n_sampel est_populasi avg_wage avg_cons rls_kab share_formal infra_kab rasio_sek umk_riil
	save "data_agregat_simulasi_kabkota.dta", replace
	
	* Statistik deskriptif tertimbang populasi antar-kabupaten/kota
	summarize avg_wage avg_cons rls_kab share_formal infra_kab rasio_sek umk_riil [aw = est_populasi]


** =============================================================================
** TAHAP 3: SISTEM PERSAMAAN SIMULTAN DENGAN INFERENSI TERTIMBANG
** =============================================================================
	use "data_agregat_simulasi_kabkota.dta", clear

	/* -------------------------------------------------------------------------
	   SPESIFIKASI STRUKTURAL & INFERENSI TERTIMBANG POPULASI:
	   - Persamaan 1 (Upah Regional):
	     avg_wage = a0 + a1*avg_cons + a2*rls_kab + a3*infra_kab + a4*umk_riil + e1
	   - Persamaan 2 (Konsumsi Regional):
	     avg_cons = b0 + b1*avg_wage + b2*share_formal + b3*rasio_sek + b4*avg_age + e2
	   
	   CATATAN DIAGNOSTIK ESTIMASI BERBOBOT:
	   - Penggunaan [pw = est_populasi] dan vce(robust) memberikan estimasi dan
	     robust standard error yang mempertimbangkan ukuran populasi tiap daerah.
	   - Pada postestimation ivregress 2sls berbobot:
	     * estat endogenous & estat overid memerlukan opsi ", forceweights"
	     * estat firststage menggunakan opsi ", all" (tanpa forceweights)
	   ------------------------------------------------------------------------- */

	// A. BASELINE OLS TERTIMBANG (MASIH MENGANDUNG SIMULTANEITY BIAS)
	reg avg_wage avg_cons rls_kab infra_kab umk_riil [pw = est_populasi], vce(robust)
	estimates store ols_w_eq1
	
	reg avg_cons avg_wage share_formal rasio_sek avg_age [pw = est_populasi], vce(robust)
	estimates store ols_w_eq2

	// B. TWO-STAGE LEAST SQUARES (2SLS) TERTIMBANG + UJI DIAGNOSTIK INSTRUMEN
	* Persamaan 1: Menginstrumenkan avg_cons dengan variabel ekslusi dari Persamaan 2
	ivregress 2sls avg_wage rls_kab infra_kab umk_riil (avg_cons = share_formal rasio_sek avg_age) [pw = est_populasi], vce(robust) first
	estimates store tsls_w_eq1
	
	* Uji Diagnostik Ekonometrika Persamaan 1:
	estat endogenous, forceweights
	estat firststage, all
	estat overid, forceweights
	
	* Persamaan 2: Menginstrumenkan avg_wage dengan variabel ekslusi dari Persamaan 1
	ivregress 2sls avg_cons share_formal rasio_sek avg_age (avg_wage = rls_kab infra_kab umk_riil) [pw = est_populasi], vce(robust) first
	estimates store tsls_w_eq2
	
	* Uji Diagnostik Ekonometrika Persamaan 2:
	estat endogenous, forceweights
	estat firststage, all
	estat overid, forceweights

	// C. THREE-STAGE LEAST SQUARES (3SLS) TERTIMBANG (FULL INFORMATION SYSTEM)
	* Mengestimasi kedua persamaan secara serentak dan memperhitungkan korelasi
	* antar-residual (Cov(e1, e2) != 0) dengan bobot populasi [aw = est_populasi]
	reg3 (eq_wage: avg_wage avg_cons rls_kab infra_kab umk_riil) (eq_cons: avg_cons avg_wage share_formal rasio_sek avg_age) [aw = est_populasi], endog(avg_wage avg_cons) exog(rls_kab infra_kab umk_riil share_formal rasio_sek avg_age) 3sls
	estimates store threesls_w_sys

	// D. GENERALIZED METHOD OF MOMENTS (IV-GMM) TERTIMBANG
	* Lebih efisien daripada 2SLS ketika terdapat heteroskedastisitas pada persamaan
	* yang overidentified dengan matriks pembobot optimal wmatrix(robust)
	ivregress gmm avg_wage rls_kab infra_kab umk_riil (avg_cons = share_formal rasio_sek avg_age) [pw = est_populasi], wmatrix(robust)
	estimates store gmm_w_eq1

	// E. TABEL PERBANDINGAN INFERENSI TERTIMBANG: OLS VS 2SLS VS 3SLS VS IV-GMM
	estimates table ols_w_eq1 tsls_w_eq1 threesls_w_sys gmm_w_eq1, b(%9.4f) se(%9.4f) stats(N r2) title("Simulasi Sistem Persamaan Simultan Tertimbang (Weighted Inference)")


** =============================================================================
** TAHAP 4: ANALISIS EKONOMETRIKA LANJUTAN (MAKRO & DESAIN SURVEI KOMPLEKS MIKRO)
** =============================================================================
	// 1. SEEMINGLY UNRELATED REGRESSION (SUR) TERTIMBANG PADA REDUCED-FORM
	* Menguji korelasi galat antar-persamaan bentuk tereduksi (Breusch-Pagan test)
	sureg (rf_wage: avg_wage rls_kab infra_kab umk_riil share_formal rasio_sek avg_age) (rf_cons: avg_cons rls_kab infra_kab umk_riil share_formal rasio_sek avg_age) [aw = est_populasi], corr

	// 2. QUANTILE REGRESSION TERTIMBANG (Q25, Q50, Q75)
	* Melihat perbedaan dampak variabel penjelas pada kabupaten/kota berupah rendah (Q25),
	* menengah (Q50/median), dan tinggi (Q75) dengan bobot populasi [pw = est_populasi]
	qreg avg_wage rls_kab infra_kab umk_riil share_formal [pw = est_populasi], quantile(0.25) vce(robust)
	estimates store q25_w
	qreg avg_wage rls_kab infra_kab umk_riil share_formal [pw = est_populasi], quantile(0.50) vce(robust)
	estimates store q50_w
	qreg avg_wage rls_kab infra_kab umk_riil share_formal [pw = est_populasi], quantile(0.75) vce(robust)
	estimates store q75_w
	
	estimates table q25_w q50_w q75_w, b(%9.4f) se(%9.4f) title("Quantile Regression Tertimbang Populasi (Q25 vs Q50 vs Q75)")

	// 3. INFERENSI SURVEI KOMPLEKS LANGSUNG DI LEVEL MIKRO (N = 50.000 INDIVIDU)
	use "data_mikro_simulasi_50k.dta", clear
	
	* CATATAN DESAIN SURVEI KOMPLEKS (SVYSET):
	* Mendeklarasikan Primary Sampling Unit (PSU) pada level kabupaten/kota (id_kab)
	* dan probability weight pada final_weight. Seluruh estimasi svy: otomatis
	* menghitung standard error yang kebal terhadap klaster wilayah & bobot survei.
	svyset id_kab [pweight = final_weight]
	
	* A. Survey-Weighted 2SLS pada Data Mikro:
	* Mengatasi endogenitas pendidikan (educ_yr) akibat unobserved ability menggunakan
	* instrumen ketersediaan sekolah daerah (rasio_sek) dengan inferensi survei penuh
	svy: ivregress 2sls ln_wage age exp exp_sq female infra_kab (educ_yr = rasio_sek)
	estimates store svy_iv_mikro
	
	* B. Survey-Weighted Probit & Logit + Average Marginal Effects (AME):
	* Memodelkan probabilitas individu bekerja di sektor formal (formal = 1)
	svy: probit formal educ_yr age exp female infra_kab
	margins, dydx(*) post
	estimates store svy_probit_ame
	
	svy: logit formal educ_yr age exp female infra_kab
	margins, dydx(*) post
	estimates store svy_logit_ame
	
	estimates table svy_probit_ame svy_logit_ame, b(%9.4f) se(%9.4f) title("Survey-Weighted Average Marginal Effects: Pekerja Formal (N=50.000)")