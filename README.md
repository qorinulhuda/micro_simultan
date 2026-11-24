# Simulasi Data Mikro-Makro & Sistem Persamaan Simultan (SEM) dengan Stata

![Stata](https://img.shields.io/badge/Software-Stata_17%2F18-1A476F?style=for-the-badge&logo=stata&logoColor=white)
![Data](https://img.shields.io/badge/Dataset-Synthetic_Monte_Carlo_(50k)-2E7D32?style=for-the-badge)
![Econometrics](https://img.shields.io/badge/Metode-2SLS_%7C_3SLS_%7C_IV--GMM_%7C_SUR-B71C1C?style=for-the-badge)

> **Penyusun:** **Qorinul Huda**  
> **Fokus Kajian:** Rekayasa *Data Generating Process* (DGP) Data Mikro Sintetis, Agregasi Tertimbang (*Weighted Collapse*), Inferensi Desain Survei Kompleks (`svy`), dan Sistem Persamaan Simultan Regional.

---

## Pendahuluan

Repositori ini merupakan eksperimen simulasi Monte Carlo berbasis **data sintetis (*synthetic data*)** yang dirancang untuk mempraktikkan alur kerja ekonometrika mikro-makro secara menyeluruh di Stata. Simulasi ini mereplikasi proses pembentukan **50.000 observasi individu** di **100 kabupaten/kota**, mengagregasinya ke tingkat daerah menggunakan penimbang survei (*survey weights*), serta menguji konsistensi berbagai estimator ketika terdapat masalah **endogenitas** dan **kausalitas dua arah (*simultaneity bias*)**.

Materi dan metode ekonometrika yang dibahas dalam repositori ini meliputi:

* **Rekayasa *Data Generating Process* (DGP) Data Mikro (`N = 50.000`):**
  Pembangkitan dataset individu sintetis dengan variabel laten tak teramati (*unobserved ability* dan heterogenitas wilayah) yang memicu bias endogenitas pada persamaan upah (*Mincer Earnings Function*), lama sekolah, status pekerjaan formal, dan pengeluaran konsumsi, lengkap dengan bobot survei heterogen (`final_weight`).
* **Agregasi Tertimbang Tingkat Kabupaten/Kota (*Weighted Collapse*):**
  Transformasi data mikro individu menjadi **100 unit agregat kabupaten/kota** menggunakan `collapse` berbobot `[aw = final_weight]`, penyimpanan estimasi total populasi (`est_populasi`), serta penyuntikan umpan balik simultan (*two-way causality*) antara rata-rata upah (`avg_wage`) dan konsumsi regional (`avg_cons`).
* **Jenis Pemodelan Ekonometrika & Uji Diagnostik yang Diterapkan:**
  1. **Baseline OLS Tertimbang (`reg ... [pw = est_populasi], vce(robust)`):** Estimasi pembanding awal untuk mendemonstrasikan bias simultanitas (*simultaneity bias*) pada persamaan struktural.
  2. **Two-Stage Least Squares (`ivregress 2sls`) & Uji Diagnostik Instrumen:** Estimasi variabel instrumen dengan pembobot populasi `[pw = est_populasi]`, disertai **Uji Endogenitas Durbin-Wu-Hausman** (`estat endogenous, forceweights`), **Uji Kekuatan Instrumen** (`estat firststage, all`), dan **Uji Validitas Overidentifikasi Hansen J** (`estat overid, forceweights`).
  3. **Three-Stage Least Squares (`reg3 ... [aw = est_populasi], 3sls`):** Estimasi sistem penuh (*Full Information System*) yang mengestimasi persamaan upah dan persamaan konsumsi secara serentak dengan memperhitungkan korelasi galat antar-persamaan ($\text{Cov}(e_1, e_2) \neq 0$).
  4. **Generalized Method of Moments (`ivregress gmm ... [pw = est_populasi], wmatrix(robust)`):** Estimasi IV-GMM dengan matriks pembobot optimal untuk efisiensi yang lebih tinggi pada model *overidentified* yang mengalami heteroskedastisitas.
  5. **Seemingly Unrelated Regression (`sureg ... [aw = est_populasi], corr`):** Estimasi sistem persamaan bentuk tereduksi (*reduced-form*) beserta uji korelasi residual lintas persamaan (*Breusch-Pagan test of independence*).
  6. **Quantile Regression Tertimbang (`qreg ... [pw = est_populasi], vce(robust)`):** Analisis heterogenitas determinan upah daerah pada **kuantil ke-25 (Q25)**, **median (Q50)**, dan **kuantil ke-75 (Q75)**.
  7. **Inferensi Survei Kompleks Tingkat Mikro (`svyset` & `svy:`):** Deklarasi klaster kabupaten/kota sebagai *Primary Sampling Unit* (`id_kab`) dan bobot survei (`final_weight`) untuk mengestimasi **Survey-Weighted 2SLS** (`svy: ivregress 2sls`) serta **Survey-Weighted Probit & Logit** beserta ***Average Marginal Effects* (`margins, dydx(*)`)** pada $50.000$ observasi individu.

---

## Struktur File Repositori

```text
├── Persamaan Simultan.do               # Skrip utama Stata (DGP mikro, collapse, SEM, & uji diagnostik)
├── data_mikro_simulasi_50k.dta         # Dataset sintetis tingkat individu (N = 50.000, 100 kab/kota)
├── data_agregat_simulasi_kabkota.dta   # Dataset hasil agregasi tertimbang tingkat kab/kota (N = 100)
└── README.md                           # Dokumentasi utama repositori
