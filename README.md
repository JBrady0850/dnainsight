DNAInsight v1.0

![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![SNPs](https://img.shields.io/badge/Bundled_SNPs-101-orange)

Personal DNA Analysis Tool — Process your raw DNA file from any major provider, annotate SNPs against ClinVar and PharmGKB, generate a Genetic Health Report and a Doctor Discussion Report, and analyze results with Grok or any AI assistant.

---

What This Tool Does

DNAInsight reads your raw DNA file (the uncompressed text file you downloaded from your DNA testing company) and:

1. Parses your 600,000+ SNPs into a local database.
2. Annotates each SNP against a bundled clinical reference (101 high-priority medical SNPs, offline, instant) and optionally the MyVariant.info API (ClinVar + PharmGKB, requires internet, no personal data transmitted).
3. Classifies findings into three tiers: Prescription-Critical, Actionable, and Informational.
4. Generates two HTML reports: a Genetic Health Report and a Doctor Discussion Report.
5. Provides a Grok-compatible AI prompt for AI-assisted clinical interpretation.

Privacy: All processing is local. Your DNA data never leaves your computer unless you explicitly enable API annotation (which sends only rsIDs, never genotypes, to myvariant.info).

---

Supported DNA Providers

| Provider | Format | Notes |
|----------|--------|-------|
| AncestryDNA | .txt (tab-delimited) | V1 and V2 arrays supported |
| 23andMe | .txt (tab-delimited) | All array versions |
| MyHeritage | .csv or .txt | |
| FamilyTreeDNA (FTDNA) | .csv | |
| LivingDNA | .txt | |
| Generic TSV | .txt / .csv | Auto-detected column layout |

Note: Upload only the uncompressed raw file. Do NOT upload the .zip file from your provider — extract it first.

---

Requirements

- Python 3.10 or newer (auto-installed by the installer if missing)
- Internet connection (optional, for extended API annotation and database updates)
- ~50 MB disk space

No other software required. All dependencies install automatically.

---

Installation

Windows

1. Extract the DNAInsight folder anywhere on your computer.
2. Double-click install.bat

The installer automatically detects whether Python is installed. If Python is not found:
- It tries Windows Package Manager (winget) first.
- If winget is unavailable, it downloads and installs Python 3.12 silently.
- If Python is auto-installed, the installer will prompt you to close and re-run it once (PATH refresh required).

All Python packages (Flask, requests) install automatically. No manual steps required.

macOS / Linux

1. Open Terminal in the DNAInsight folder.
2. Run: bash install.sh

The installer automatically detects whether Python 3 is installed and uses the appropriate
package manager (Homebrew, apt, dnf, yum, or pacman) to install it. If the system pip
conflicts with PEP 668 restrictions (Ubuntu 23.04+, Debian 12+), the installer creates
a virtual environment automatically.

---

Running DNAInsight

Windows: Double-click launch.bat

macOS/Linux: ./launch.sh or python3 app.py

Any platform: python app.py

DNAInsight opens automatically in your browser at http://127.0.0.1:5050

Press Ctrl+C in the terminal to stop the server.

---

User Walkthrough

Step 1: Download Your Raw DNA File

AncestryDNA:
1. Log in at ancestry.com
2. Click your name > DNA > Settings
3. Scroll to "Download DNA Data" > click Download
4. Extract the zip; use the file named AncestryDNA.txt

23andMe:
1. Log in at 23andme.com
2. Browse to Tools > Browse Raw Data > Download
3. Request download; wait for the email link
4. Extract and use the file named genome_*.txt

MyHeritage:
1. Log in > DNA > Manage DNA Kits > Download Raw Data

FamilyTreeDNA:
1. Log in > myDNA > Chromosome Browser > Download Raw Data

Step 2: Add a Profile

1. Open DNAInsight in your browser.
2. Click "+ Add Profile" in the left sidebar.
3. Enter your name, date of birth, and biological sex.
4. Upload your raw DNA file by clicking or dragging it onto the upload area.
5. Click Next then Upload and Create Profile.

DNAInsight will parse the file and report how many SNPs were loaded (typically 600,000-700,000).

Step 3: Run a Scan

1. Click Run Scan in the navigation.
2. Choose whether to use the online API (recommended for more complete results).
3. Click Start Scan.

The scan has two phases:
- Bundled reference (instant, offline): Annotates 91 high-priority medical SNPs immediately.
- API annotation (if enabled): Queries myvariant.info in batches. This may take several minutes for a full file, but DNAInsight will show progress.

Step 4: Review Findings

Click Findings in the navigation to view your annotated SNPs.

| Tier | Meaning |
|------|---------|
| Prescription-Critical | Variants that affect drug metabolism, dosing, or contraindications. Show these to your prescriber before any medication changes. |
| Actionable | Variants with lifestyle, supplement, or monitoring implications. |
| Informational | Background genetic information; no immediate action required. |

Use the search bar and gene filter to explore findings by gene, rsID, or keyword.

Step 5: Generate Reports

Click Reports in the navigation.

Genetic Health Report: A comprehensive annotated report covering all findings organized by tier and gene pathway. Designed for personal use and general health review.

Doctor Discussion Report: A clinical-style document designed to bring to a physician, pharmacist, or genetic counselor. Includes:
- Prescription-critical variant table
- Drug class interaction summary
- Recommended lab follow-ups
- AI analysis prompt block (ready to paste into Grok or Claude)

Both reports open in your browser and can be printed or saved as PDF using your browser's print function (Ctrl+P > Save as PDF).

Step 6: AI-Assisted Analysis with Grok

1. Open your Doctor Discussion Report.
2. Copy the text in the AI-Assisted Analysis Prompt section.
3. Open Grok at https://x.ai/grok or the Grok desktop app.
4. Paste the prompt into Grok.

Alternatively, use the standalone prompt in grok/GROK_SYSTEM_PROMPT.md with any AI assistant.

---

Multiple Family Members

DNAInsight supports multiple profiles. Each profile is stored separately in the local database. Add as many family members as needed by clicking "+ Add Profile" and uploading their raw DNA file.

---

Updating the SNP Reference Database

The bundled SNP reference can be refreshed directly from DNAInsight without any command-line steps.

In-app update (recommended, monthly):

1. Open DNAInsight.
2. Click Database in the sidebar or top navigation.
3. Click Update Databases Now.

DNAInsight re-fetches the latest ClinVar clinical significance ratings and associated
condition names for all 101 bundled SNPs from MyVariant.info. The update takes 1-3 minutes
and requires internet access. Only rsIDs are transmitted -- no genotype data leaves your
computer.

DNAInsight displays a warning banner if the reference has not been updated in 30 or more days.

Manual rebuild (advanced):

To regenerate the reference from the curated source list in build_reference.py:
  python data/build_reference.py

This resets clinical significance and interpretation back to the curated static values.

---

Contributing & Extending the SNP Reference

Evidence-based SNP additions are welcome. To add entries:

1. Edit data/build_reference.py and add a new tuple to the REFERENCE list using the format: (rsID, gene, category, clinical_sig, interpretation)
2. Run: python data/build_reference.py
3. Verify by restarting DNAInsight and scanning a test DNA file.
4. Submit a Pull Request with supporting references (CPIC, PharmGKB, ClinVar, or GWAS catalog).

Criteria for inclusion:
- Covered on major consumer arrays (23andMe v4/v5, AncestryDNA v2)
- CPIC Level A or B evidence, high ClinVar significance, or replicated GWAS hits
- Actionable (lifestyle, supplement, or physician discussion implication)
- Plain-English interpretation targeted at non-experts

---

Bundled SNP Reference (101 High-Priority SNPs)

The file data/snp_reference.json contains 101 carefully curated SNPs focused on maximum actionability for consumer DNA arrays. These are prioritized for clear lifestyle, supplement, or physician discussion implications.

| Category | Focus | Key Genes |
|----------|-------|-----------|
| Pharmacogenomics (PHARM) | CPIC-level drug response: warfarin, statins, antidepressants, opioids, antiplatelet agents, thiopurines | CYP2D6, CYP2C19, CYP2C9, CYP4F2, VKORC1, SLCO1B1, ABCB1, TPMT, APOE |
| Metabolic Health (METAB) | Obesity, T2D, nutrient processing, iron overload | FTO, TCF7L2, PPARG, MTHFR, HFE, SLC30A8, KCNJ11 |
| Inflammation (INFLAM) | Chronic inflammation, autoimmune susceptibility | IL6, TNFA, IL10, CRP, CTLA4, IL6R |
| Neurological (NEURO) | Mood, folate cycle, stress response, neurotransmitters, social behavior | MTHFR, COMT, BDNF, MAOA, SLC6A4, FKBP5, CLOCK, OXTR |
| Detox & Cardio (DETOX/CARDIO) | Oxidative stress, clotting risk, alcohol metabolism, nicotine dependence | GSTP1, SOD2, NQO1, ALDH2, F5, F2, CYP1A2, CHRNA3 |

**Why only 101?** Consumer arrays have limited SNP coverage compared to clinical sequencing. This reference focuses on well-covered, high-evidence SNPs with CPIC Level A/B support or strong replicated GWAS associations.

To view or extend the full curated list, open data/build_reference.py.
To rebuild the reference after edits: python data/build_reference.py

---

Privacy and Data Security

- All DNA data is stored locally in dnainsight.db (SQLite) on your computer.
- No data is uploaded to any server.
- When API annotation is enabled, only rsIDs (not genotypes) are sent to myvariant.info. rsIDs are not personally identifiable.
- To delete all data for a profile, use the Delete button in the Profiles view.

---

Disclaimer

DNAInsight is not a medical device and does not provide medical advice. Consumer DNA arrays are not clinical-grade tests and have significantly limited coverage compared to clinical exome or genome sequencing. Negative results do not rule out genetic risks. Results must not be used to make prescribing or diagnostic decisions without consultation with a licensed healthcare provider or clinical geneticist. All findings require clinical validation before any medical action.

For questions about your genetic results, consider consulting a genetic counselor. The American Board of Genetic Counseling maintains a directory at findageneticcounselor.com.

The information provided is for personal educational use only.

---

File Structure


DNAInsight/
├── app.py                    Entry point (Flask server)
├── requirements.txt          Python dependencies
├── install.bat               Windows installer
├── install.sh                macOS/Linux installer
├── launch.bat                Windows launcher (created by installer)
├── README.md                 This file
├── dnainsight.db             SQLite database (created on first run)
├── backend/
│   ├── __init__.py
│   ├── parsers.py            DNA file format parsers
│   ├── scanner.py            SNP annotation engine
│   ├── database.py           SQLite data access layer
│   ├── genetic_report.py     Genetic Health Report generator
│   ├── doctor_report.py      Doctor Discussion Report generator
│   └── routes.py             Flask API endpoints
├── data/
│   ├── build_reference.py    Script to regenerate bundled SNP reference
│   └── snp_reference.json    Bundled SNP reference (91 SNPs)
├── frontend/
│   └── index.html            Single-page web application
├── grok/
│   └── GROK_SYSTEM_PROMPT.md  AI analysis prompt for Grok/Claude/ChatGPT
├── uploads/                  Uploaded DNA files (created on first run)
└── reports_output/           Generated HTML reports (created on first run)


---

Troubleshooting

"No module named flask"
Run: pip install flask or python3 -m pip install flask

"Parsed 0 valid SNPs"
The file may be compressed. Open the zip from your DNA provider and use the .txt file inside.

Browser does not open automatically
Navigate to http://127.0.0.1:5050 manually.

Port 5050 already in use
Run: python app.py --port 8080 and open http://127.0.0.1:8080

Scan takes a long time
The API phase processes SNPs in batches of 200 with 4 parallel threads. For a full 650,000-SNP file, a full API scan can take 30-60 minutes. The bundled reference scan (no API) is instant.

---

Screenshots

(Screenshots coming soon -- UI, findings view, and sample report.)

---

Roadmap

| Version | Planned Features |
|---------|-----------------|
| v1.1 | Expanded reference (120+ SNPs), PDF export, multi-profile comparison charts |
| v1.2 | Basic polygenic risk score (PRS) support for T2D and CAD using public weight files |
| v1.3 | CYP2D6 star allele calling, Docker support, improved risk visualizations (Chart.js heatmaps) |
| Future | Community-curated reference updates, auto-update mechanism for new app versions |

Contributions and feedback are welcome. Open an issue or submit a pull request.

---

License

MIT License. Free to use, modify, and distribute for personal and non-commercial use.
