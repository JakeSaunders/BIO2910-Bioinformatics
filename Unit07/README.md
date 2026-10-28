# Unit 07 — Multiple Sequence Alignments (in Python, on Google Colab)

This unit runs in **Google Colab** instead of Posit Cloud, and in **Python** instead of R.
There is nothing to install: click a button below and the notebook opens in your browser.
You will need a Google account (a personal Gmail is fine).

| Assignment | What you'll do | Open |
|---|---|---|
| **07-01 In class** | Python crash course · catch seafood fraud with DNA barcodes · K⁺ channel sequence logo | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/JakeSaunders/BIO2910-Bioinformatics/blob/main/Unit07/BIO2910_Assignment_07-01_MSA_inclass.ipynb) |
| **07-02 Homework** | Opsins, optogenetics and the 2026 Nobel Prize · hydropathy vs. UniProt · are animal and microbial rhodopsins related? | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/JakeSaunders/BIO2910-Bioinformatics/blob/main/Unit07/BIO2910_Assignment_07-02_MSA_homework.ipynb) |

## Every time you open a notebook
1. Click the **Open in Colab** button above.
2. **File → Save a copy in Drive.** Otherwise your work is lost when you close the tab.
3. Run cells from top to bottom with **Shift + Enter**.
4. Figures you save appear in the **Files** panel (📁 on the left). Right-click a file → **Download**.

If something breaks: **Runtime → Restart session and run all**.

## Data (in [`/data`](../data))
| File | Contents | Source |
|---|---|---|
| `seafood_COI_references.fasta` | COI barcodes, 14 commercial fish species | GenBank |
| `seafood_COI_samples.fasta` | 7 "restaurant/market" samples to identify | GenBank |
| `K_channel_pore.fasta` | Pore region (pore helix + selectivity filter) of 21 K⁺ channels | UniProt |
| `opsins.fasta` | 8 animal opsins + 3 microbial rhodopsins | UniProt |
| `opsins_TM_uniprot.tsv` | UniProt transmembrane annotations (backup if UniProt is down) | UniProt |
