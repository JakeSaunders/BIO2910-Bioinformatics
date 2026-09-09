# BIO 2910 — Bioinformatics

### Kean University · Department of Biology · Saunders Lab

Course materials for BIO 2910. Every coding activity in this course runs in
**[Posit Cloud](https://posit.cloud)** — RStudio in your web browser.
**You do not need to install anything, on any computer.**

Advanced users can chose to install both R and RStudio on their comupters if they wish. Most of the code should work in the native installation with a few minimal change to the pre-written working directory. 

---

## :computer: Set up your account (before the first R class)

1. Go to **[posit.cloud](https://posit.cloud)** and click **Sign Up**. Choose the **free** plan
   and **use your `kean.edu` email address.**
2. Confirm the email Posit sends you.
3. Once you login to posit.cloud. Click on **New Project** and select **New RStudio Project**
4. Once Rstudio Cloud opens, I recommend changing the name of the project from **Untitled Project** to something like **BIO2910**.
5. Before each Unit you will be given code to run that will download the required files for each unit. 

If you have issues getting this to work please message me **before** class rather than waiting for class to trouble shoot the problem.

---

## :arrow_down: Get the files for a unit

In the **Console** (bottom left panel), run these three lines:

```r
source("https://raw.githubusercontent.com/JakeSaunders/BIO2910-Bioinformatics/main/getUnit.R")
getUnit("Unit02")
```

Change `"Unit02"` to whichever unit we are on. The unit's files then appear in the **Files** panel
(bottom right) — click a `.R` file to open it in the editor.

**To start over.** If you want to throw away your changes and go back to the original files:

```r
getUnit("Unit02", overwrite = TRUE)
```

This deletes your copy of that unit, so only run it if you mean it.

---

## Units

| Unit | Topic | Get the files |
|---|---|---|
| **[Unit02](Unit02)** | Basics of R and RStudio | `getUnit("Unit02")` |

---

## Three things that will save you a headache

**Type your code in the `.R` file, not the Console.** The Console forgets. The file is what gets
saved and what you turn in. Run the current line with **Ctrl+Enter** (**Cmd+Enter** on a Mac).

**Your work is saved and it persists.** Close the tab, come back next week, and your project is
still there — your files, your variables, your packages. Nothing to re-download.

**Close the tab when you are done.** After 15 minutes of sitting idle your session goes to sleep.
Nothing is lost, but an open project keeps spending the course's computing budget. Don't install
packages unless an assignment tells you to — everything you need is already there.

---

## :outbox_tray: Turning in your work

1. In the **Files** panel, check the box next to the file you want to submit.
2. **More → Export.** The file downloads to your computer.
3. Upload it to the assignment in **Canvas**.

Canvas is the system of record. A file sitting in Posit Cloud has not been submitted.

---

## Data

Shared data files live in [`data/`](data). You do not need to download them separately —
`getUnit()` copies them into the unit folder for you, so a script can read
`data/Bio2910-Finches.csv` once you have run `setwd()`.

## License

Course materials © Cecil Jake Saunders. Released under the MIT License.
