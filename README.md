# Codebase for the "Phylogenetic Reconstruction of Aging Rates in the Primate Lineage" manuscript.

## Code Organizational Structure

## Directory Structure
- data/ -  input datasets, including raw survival lifespan, phylogenetic tree, and common names
- plots/ - generated figures
- outputs/ - generate text files

## Instructions for Reproducibility
To reproduce the analyses and figures in the manuscript, follow the steps below:

1. Clone this repository to your local machine.
2. Request subject-level data from the PAD (see below)
    - Extract "PrimateAgingRates.zip" and place tsv files in the `data/` directory.
3. Start R (version 4.4 or later)
4. Install the required packages (see `requirements.txt` for a list of R packages and their versions).
5. execute `source("PIPELINE.R")` to run the full analysis pipeline.

## How to request subject-level data from PAD
To request access to the data, please contact administrators at PAD (Primate Aging Database) via their website: https://primatedatabase.org/authentication/registration
and request access to the "PrimateAgingRates.zip" dataset. 

### data/PAD_subjects.tsv 

| Column | Description |
|---|---|
| `subject` | Unique subject ID, encoding institution + species abbreviation + number (e.g., `Al-Ch-1`) |
| `birthdate` | Date of birth (M/D/YYYY format) |
| `deathdate` | Date of death (M/D/YYYY); empty if still alive or lost |
| `status` | Vital status at last observation: `alive`, `died`, or `lost` |
| `experimental` | Whether the subject was part of an experimental treatment: `YES` or `NO` |
| `species` | Common species name (44 species, e.g., Rhesus Indian-derived, Vervet, Chimpanzee, Common Marmoset, various lemurs) |
| `sex` | `Female` or `Male` |
| `housing` | Housing environment: `Inside`, `Outside`, `Mixed`, or `Unknown` |
| `diet` | Diet type: `Chow` (standard; ~95%) or `Other` |
| `last_measurement_date` | Date of the last recorded measurement for the subject |
| `calc_lifespan` | Calculated lifespan in decimal years (0–63.27); age at death or last observation |
| `last_transfer_date` | Date of last institutional transfer (mostly empty) |
| `last_observation_date` | Date of last known observation (used for censoring in survival analysis) |
| `transfer_comments` | Free-text notes on transfers or losses (mostly empty) |


### data/PAD_subjects_measurements_quarterly_50th_BodyWeight.tsv

| Column | Description |
|---|---|
| `subject` | Unique subject ID |
| y0 | body weight measurement at birth (0 years) |
| y0.25 | body weight measurement at 0.25 years (3 months) |
etc. 




## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
