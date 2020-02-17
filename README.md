# ChAT Project

## Introduction

These are the analysis files used in the following publication:

- Taylor R Fore et al. (2020) "Acetylcholine modulates cerebellar granule cell spiking by regulating the balance of synaptic excitation and inhibition" Journal of Neuroscience

## How to Use

### Requirements

- R language: [R] (R >= 3.2.0)
- R package: [matools]
- R package: [googlesheets4] (required if running the scripts as CRON jobs)
- R package: [utils] (core package - required if running the scripts as CRON jobs)
- RStudio: [RStudio] (IDE for interactive analysis)

### Running Workflows

1. The workflow scripts in the repository are fully automated and expect the directory structure as follows:

```sh
├── analysis_scripts
│   ├── evoked_ap_analysis.R
│   └── evoked_psp_analysis.R
├── data
│   ├── cell_attached_grc_evoked_ap
│   │   ├── chat_evoked_ap_parameters.csv
│   │   └── data
│   │       ├── rec000000.ASC
│   │       ├── rec000000.txt
│   │       ...
│   │       ├── rec00000x.ASC
│   │       └── rec00000x.txt
│   └── evoked_psp
│       ├── chat_evoked_psp_parameters.csv
│       └── data
│           ├── rec000000.ASC
│           ...
│           └── rec00000x.ASC
```

2. Clone the repository and change into the `analysis_scripts` directory

```sh
$ git clone https://github.com/trfore/chat-project.git
$ cd chat-project
```

3. Run the scripts. PDF figures are placed within the parent directory, `chat-project`.

```R
# Evoked Granule Cell Action Potential Analysis
Rscript analysis_scripts/evoked_ap_analysis.R

# Evoked Post-Synaptic Potential Analysis
Rscript analysis_scripts/evoked_psp_analysis.R
```

4. Alternatively, interactively run these scripts within [RStudio].

### Running Workflows as CRON jobs

1. Modify a workflow by changing the following:

```R
# Project Data Processing ----
# Local data files
data_path <- paste(getwd(), "data/evoked_psp", sep = "/")
data_folder <- paste(data_path, "data", sep = "/")

# Option 1: Using local parameter sheet
data_parameters <- paste(data_path, "chat_evoked_psp_parameters.csv", sep = "/")

# Option 2: Using publicly shared Google Sheet
googlesheets4::gs4_deauth()
public_sheet <- "https://docs.google.com/spreadsheets/d/YOUR_SHEET_HASH_HERE/edit?usp=sharing"
data_parameters <-
  googlesheets4::read_sheet(public_sheet)

# Import the parameters
import_experiment_parameters(data_parameters)
```

2. Create a cron job

```sh
# crontab file
# run daily at midnight
0 0 * * * cd /home/chat-project/analysis_scripts; Rscript evoked_ap_analysis.R >/dev/null 2>&1
# run daily at 12:30 am
30 0 * * * cd /home/chat-project/analysis_scripts; Rscript evoked_psp_analysis.R >/dev/null 2>&1
```

# Authors

- Taylor Fore (https://github.com/trfore)

# References

- Fore et al. 2020
- https://github.com/trfore/matools
- https://github.com/trfore/chat-project
- https://github.com/trfore/chatmodel
- https://www.r-project.org/
- https://www.rstudio.com/

## R packages

- https://github.com/tidyverse/googlesheets4/
- https://github.com/HenrikBengtsson/R.utils

[googlesheets4]: https://github.com/tidyverse/googlesheets4/
[matools]: https://github.com/trfore/matools
[R]: https://www.r-project.org/
[RStudio]: https://www.rstudio.com/
[utils]: https://github.com/HenrikBengtsson/R.utils
