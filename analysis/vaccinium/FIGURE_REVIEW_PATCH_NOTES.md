# Figure review patch — server rerender

The archived `make_main_figures_v3.R` now writes Figure 1 as a three-panel
study-aware overview:

- panel A: seven paired cohorts → four discovery cohorts → study-held-out benchmark → three independent external cohorts;
- panel B: paired-sample QC retention;
- panel C: high-confidence bacterial and fungal read depth.

On the server, replace the existing source with this archived version and rerender:

```bash
cd ~/projects/blueberry_meta

cp -p make_main_figures_v3.R make_main_figures_v3.before_figure1_workflow.R
# Upload/copy this package's:
# 06_Data_Code_Archive/scripts/make_main_figures_v3.R
# to ~/projects/blueberry_meta/make_main_figures_v3.R before the next command.

export R_LIBS_USER="$HOME/R/library"
Rscript --vanilla make_main_figures_v3.R \
  2>&1 | tee figure1_workflow_rerender.log

find global_harmonized/figures_nature -maxdepth 1 -type f \
  -name 'Figure1_study_design_QC_nature.*' -printf '%TY-%Tm-%Td %TH:%TM:%TS  %s  %p\n' | sort
sha256sum \
  global_harmonized/figures_nature/Figure1_study_design_QC_nature.pdf \
  global_harmonized/figures_nature/Figure1_study_design_QC_nature.png
```

Copy the regenerated PDF/PNG into `02_Main_Figures/` after checking that the
workflow panel is legible at 183 mm width. Figure 4 and Figure 6 titles/captions
have already been synchronized in the manuscript, legends, and
`figure_source_data/Figure_Captions.txt`; Figure 4 is explicitly labelled
exploratory and Figure 6 is explicitly labelled study-held-out external validation.
