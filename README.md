# External validation of the KTFS 

R scripts from the paper "A Canado-European external validation of the Kidney Transplant Failure Score" by Chatton et al. (2024, submitted).

* ScriptKTFSExtVal.r: Merge the cleaned databases, fit the cause-specific Cox model, and run bootstrap to compute predictive performance
* ScriptPerformanceAndFigure.r: CI 95% from the bootstrap of the previous script, survival curves, flexible calibration curves, and decision curves.
