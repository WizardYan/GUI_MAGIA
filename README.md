
# MAGIA GUI Tutorial




## 1. Introduction

PET kinetic modeling. It allows researchers to easily configure subject data, tracer parameters, and modeling options without manual scripting. Through interactive menus, users can set environment variables, validate inputs, and launch modeling tasks on local or cluster environments. The GUI also supports batch processing of multiple subjects, dynamic parameter updates, and visualization of modeling outputs. By lowering technical barriers, the MAGIA GUI enhances reproducibility, efficiency, and accessibility, enabling neuroscientists and clinicians to focus on interpreting PET imaging results rather than managing complex code.

The **MAGIA GUI** provides a user-friendly interface for running PET kinetic modeling with the MAGIA [https://github.com/tkkarjal/magia] pipeline.  
The new version introduces:
- Multi-subject support
- Environment variable customization
- Automatic parameter loading from options/spec files
- Improved error handling and logging

**Prerequisites**
- MATLAB (R2022 or later)
- SPM, FreeSurfer, FSL installed and configured
- Preprocessed PET and MRI data]

---

## 2. Launching the GUI
1. Open MATLAB.
2. Navigate to the MAGIA installation folder.
3. Run:
   ```matlab
   run_magia_gui

## GUI Screenshot

![MAGIA GUI](main_figure.png)
