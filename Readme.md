<div style="text-align: center; line-height: 0;">
  <!-- link: https://readme-typing-svg.demolab.com/demo/ and https://devicon.dev/ -->
  <a href="https://git.io/typing-svg"><img src="https://readme-typing-svg.demolab.com?font=Fira+Code&size30&pause=1000&center=true&multiline=true&width=1500&height=100&lines=Targeted+memory+reactivation+elicits+temporally+compressed+reactivation+linked+to+spindles" alt="Typing SVG" /></a>
<div align="center">  
  <img src="https://i.giphy.com/media/v1.Y2lkPTc5MGI3NjExdDg3b3psajducG85OWx3ZWoyNDgwNTdzMHV3Zjh2OGhqcGc1am5haCZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/M3qFDl8PLYHZMFmmdZ/giphy.gif" width="600" alt="Example GIF" style="margin-top: -50px; margin-bottom: 50px;">
</div>
</div>

<div align="center">

## 🎇 Thanks for exploring the scripts of our work. We hope you find them useful and share our enthusiasm for this research

### 📜 Main Pipeline

The core of our analysis is in `jitteredCues_mclass_classification_pipeline.m`, with figure labels in the comments to reproduce any figure with ease.

### 📂 Data

Data is publicly available on OSF here: [OSF link](https://osf.io/byvcg/?view_only=9b149e0387814bf1a6fca692f90e9167). Please put all folders (code, data, fieldtrip-20190419) in the same main folder, as this is going to be used in the code.

### 🔧 Supporting Files

Files starting with `lv` are used by the main pipeline. They're part of the Lively Vectors (LV) MATLAB toolbox by Mahmoud Abdellahi.

  <a href="https://github.com/MahmoudAbdellahi/Lively_Vectors" target="_blank">
    <img src="https://img.shields.io/badge/GitHub-Lively_Vectors-blue?style=for-the-badge&logo=github" alt="Lively Vectors GitHub">
  </a>

</div>

---

## 🚀 How to Use This Code

1. **Choose Your Analysis**: Set the `analysis` variable to one of:

   - `'erp_tf'` for erp and time frequence analyses (figure 2b)
   - `'classification'` multi-class classification of reactivations (figure 2c)
   - `'high_spindle_pw'` high spindle power classification (figure 4b)
   - `'low_spindle_pw'` low spindle power classification (supplementary figure 2)
   - `'temporal_compression'` fast memory reactivation in sleep (figure 3)
   - `'behavioural_analyses'` analyses of behavioural performance and correlation with classification (figure 2a, figure 2d, figure 4a)

2. **Find Figures**: Figure names produced from the code match those in the article.

3. **Two-Part Structure**:

   - Part 1: Loops through participants for the selected analysis.
   - Part 2: Performs statistical analyses and plots results.

   💡 **Tip**: Set a breakpoint if you like to have a look at different segments.

---

## 💻 System Requirements

### Minimum Specs:

- CPU: i7 6700K (4 GHz, 4 Core, 8 Threads, 8 MB cache)
- RAM: 16 GB
- GPU: Nvidia GTX 980M
- OS: Windows 8, Linux, OSX
- MATLAB: 2018a
- ⏱️ Run Time: &lt; 20 hours

### Recommended Specs:

- CPU: i7 11700K (3.6 GHz, 4 Core, 16 Threads, 16 MB cache)
- RAM: 32 GB
- GPU: Nvidia RTX 3060 or higher
- OS: Windows 10/11, NVIDIA-supported Linux distros, OSX
- MATLAB: 2023a
- ⏱️ Run Time: &lt; 6 hours

### ⚠️ Notes:

- Run times are for all analyses back-to-back.
- Individual analyses vary: ERP analysis (&lt; 5 minutes) vs. Temporal compression (several hours).

### 📊 Data Availability:

All generated data is also provided for convenience. You could load them directly or use the code to regenerate if you like.
