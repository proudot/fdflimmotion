# fdflimmotion
fdflimmotion enables live cell fluorescence lifetime estimation through the reconstruction lifetime map corrupted by motion artifacts.

Artifact origin in live FD-FLIM
![alt text](https://raw.githubusercontent.com/proudot/fdflimmotion/master/img/FD-FLIM-postproc-01.png)

Model and results
![alt text](https://raw.githubusercontent.com/proudot/fdflimmotion/master/img/model-results-01.png)

# Usage

Build 

```bash
mkdir build
cd build
cmake ..
make
```

Usage


```bash
fdflimmotion -h
```

# Manuscripts

Fluorescence lifetime is usually defined as the average nanosecond-scale delay between excitation and emission of fluorescence. It has been established that lifetime measurements yield numerous indications on cellular processes such as interprotein and intraprotein mechanisms through fluorescent tagging and Förster resonance energy transfer. In this area, frequency-domain fluorescence lifetime imaging microscopy is particularly appropriate to probe a sample noninvasively and quantify these interactions in living cells. The aim is then to measure the fluorescence lifetime in the sample at each location in space from fluorescence variations observed in a temporal sequence of images obtained by phase modulation of the detection signal. This leads to a sensitivity of lifetime determination to other sources of fluorescence variations such as intracellular motion. In this paper, we propose a robust statistical method for lifetime estimation for both background and small moving structures with a focus on intracellular vesicle trafficking.

Roudot, P., C. Kervrann, and F. Waharte. “Lifetime Estimation of Moving Vesicles in Frequency-Domain Fluorescence Lifetime Imaging Microscopy.” In IEEE International Symposium on Biomedical Imaging: Nano to Macro, 668–671, 2012.

Roudot, P., C. Kervrann, F. Waharte, and J. Boulanger. “Lifetime Map Reconstruction in Frequency-Domain Fluorescence Lifetime Imaging Microscopy.” In 2012 19th IEEE International Conference on Image Processing, 2537–40, 2012. https://doi.org/10.1109/ICIP.2012.6467415.

Roudot, P., C. Kervrann, J. Boulanger, and F. Waharte. “Noise Modeling for Intensified Camera in Fluorescence Imaging: Application to Image Denoising.” In IEEE ISBI 2013, 600–603, 2013. https://doi.org/10.1109/ISBI.2013.6556546.

Roudot, Philippe, Charles Kervrann, Cedric M. Blouin, and Francois Waharte. “Lifetime Estimation of Moving Subcellular Objects in Frequency-Domain Fluorescence Lifetime Imaging Microscopy.” Journal of the Optical Society of America A 32, no. 10 (October 1, 2015): 1821. https://doi.org/10.1364/JOSAA.32.001821.
