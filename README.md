# Post-PIV-Analysis-Tools-for-Moving-Airfoil

Analyzing PIV data of moving objects (e.g., airfoils, wings) requires **Coordinate transformation** to shift the origin to the object's reference frame to compute velocity fields (average/RMS) and fluxes at a specific distance with respect to the airfoil.

This toolkit automates these steps for **experimental setups with moving airfoils**, where the object enters/exits the camera frame.

Based on the location of the trailing edge, the coordinate system of the images is transposed. This approach only shifts the origin. The values of the velocity fields continue to be with respect to an external lab reference frame. The velocity values themselves are not transformed to the airfoil reference frame.

The end result is the average velocity and fluxes with respect to the airfoil coordinate system.

## Usage

1. Clone the repository:
  
  ```bash
  git clone https://github.com/rakshithajoshi/Post-PIV-Analysis-Tools-for-Moving-Airfoil
  cd Post-PIV-Analysis-Tools-for-Moving-Airfoil
  ```
  
2. Load
   - **MATLAB** (tested on R2020a+)
   - **PIVlab** (for initial PIV processing)
   - **Image Processing Toolbox** (MATLAB)
     
3. Prerequisites:
   - PIVlab velocity fields (`.mat` files) in workspace
   - Mean path of the airfoil stored in `.dat` file
   - The location of the trailing edge in `.dat` file
     
5. Instructions:
   - Run Post-PIV-Analysis-Tools-for-Moving-Airfoil.m
   - Input experimental parameters (pitching point, frequency, amplitude, trial number)
   - It will fetch data from mean path of the airfoil and the location of the trailing edge. 
    
6. Output:
   - Average velocity and fluxes with respect to the airfoil coordinate system.

## References

> [PIVlab Documentation](https://pivlab.blogspot.com/) or https://www.pivlab.de/

> Thielicke, W. and Stamhuis, E.J. (2014): PIVlab – Towards User-friendly, Affordable and Accurate Digital Particle Image Velocimetry in MATLAB. Journal of Open Research Software 2(1):e30, DOI: [https://doi.org/10.5334/jors.bl](https://doi.org/10.5334/jors.bl)

## How to Cite

If you use this toolkit in your research, please cite it as follows:

> Joshi, R. U. (2026). Post PIV Analysis Tools for Moving Airfoil [Computer software]. https://github.com/rakshithajoshi/Post-PIV-Analysis-Tools-for-Moving-Airfoil

###
