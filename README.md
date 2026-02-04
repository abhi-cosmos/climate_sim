# Climate Lab

C++ application for interactive visualization of NASA GISTEMP global temperature anomaly data.

![Climate Lab Screenshot](screenshot.png)

## Features

- **Real-time NetCDF4 processing** of 140+ years of gridded temperature data
- **Interactive world map** with click-to-explore cell analysis
- **Time-series visualization** with ImPlot (global and per-cell)
- **OpenGL texture rendering** for smooth map interaction
- **Custom color mapping** with adjustable anomaly scales

## Tech Stack

- C++17
- NetCDF4 (climate data I/O)
- ImGui + ImPlot (UI and plotting)
- GLFW + OpenGL2 (windowing and rendering)

## Building

```bash
# Install dependencies (Ubuntu/Debian)


# Build

# Run (requires data/gistemp1200_GHCNv4_ERSSTv5.nc)
cmake --build build                                                          
./build/climate_lab 


## Data
Uses NASA GISTEMP v4 surface temperature analysis. Download from:
https://data.giss.nasa.gov/gistemp/

Place in data/ directory.


##License 
MIT License 