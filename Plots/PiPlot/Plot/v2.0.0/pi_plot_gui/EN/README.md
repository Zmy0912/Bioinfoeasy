# Pi Plot - Sliding Window Gene Annotation Tool

## Version Information

**Version**: v2.0 (GUI)
**Release Date**: 2026-03-07

## New Features

This GUI version includes the following new features based on the original command-line version:

### 1. User-Friendly Graphical Interface
- Modern interface built with tkinter
- Intuitive parameter settings panel
- Real-time run log display
- Progress indicator

### 2. Multiple Output Format Support
- **PNG Format**: Bitmap format, suitable for screen viewing
- **SVG Format**: Vector format, lossless scaling, suitable for printing and publishing
- **PDF Format**: Vector format, suitable for document embedding and printing
- **All Formats**: Generate PNG, SVG, and PDF formats simultaneously

### 3. All Original Features
- Parse Pi sliding window data files
- Parse GFF3 gene annotation files
- Three significance criteria (relative, absolute, lower threshold)
- Smart label layout
- Multi-condition gene priority sorting
- Generate detailed analysis reports

## System Requirements

- Python 3.6 or higher
- Operating System: Windows, Linux, or macOS
- tkinter (usually installed with Python)

## Installation

### 1. Install Dependencies

```bash
pip install -r requirements.txt
```

Or install manually:

```bash
pip install pandas matplotlib numpy
```

### 2. Check tkinter

tkinter is usually installed with Python. If you encounter issues, install tkinter:

- **Windows**: Usually no additional installation needed
- **Linux**: `sudo apt-get install python3-tk`
- **macOS**: Usually no additional installation needed

## Usage

### Start Program

```bash
python pi_plot_gui.py
```

### Interface Workflow

1. **Select Files**
   - Click "Browse..." button to select Pi file
   - Click "Browse..." button to select GFF3 file
   - Click "Browse..." button to select output directory

2. **Set Output Format**
   - Choose desired output format (PNG/SVG/PDF/All)
   - SVG or PDF recommended for academic publishing

3. **Adjust Plot Parameters**
   - Output filename: Default is `pi_sliding_window_plot`
   - Font size: 8-18 (recommended 10-12)
   - Max genes: 5-50 (recommended 10-15)
   - Window size: 1000-100000 bp (recommended 10000-20000)
   - Relative significance: 0.5-5.0 (recommended 1.5-3.0)
   - Absolute threshold: Optional, for absolute significance filtering
   - Lower threshold: Optional, for detecting lowest Pi value genes

4. **Start Analysis**
   - Click "Start Analysis" button
   - View run log for processing progress
   - Popup message will appear when analysis is complete

5. **View Results**
   - View generated chart files in output directory
   - Check `_info.txt` file for detailed analysis results

## Output Format Comparison

### PNG Format
- **Pros**: Small file size, easy to view
- **Cons**: Distortion after scaling
- **Use Cases**: Screen viewing, quick preview

### SVG Format (Vector Graphics)
- **Pros**: Lossless scaling, suitable for large format printing
- **Cons**: Larger file size
- **Use Cases**: Academic publishing, posters, high-quality printing

### PDF Format
- **Pros**: Vector format, suitable for document embedding
- **Cons**: Requires professional software for editing
- **Use Cases**: Academic documents, reports, printing

### All Formats
- Generate all three formats simultaneously
- Suitable for scenarios requiring multiple uses

## Parameter Tuning Recommendations

### Window Size (--window-size)
- **Small window (5000-10000 bp)**: Suitable for detecting local changes
- **Medium window (10000-20000 bp)**: Balance local and global features (recommended)
- **Large window (20000-50000 bp)**: Detect large-scale trends

### Significance Standard (--significance-std)
- **1.5**: Lenient standard, annotate more genes
- **2.0**: Standard balance (default)
- **2.5-3.0**: Strict standard, only annotate most significant genes

### Output Format Selection Recommendations

| Use Case | Recommended Format | Reason |
|----------|------------------|--------|
| Screen viewing | PNG | Fast loading, small file |
| Academic papers | SVG or PDF | Vector format, lossless scaling |
| Poster presentation | SVG | Vector format, high print quality |
| Document embedding | PDF | Good compatibility |
| Multiple uses | All formats | Meet various needs |

## Common Issues

### Q1: Program fails to start, tkinter not found

**Solution**:
- Windows: Reinstall Python, ensure "tcl/tk and IDLE" is checked
- Linux: `sudo apt-get install python3-tk`
- macOS: Install using Homebrew: `brew install python-tk`

### Q2: Chart text displays as boxes

**Solution**:
- Windows: Chinese font support built-in, should display normally
- Linux: Install Chinese font package `sudo apt-get install fonts-wqy-zenhei`

### Q3: Output files are very large after selecting SVG or PDF

**Solution**:
- Reduce number of annotated genes
- Reduce font size
- This is normal, it's a characteristic of vector formats

### Q4: Interface freezes, no response

**Solution**:
- Program runs analysis in background thread, interface should not freeze
- If issues occur, check file size and data format
- Large genomes may require longer processing time

### Q5: How to get high-quality charts for publication?

**Solution**:
- Choose SVG or PDF format
- Adjust parameters for clear charts (reduce gene count, appropriate font size)
- Use vector graphics editing software (e.g., Adobe Illustrator, Inkscape) for further editing
- PDF format can be directly embedded in LaTeX documents

## Comparison with Command Line Version

| Feature | CLI Version (v1.1.0) | GUI Version (v2.0) |
|---------|---------------------|-------------------|
| Pi file parsing | ✓ | ✓ |
| GFF3 file parsing | ✓ | ✓ |
| Three significance criteria | ✓ | ✓ |
| Smart label layout | ✓ | ✓ |
| PNG output | ✓ | ✓ |
| SVG output (vector) | ✗ | ✓ |
| PDF output | ✗ | ✓ |
| Graphical interface | ✗ | ✓ |
| Real-time logging | ✗ | ✓ |
| Progress indication | ✗ | ✓ |
| Batch processing | ✓ (script) | ✗ |

**Recommendations**:
- Batch process multiple files: Use CLI version, write batch scripts
- Single analysis and visualization: Use GUI version, more intuitive and convenient

## Technical Details

### Vector Graphics Advantages

SVG and PDF formats are vector-based graphic formats with the following advantages:

1. **Lossless scaling**: Can zoom in arbitrarily without distortion
2. **Small file size**: Compared to bitmaps of equivalent resolution
3. **Editable**: Can be modified in vector graphics editing software
4. **High print quality**: Suitable for high-quality printing
5. **Publication-ready**: Academic journals typically require vector formats

### Code Structure

- `PiPlotProcessor`: Core data processing and plotting class
  - `parse_pi_file`: Parse Pi file
  - `parse_gff3_file`: Parse GFF3 file
  - `find_significant_genes`: Find significant genes
  - `plot_pi_with_genes`: Plot chart

- `PiPlotGUI`: Graphical user interface class
  - `setup_ui`: Setup interface layout
  - `browse_*_file`: File browsing functionality
  - `start_analysis`: Start analysis
  - `run_analysis`: Run analysis (background thread)

## Citation

If you use this program in your research, please cite:

```
PiPlot GUI: A Python Tool for Visualizing Pi Values in Genome Regions
https://github.com/Zmy0912/Bioinfo/tree/master/PiPlot
```

## License

This program is licensed under the MIT License. See LICENSE file in the original project for details.

## Changelog

### v2.0 (2026-03-07)
- Added graphical user interface
- Added SVG (vector graphics) format output
- Added PDF format output
- Support for generating multiple formats simultaneously
- Added real-time logging and progress indication
- Improved user experience

### v1.1.0 (2026-01-04)
- Added lower threshold parameter
- Optimized gene sorting logic
- Support for multiple significance type combinations
- Updated color coding system

### v1.0.0 (2026-01-03)
- Initial release
- Support for relative and absolute significance detection
- Smart label layout algorithm

## Contact

For questions or suggestions, please contact via:
- Submit an Issue
- Send email to myzhang0726@foxmail.com

## Acknowledgments

Thanks to the following open source projects:
- pandas
- matplotlib
- numpy
- tkinter
