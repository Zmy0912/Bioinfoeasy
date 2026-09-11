# Time Scale Generator

A powerful graphical time scale generation tool, especially suitable for creating geological time scales in geological research.

## 🚀 Quick Start

### Install Dependencies

```bash
pip install -r requirements.txt
```

### Run the Program

```bash
python time_scale_generator.py
```

### System Requirements

- **Python**: 3.6 or higher
- **Operating System**: Windows / macOS / Linux

## ✨ Main Features

### 1. Time Range Settings
- Support custom minimum and maximum time
- Multiple time unit support:
  - Ma (Million Years)
  - ka (Thousand Years)
  - yr (Years)
  - Ga (Billion Years)
  - Custom time unit

### 2. Scale Configuration
- **Independent scale direction control**
  - Main axis direction: top, bottom, left, right
  - Tick line direction: outward, inward, up, down, left, right
- Adjustable tick interval
- Adjustable tick length (major and minor ticks)
- **Time axis reversal**: Support starting from opposite direction with 0

### 3. Display Options
- Show/hide time label numbers
- Show/hide geological time scale
- Show/hide geological era end times
- Show/hide scale end times
- Show/hide time unit labels
- **Trim trailing zeros**: Smart optimization for number display
- Geological time type selection:
  - Era: Cenozoic, Mesozoic, Paleozoic, Precambrian
  - Period: Detailed geological period division

### 4. Color Customization
- Line color
- Tick color
- Label color
- Background color

### 5. Size Adjustment
- Scale length
- Scale thickness

### 6. Quick Presets
- Cenozoic (0-66 Ma)
- Mesozoic (66-252 Ma)
- Paleozoic (252-541 Ma)
- Phanerozoic (0-541 Ma)

### 7. Element Editing
- Click to select geological era blocks
- Edit element visibility, color, font size, etc.
- Each element can be independently controlled

### 8. Export Functions
- **PDF format**: High-quality vector format
- **SVG format**: Infinitely scalable vector format
- **PNG format**: Easy-to-share bitmap format

### 9. Real-time Preview
- All parameter modifications are reflected in the preview area in real-time
- Automatic updates, no manual refresh required
- Provides "Refresh Preview" button for manual update

### 10. Time Display Precision
- **4 decimal places precision**: Meet precise time labeling requirements
- **Smart formatting**: Optional trimming of trailing zeros after decimal point
- **Full unit names**: Display complete English names and abbreviations

## 📖 Usage

### Basic Operation Flow

1. **Set Time Range**
   - Enter minimum time and maximum time in the "Time Range" area
   - Select appropriate time unit

2. **Configure Scale**
   - Adjust tick interval in "Tick Settings"
   - Use sliders to adjust tick length
   - Select main axis direction and tick line direction in "Scale Direction"
   - Optionally check "Reverse Time Axis (0 starts from opposite direction)"

3. **Display Options**
   - Check "Show Time Labels" to show/hide time numbers
   - Check "Show Geological Scale" to show/hide geological eras
   - Check "Show Era End Times" to show geological era end labels
   - Check "Show Scale End Times" to show scale start/end times
   - Check "Show Time Unit Label" to show complete unit name
   - Check "Trim trailing zeros" to optimize number display
   - Select geological time type (Era or Period)

4. **Customize Appearance**
   - Click color buttons to modify element colors
   - Adjust scale length and thickness

5. **Use Presets**
   - Click preset buttons to quickly load common geological time ranges

6. **Edit Elements**
   - Click to select geological era blocks in the preview area
   - Click "Edit Selected Element" button
   - Modify element properties in the popup dialog

7. **Export Results**
   - Click export buttons to select output format (PDF/SVG/PNG)
   - Choose save location and file name
   - Confirm export

## Technical Features

### Time Unit Conversion
The program internally uses million years (Mya) as the base unit, supporting automatic conversion between multiple time units:
- Ma (Million Years) ← 1 Ma = 1 Mya
- ka (Thousand Years) ← 1 ka = 0.001 Mya
- yr (Years) ← 1 yr = 0.000001 Mya
- Ga (Billion Years) ← 1 Ga = 1000 Mya

### Time Display Precision
- Default display of 4 decimal places precision
- Support smart trimming of trailing zeros after decimal point
- Time unit labels display complete English names and abbreviations
- All export formats (PDF/SVG/PNG) maintain consistent display precision

### Independent Direction Control
- Scale main axis direction and tick line direction are completely independent
- Support main axis position selection: top, bottom, left, right
- Support tick line direction selection: outward, inward, up, down, left, right
- Support time axis reversal: 0 can start from any direction

### Geological Time Data
The program includes complete geological time data:
- Four major geological eras: Cenozoic, Mesozoic, Paleozoic, Precambrian
- Detailed geological periods: Complete division from Quaternary to Hadean
- Each era has corresponding color identification

### Element System
Each rendered element is an independent object, supporting:
- Independent visibility control
- Custom colors
- Font size adjustment
- Position and size editing

## System Requirements

- Python 3.6+
- PyQt5
- reportlab
- svgwrite
- Pillow

## 🛠️ Common Issues

### Program Won't Start

**Issue 1: Prompt "Cannot import PyQt5"**
- Solution: Install PyQt5
```bash
pip install PyQt5
```

**Issue 2: Prompt "Cannot import reportlab"**
- Solution: Install reportlab
```bash
pip install reportlab
```

**Issue 3: Prompt "Cannot import svgwrite"**
- Solution: Install svgwrite
```bash
pip install svgwrite
```

**Issue 4: Prompt "Cannot import Pillow"**
- Solution: Install Pillow
```bash
pip install Pillow
```

### Export Function Issues

**Issue: PDF or SVG export failed**
- Solution: Ensure all dependency packages are installed
- Ensure write permission
- Check disk space

### Display Issues

**Issue: Preview area display abnormal**
- Solution:
  1. Check if time range settings are reasonable
  2. Adjust scale length and thickness
  3. Restart program

### Number Display Issues

**Issue: Time values display too many decimals**
- Solution: Check "Trim trailing zeros" option

**Issue: Time values need to display full precision**
- Solution: Uncheck "Trim trailing zeros" option

## 💡 Usage Tips

1. **Set Tick Interval Reasonably**
   - Use large intervals for large time ranges (e.g., 100 Ma)
   - Use small intervals for small time ranges (e.g., 1 Ma)

2. **Effectively Use Presets**
   - Preset buttons can quickly load common geological time ranges
   - Can make fine adjustments on this basis

3. **Element Editing**
   - Clicking geological era blocks allows individual editing of their properties
   - Can hide unwanted geological eras

4. **Export Format Selection**
   - PDF: Suitable for printing and high-quality documents
   - SVG: Vector format, infinitely scalable
   - PNG: Bitmap format, easy to share and view

5. **Direction Control**
   - Using "outward" tick direction provides best visual effect
   - Time axis reversal can be used for special scenario time axis display

6. **Number Display Optimization**
   - Checking "Trim trailing zeros" makes integer display more concise
   - Can uncheck to maintain full precision when accurate comparison is needed

## Notes

1. Time range inputs must be valid numbers
2. Minimum time must be less than maximum time
3. Tick interval cannot be 0
4. Ensure write permission when exporting files
5. Large time ranges may require adjusting tick interval for better display effect

## Troubleshooting

### Import Errors
If you encounter module import errors, please confirm all dependencies are correctly installed:
```bash
pip install -r requirements.txt --upgrade
```

### Display Issues
If the preview area displays abnormally, try:
1. Adjust scale length and thickness
2. Check if time range settings are reasonable
3. Restart program

### Export Failed
If export function fails, please check:
1. Whether target path has write permission
2. Whether disk space is sufficient
3. Whether filename contains special characters

### Encoding Issues
Windows users may see encoding warnings in the console, this is normal and does not affect program use.

---

## 📝 Dependency Package Description

- **PyQt5**: For building graphical user interface
- **reportlab**: For PDF export functionality
- **svgwrite**: For SVG export functionality
- **Pillow**: For PNG export functionality

## Future Improvement Directions

- Add more geological era presets
- Support custom geological time data
- Add scale style templates
- Support batch export
- Add Undo/Redo functionality
- Support saving and loading project files

## License

This project is for learning and research purposes only.

## Contact

For questions or suggestions, please contact via:
- Create Issue
- Send email feedback

---

**Enjoy using Time Scale Generator!**
