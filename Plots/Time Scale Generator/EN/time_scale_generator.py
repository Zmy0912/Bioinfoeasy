#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Time Scale Generator
Supports geological time scale generation with multiple time units and multi-direction output
"""

import sys
import io

# 设置标准输出编码为UTF-8（修复Windows控制台编码问题）
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout,
                             QHBoxLayout, QLabel, QLineEdit, QPushButton,
                             QComboBox, QCheckBox, QGroupBox, QSlider,
                             QColorDialog, QSplitter, QFileDialog,
                             QMessageBox, QScrollArea, QSpinBox, QDialog)
from PyQt5.QtCore import Qt, QTimer, pyqtSignal, QPointF, QRectF
from PyQt5.QtGui import QPainter, QPen, QColor, QFont, QBrush

# Geological time scale data (in millions of years ago)
GEOLOGICAL_ERAS = [
    {"name": "Cenozoic", "start": 0, "end": 66, "color": "#90EE90"},
    {"name": "Mesozoic", "start": 66, "end": 252, "color": "#FFD700"},
    {"name": "Paleozoic", "start": 252, "end": 541, "color": "#FFA07A"},
    {"name": "Precambrian", "start": 541, "end": 4600, "color": "#DDA0DD"},
]

# Geological period data
GEOLOGICAL_PERIODS = [
    {"name": "Quaternary", "start": 0, "end": 2.588, "color": "#98FB98", "era": "Cenozoic"},
    {"name": "Neogene", "start": 2.588, "end": 23.03, "color": "#8FBC8F", "era": "Cenozoic"},
    {"name": "Paleogene", "start": 23.03, "end": 66, "color": "#7CFC00", "era": "Cenozoic"},
    {"name": "Cretaceous", "start": 66, "end": 145, "color": "#F0E68C", "era": "Mesozoic"},
    {"name": "Jurassic", "start": 145, "end": 201, "color": "#BDB76B", "era": "Mesozoic"},
    {"name": "Triassic", "start": 201, "end": 252, "color": "#DAA520", "era": "Mesozoic"},
    {"name": "Permian", "start": 252, "end": 299, "color": "#FF6347", "era": "Paleozoic"},
    {"name": "Carboniferous", "start": 299, "end": 359, "color": "#FF4500", "era": "Paleozoic"},
    {"name": "Devonian", "start": 359, "end": 419, "color": "#DC143C", "era": "Paleozoic"},
    {"name": "Silurian", "start": 419, "end": 444, "color": "#FF7F50", "era": "Paleozoic"},
    {"name": "Ordovician", "start": 444, "end": 485, "color": "#FFDAB9", "era": "Paleozoic"},
    {"name": "Cambrian", "start": 485, "end": 541, "color": "#FFE4B5", "era": "Paleozoic"},
    {"name": "Ediacaran", "start": 541, "end": 635, "color": "#D8BFD8", "era": "Precambrian"},
    {"name": "Cryogenian", "start": 635, "end": 850, "color": "#DA70D6", "era": "Precambrian"},
    {"name": "Tonian", "start": 850, "end": 1000, "color": "#BA55D3", "era": "Precambrian"},
    {"name": "Stenian", "start": 1000, "end": 1200, "color": "#9370DB", "era": "Precambrian"},
    {"name": "Orosirian", "start": 1200, "end": 1400, "color": "#8A2BE2", "era": "Precambrian"},
    {"name": "Rhyacian", "start": 1400, "end": 1600, "color": "#800080", "era": "Precambrian"},
    {"name": "Calymmian", "start": 1600, "end": 1800, "color": "#7B68EE", "era": "Precambrian"},
    {"name": "Statherian", "start": 1800, "end": 2050, "color": "#6A5ACD", "era": "Precambrian"},
    {"name": "Orosirian", "start": 2050, "end": 2300, "color": "#483D8B", "era": "Precambrian"},
    {"name": "Siderian", "start": 2300, "end": 2500, "color": "#4169E1", "era": "Precambrian"},
    {"name": "Archean", "start": 2500, "end": 4000, "color": "#6495ED", "era": "Precambrian"},
    {"name": "Hadean", "start": 4000, "end": 4600, "color": "#87CEEB", "era": "Precambrian"},
]

class TimeScaleElement:
    """Time scale element class, supports segmented editing"""
    def __init__(self, element_type, position, size, text="", color=QColor(0, 0, 0)):
        self.element_type = element_type  # 'line', 'tick', 'label', 'era_block'
        self.position = position  # QPointF
        self.size = size  # QSizeF
        self.text = text
        self.color = color
        self.visible = True
        self.font_size = 10
        self.font_family = "Arial"
        self.font_bold = False

class TimeScaleRenderer(QWidget):
    """Time scale rendering component"""
    
    elementsChanged = pyqtSignal()
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.elements = []
        self.min_time = 0
        self.max_time = 100
        self.time_unit = "Ma"  # Ma (million years), ka (thousand years), yr (years), Ga (billion years)
        # Time unit name mapping
        self.time_unit_names = {
            "Ma": "Million years (Ma)",
            "ka": "Thousand years (ka)",
            "yr": "Years (yr)",
            "Ga": "Billion years (Ga)"
        }
        self.start_direction = "bottom"  # Scale main axis direction: top, bottom, left, right
        self.tick_direction = "outward"  # Tick line direction: inward, outward, up, down, left, right
        self.tick_interval = 10
        self.show_labels = True
        self.show_geological = False
        self.show_era_end_labels = True  # Show geological period end time
        self.show_scale_end_time = True  # Show scale end time
        self.show_time_unit_label = True  # Show time unit label
        self.reverse_time_axis = False  # Reverse time axis (0 starts from opposite direction)
        self.trim_trailing_zeros = False  # Omit trailing zeros after decimal point
        self.scale_length = 800
        self.scale_thickness = 50
        self.tick_length = 20
        self.major_tick_length = 30
        self.label_offset = 15
        self.line_color = QColor(0, 0, 0)
        self.tick_color = QColor(0, 0, 0)
        self.label_color = QColor(0, 0, 0)
        self.background_color = QColor(255, 255, 255)
        self.selected_element = None
        self.geological_scale_type = "era"  # era, period, or both
        self.custom_elements = []

        self.setMouseTracking(True)
        self.setMinimumSize(800, 600)
        self.generate_scale()

    def convert_to_mya(self, value):
        """Convert different time units to millions of years ago"""
        conversions = {
            "Ma": 1.0,
            "ka": 0.001,
            "yr": 0.000001,
            "Ga": 1000.0
        }
        return value * conversions.get(self.time_unit, 1.0)

    def convert_from_mya(self, value_mya):
        """Convert millions of years ago to current time unit"""
        conversions = {
            "Ma": 1.0,
            "ka": 1000.0,
            "yr": 1000000.0,
            "Ga": 0.001
        }
        return value_mya * conversions.get(self.time_unit, 1.0)

    def format_time_value(self, value):
        """Format time value, decide whether to omit trailing zeros based on options"""
        if self.trim_trailing_zeros:
            # Omit trailing zeros after decimal point
            if value == int(value):
                return f"{int(value)}"
            else:
                # Keep up to 4 decimal places, but omit trailing zeros
                formatted = f"{value:.4f}"
                # Remove trailing zeros
                formatted = formatted.rstrip('0').rstrip('.') if '.' in formatted else formatted
                return formatted
        else:
            # Always show 4 decimal places
            return f"{value:.4f}"

    def generate_scale(self):
        """Generate time scale elements"""
        self.elements = []

        # Calculate main axis position based on scale direction
        if self.start_direction == "top":
            # Main axis at top
            main_line_start = QPointF(50, 300)
            main_line_end = QPointF(50 + self.scale_length, 300)
            scale_center_y = 300
            scale_center_x = 50 + self.scale_length / 2
        elif self.start_direction == "bottom":
            # Main axis at bottom
            main_line_start = QPointF(50, 300)
            main_line_end = QPointF(50 + self.scale_length, 300)
            scale_center_y = 300
            scale_center_x = 50 + self.scale_length / 2
        elif self.start_direction == "left":
            # Main axis at left
            main_line_start = QPointF(400, 50)
            main_line_end = QPointF(400, 50 + self.scale_length)
            scale_center_x = 400
            scale_center_y = 50 + self.scale_length / 2
        else:  # right
            # Main axis at right
            main_line_start = QPointF(400, 50)
            main_line_end = QPointF(400, 50 + self.scale_length)
            scale_center_x = 400
            scale_center_y = 50 + self.scale_length / 2

        # Main scale line
        main_line = TimeScaleElement('line', main_line_start, main_line_end)
        self.elements.append(main_line)

        # Calculate number of ticks
        min_mya = self.convert_to_mya(self.min_time)
        max_mya = self.convert_to_mya(self.max_time)
        tick_interval_mya = self.convert_to_mya(self.tick_interval)

        if tick_interval_mya == 0:
            tick_interval_mya = 1.0

        total_ticks = int((max_mya - min_mya) / tick_interval_mya) + 1

        for i in range(total_ticks):
            tick_value_mya = min_mya + i * tick_interval_mya
            tick_value_current = self.convert_from_mya(tick_value_mya)

            # Calculate tick position on main axis
            progress = (tick_value_mya - min_mya) / (max_mya - min_mya) if (max_mya - min_mya) > 0 else 0

            # Reverse time axis direction
            if self.reverse_time_axis:
                progress = 1.0 - progress

            # Check if it's a major tick
            is_major = i % 5 == 0
            tick_len = self.major_tick_length if is_major else self.tick_length

            # Calculate tick line based on scale direction and tick direction
            if self.start_direction in ["top", "bottom"]:
                # Horizontal main axis
                tick_x = 50 + progress * self.scale_length
                tick_y = scale_center_y

                # Determine tick line end point based on tick direction
                if self.tick_direction == "outward":
                    # Outward (relative to scale main axis)
                    tick_y_end = tick_y - tick_len if self.start_direction == "top" else tick_y + tick_len
                elif self.tick_direction == "inward":
                    # Inward
                    tick_y_end = tick_y + tick_len if self.start_direction == "top" else tick_y - tick_len
                elif self.tick_direction == "up":
                    tick_y_end = tick_y - tick_len
                elif self.tick_direction == "down":
                    tick_y_end = tick_y + tick_len
                else:  # left, right - default downward for horizontal scale
                    tick_y_end = tick_y + tick_len

                # Tick line
                tick = TimeScaleElement('tick', QPointF(tick_x, tick_y),
                                       QPointF(tick_x, tick_y_end),
                                       str(tick_value_current), self.tick_color)
                tick.is_major = is_major
                self.elements.append(tick)

                # Label
                if self.show_labels:
                    label_y = tick_y_end - self.label_offset if self.tick_direction == "up" else tick_y_end + self.label_offset
                    label = TimeScaleElement('label', QPointF(tick_x, label_y),
                                           QPointF(tick_x, label_y),
                                           self.format_time_value(tick_value_current), self.label_color)
                    label.is_major = is_major
                    self.elements.append(label)

            else:  # left, right - vertical main axis
                tick_y = 50 + progress * self.scale_length
                tick_x = scale_center_x

                # Determine tick line end point based on tick direction
                if self.tick_direction == "outward":
                    # Outward (relative to scale main axis)
                    tick_x_end = tick_x - tick_len if self.start_direction == "left" else tick_x + tick_len
                elif self.tick_direction == "inward":
                    # Inward
                    tick_x_end = tick_x + tick_len if self.start_direction == "left" else tick_x - tick_len
                elif self.tick_direction == "left":
                    tick_x_end = tick_x - tick_len
                elif self.tick_direction == "right":
                    tick_x_end = tick_x + tick_len
                else:  # up, down - default rightward for vertical scale
                    tick_x_end = tick_x + tick_len

                # Tick line
                tick = TimeScaleElement('tick', QPointF(tick_x, tick_y),
                                       QPointF(tick_x_end, tick_y),
                                       str(tick_value_current), self.tick_color)
                tick.is_major = is_major
                self.elements.append(tick)

                # Label
                if self.show_labels:
                    label_x = tick_x_end - self.label_offset if self.tick_direction == "left" else tick_x_end + self.label_offset
                    label = TimeScaleElement('label', QPointF(label_x, tick_y),
                                           QPointF(label_x, tick_y),
                                           self.format_time_value(tick_value_current), self.label_color)
                    label.is_major = is_major
                    self.elements.append(label)

        # Add scale start and end time labels
        if self.show_scale_end_time:
            if self.start_direction in ["top", "bottom"]:
                # Horizontal scale
                # Start time
                start_value = self.max_time if self.reverse_time_axis else self.min_time
                start_label = TimeScaleElement('label',
                                              QPointF(50 - 60, scale_center_y - 20),
                                              QPointF(50 - 60, scale_center_y - 20),
                                              self.format_time_value(start_value), QColor(0, 0, 100))
                self.elements.append(start_label)

                # End time
                end_value = self.min_time if self.reverse_time_axis else self.max_time
                end_label = TimeScaleElement('label',
                                            QPointF(50 + self.scale_length + 60, scale_center_y - 20),
                                            QPointF(50 + self.scale_length + 60, scale_center_y - 20),
                                            self.format_time_value(end_value), QColor(0, 0, 100))
                self.elements.append(end_label)
            else:  # left, right - vertical scale
                # Start time
                start_value = self.max_time if self.reverse_time_axis else self.min_time
                start_label = TimeScaleElement('label',
                                              QPointF(scale_center_x - 80, 50 - 15),
                                              QPointF(scale_center_x - 80, 50 - 15),
                                              self.format_time_value(start_value), QColor(0, 0, 100))
                self.elements.append(start_label)

                # End time
                end_value = self.min_time if self.reverse_time_axis else self.max_time
                end_label = TimeScaleElement('label',
                                            QPointF(scale_center_x - 80, 50 + self.scale_length + 15),
                                            QPointF(scale_center_x - 80, 50 + self.scale_length + 15),
                                            self.format_time_value(end_value), QColor(0, 0, 100))
                self.elements.append(end_label)

        # Add time unit label
        if self.show_time_unit_label:
            # Get full unit name
            unit_display_name = self.time_unit_names.get(self.time_unit, self.time_unit)
            unit_label = TimeScaleElement('label',
                                         QPointF(scale_center_x, scale_center_y - 60) if self.start_direction in ["top", "bottom"] else QPointF(scale_center_x + 100, scale_center_y),
                                         QPointF(scale_center_x, scale_center_y - 60) if self.start_direction in ["top", "bottom"] else QPointF(scale_center_x + 100, scale_center_y),
                                         f"Unit: {unit_display_name}", QColor(0, 0, 139))
            unit_label.font_size = 12
            unit_label.font_bold = True
            self.elements.append(unit_label)

        # Geological time scale
        if self.show_geological:
            self.generate_geological_scale(min_mya, max_mya)

        self.elementsChanged.emit()
        self.update()

    def generate_geological_scale(self, min_mya, max_mya):
        """Generate geological time scale"""
        # Determine whether to show eras and periods based on geological time scale type
        if self.geological_scale_type == "era":
            show_era = True
            show_period = False
        elif self.geological_scale_type == "period":
            show_era = False
            show_period = True
        else:  # both
            show_era = True
            show_period = True

        # Determine geological time scale position based on scale direction
        if self.start_direction == "top":
            base_scale_pos = 300
            is_horizontal = True
            scale_start = 50
        elif self.start_direction == "bottom":
            base_scale_pos = 300 + self.scale_thickness
            is_horizontal = True
            scale_start = 50
        elif self.start_direction == "left":
            base_scale_pos = 400
            is_horizontal = False
            scale_start = 50
        else:  # right
            base_scale_pos = 400 + self.scale_thickness
            is_horizontal = False
            scale_start = 50

        # Show eras (outer layer)
        if show_era:
            era_scale_size = 30
            if is_horizontal:
                scale_pos = base_scale_pos + 10
            else:
                scale_pos = base_scale_pos + 10
            
            self._generate_single_layer(GEOLOGICAL_ERAS, min_mya, max_mya, scale_pos,
                                        scale_start, era_scale_size, is_horizontal, is_era=True)

        # Show periods (inner layer, inside eras)
        if show_period:
            period_scale_size = 25
            if show_era:
                # If showing both eras and periods, periods are outside eras
                if is_horizontal:
                    scale_pos = base_scale_pos + 30 + 10 + 5  # era height + spacing
                else:
                    scale_pos = base_scale_pos + 30 + 10 + 5
            else:
                # If only showing periods
                if is_horizontal:
                    scale_pos = base_scale_pos + 10
                else:
                    scale_pos = base_scale_pos + 10
            
            self._generate_single_layer(GEOLOGICAL_PERIODS, min_mya, max_mya, scale_pos,
                                        scale_start, period_scale_size, is_horizontal, is_era=False)

    def _generate_single_layer(self, geo_data, min_mya, max_mya, scale_pos, scale_start, scale_size, is_horizontal, is_era):
        """Generate single layer geological time scale"""
        for period in geo_data:
            # Check if within range
            if period["end"] < min_mya or period["start"] > max_mya:
                continue

            # Calculate display range
            start_pos = max(period["start"], min_mya)
            end_pos = min(period["end"], max_mya)

            progress_start = (start_pos - min_mya) / (max_mya - min_mya) if (max_mya - min_mya) > 0 else 0
            progress_end = (end_pos - min_mya) / (max_mya - min_mya) if (max_mya - min_mya) > 0 else 0

            # Calculate position and size of geological time block
            if is_horizontal:
                x1 = scale_start + progress_start * self.scale_length
                y1 = scale_pos
                block_length = (progress_end - progress_start) * self.scale_length
                block_width = block_length
                block_height = scale_size

                era_block = TimeScaleElement('era_block', QPointF(x1, y1),
                                            QPointF(block_width, block_height),
                                            period["name"], QColor(period["color"]))
                era_block.start_time = period["start"]
                era_block.end_time = period["end"]
                era_block.is_era = is_era  # Mark whether it's an era
                self.elements.append(era_block)

                # Show end time label
                if self.show_era_end_labels and end_pos <= max_mya:
                    end_x = scale_start + progress_end * self.scale_length
                    label_y = scale_pos + scale_size + 12
                    label = TimeScaleElement('label', QPointF(end_x, label_y),
                                           QPointF(end_x, label_y),
                                           f"{period['end']:.1f} Ma", QColor(0, 0, 0))
                    label.is_major = True
                    self.elements.append(label)
            else:  # vertical
                y1 = scale_start + progress_start * self.scale_length
                x1 = scale_pos
                block_height = (progress_end - progress_start) * self.scale_length
                block_width = scale_size

                era_block = TimeScaleElement('era_block', QPointF(x1, y1),
                                            QPointF(block_width, block_height),
                                            period["name"], QColor(period["color"]))
                era_block.start_time = period["start"]
                era_block.end_time = period["end"]
                era_block.is_era = is_era  # Mark whether it's an era
                self.elements.append(era_block)

                # Show end time label
                if self.show_era_end_labels and end_pos <= max_mya:
                    end_y = scale_start + progress_end * self.scale_length
                    label_x = scale_pos + scale_size + 30
                    label = TimeScaleElement('label', QPointF(label_x, end_y),
                                           QPointF(label_x, end_y),
                                           f"{period['end']:.1f} Ma", QColor(0, 0, 0))
                    label.is_major = True
                    self.elements.append(label)

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        super().paintEvent(event)
        
        # Draw background
        painter.fillRect(self.rect(), self.background_color)
        
        # Draw all elements
        for element in reversed(self.elements):  # Reverse traverse, ensure labels are on top
            if not element.visible:
                continue

            if element.element_type == 'line':
                pen = QPen(element.color, 2)
                painter.setPen(pen)
                painter.drawLine(element.position, element.size)

            elif element.element_type == 'tick':
                pen = QPen(element.color, 1 if not getattr(element, 'is_major', False) else 2)
                painter.setPen(pen)
                painter.drawLine(element.position, element.size)
            
            elif element.element_type == 'label':
                painter.setPen(element.color)
                font = QFont(element.font_family, element.font_size)
                font.setBold(element.font_bold)
                painter.setFont(font)

                # Calculate rectangle size based on text length
                font_metrics = painter.fontMetrics()
                text_width = font_metrics.width(element.text)
                text_height = font_metrics.height()

                # Horizontal scale label: adaptive text width, fixed height
                if self.start_direction in ["top", "bottom"]:
                    text_rect = QRectF(element.position.x() - text_width/2 - 10,
                                     element.position.y() - text_height/2,
                                     text_width + 20,
                                     text_height + 5)
                else:
                    # Vertical scale label
                    text_rect = QRectF(element.position.x() - text_width/2 - 10,
                                     element.position.y() - text_height/2,
                                     text_width + 20,
                                     text_height + 5)

                painter.drawText(text_rect, Qt.AlignCenter, element.text)
            
            elif element.element_type == 'era_block':
                painter.setBrush(QBrush(element.color))
                painter.setPen(QPen(Qt.black, 1))
                rect = QRectF(element.position.x(), element.position.y(),
                             element.size.x(), element.size.y())
                painter.drawRect(rect)
                
                # Draw period name
                painter.setPen(Qt.black)
                painter.setFont(QFont("Arial", 8))
                text_rect = QRectF(rect)
                painter.drawText(text_rect, Qt.AlignCenter, element.text)

    def mousePressEvent(self, event):
        """Mouse click to select element"""
        pos = event.pos()
        for element in reversed(self.elements):
            if not element.visible:
                continue
                
            if element.element_type == 'era_block':
                rect = QRectF(element.position.x(), element.position.y(),
                             element.size.x(), element.size.y())
                if rect.contains(pos):
                    self.selected_element = element
                    self.update()
                    return

class ElementEditDialog(QDialog):
    """Element edit dialog"""
    
    def __init__(self, element, parent=None):
        super().__init__(parent)
        self.element = element
        self.init_ui()
    
    def init_ui(self):
        layout = QVBoxLayout()
        
        # Visibility
        visible_layout = QHBoxLayout()
        visible_label = QLabel("Visible:")
        visible_checkbox = QCheckBox("Show")
        visible_checkbox.setChecked(self.element.visible)
        visible_checkbox.toggled.connect(self.set_visible)
        visible_layout.addWidget(visible_label)
        visible_layout.addWidget(visible_checkbox)
        layout.addLayout(visible_layout)
        
        # Color
        color_layout = QHBoxLayout()
        color_label = QLabel("Color:")
        color_button = QPushButton("Choose Color")
        color_button.clicked.connect(self.choose_color)
        color_layout.addWidget(color_label)
        color_layout.addWidget(color_button)
        layout.addLayout(color_layout)
        
        # Font size
        if self.element.element_type in ['label', 'era_block']:
            font_layout = QHBoxLayout()
            font_label = QLabel("Font Size:")
            font_spin = QSpinBox()
            font_spin.setRange(6, 72)
            font_spin.setValue(self.element.font_size)
            font_spin.valueChanged.connect(self.set_font_size)
            font_layout.addWidget(font_label)
            font_layout.addWidget(font_spin)
            layout.addLayout(font_layout)
        
        # OK and Cancel buttons
        button_layout = QHBoxLayout()
        ok_button = QPushButton("OK")
        ok_button.clicked.connect(self.accept)
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        button_layout.addWidget(ok_button)
        button_layout.addWidget(cancel_button)
        layout.addLayout(button_layout)
        
        self.setLayout(layout)
        self.setWindowTitle("Edit Element")
    
    def set_visible(self, visible):
        self.element.visible = visible
    
    def choose_color(self):
        color = QColorDialog.getColor(self.element.color, self)
        if color.isValid():
            self.element.color = color
    
    def set_font_size(self, size):
        self.element.font_size = size

class TimeScaleGenerator(QMainWindow):
    """Time scale generator main window"""
    
    def __init__(self):
        super().__init__()
        self.init_ui()
        
        # Real-time update timer
        self.update_timer = QTimer()
        self.update_timer.timeout.connect(self.update_preview)
        self.update_timer.start(300)  # Update every 300ms

    def init_ui(self):
        self.setWindowTitle("Time Scale Generator")
        self.setGeometry(100, 100, 1400, 900)
        
        # Main window widget
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        
        # Main layout
        main_layout = QHBoxLayout()
        main_widget.setLayout(main_layout)
        
        # Left control panel
        control_panel = self.create_control_panel()
        
        # Right preview area
        preview_panel = self.create_preview_panel()
        
        # Use splitter
        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(control_panel)
        splitter.addWidget(preview_panel)
        splitter.setStretchFactor(1, 2)
        
        main_layout.addWidget(splitter)
    
    def create_control_panel(self):
        """Create control panel"""
        panel = QScrollArea()
        panel.setWidgetResizable(True)
        panel.setFixedWidth(400)
        
        container = QWidget()
        layout = QVBoxLayout()
        container.setLayout(layout)
        
        # Time range settings
        time_group = QGroupBox("Time Range")
        time_layout = QVBoxLayout()
        
        min_time_layout = QHBoxLayout()
        min_time_layout.addWidget(QLabel("Min Time:"))
        self.min_time_input = QLineEdit("0")
        self.min_time_input.textChanged.connect(self.on_parameter_change)
        min_time_layout.addWidget(self.min_time_input)
        time_layout.addLayout(min_time_layout)
        
        max_time_layout = QHBoxLayout()
        max_time_layout.addWidget(QLabel("Max Time:"))
        self.max_time_input = QLineEdit("100")
        self.max_time_input.textChanged.connect(self.on_parameter_change)
        max_time_layout.addWidget(self.max_time_input)
        time_layout.addLayout(max_time_layout)
        
        time_group.setLayout(time_layout)
        layout.addWidget(time_group)
        
        # Time unit settings
        unit_group = QGroupBox("Time Unit")
        unit_layout = QVBoxLayout()
        
        self.time_unit_combo = QComboBox()
        self.time_unit_combo.addItems(["Ma (Million Years)", "ka (Thousand Years)", "yr (Years)", "Ga (Billion Years)", "Custom"])
        self.time_unit_combo.currentTextChanged.connect(self.on_unit_change)
        unit_layout.addWidget(self.time_unit_combo)
        
        self.custom_unit_input = QLineEdit()
        self.custom_unit_input.setPlaceholderText("Enter custom unit name")
        self.custom_unit_input.textChanged.connect(self.on_parameter_change)
        unit_layout.addWidget(self.custom_unit_input)
        self.custom_unit_input.hide()
        
        unit_group.setLayout(unit_layout)
        layout.addWidget(unit_group)
        
        # Tick interval settings
        tick_group = QGroupBox("Tick Settings")
        tick_layout = QVBoxLayout()
        
        tick_interval_layout = QHBoxLayout()
        tick_interval_layout.addWidget(QLabel("Tick Interval:"))
        self.tick_interval_input = QLineEdit("10")
        self.tick_interval_input.textChanged.connect(self.on_parameter_change)
        tick_interval_layout.addWidget(self.tick_interval_input)
        tick_layout.addLayout(tick_interval_layout)
        
        tick_length_layout = QHBoxLayout()
        tick_length_layout.addWidget(QLabel("Tick Length:"))
        self.tick_length_slider = QSlider(Qt.Horizontal)
        self.tick_length_slider.setRange(5, 50)
        self.tick_length_slider.setValue(20)
        self.tick_length_slider.valueChanged.connect(self.on_parameter_change)
        tick_length_layout.addWidget(self.tick_length_slider)
        tick_layout.addLayout(tick_length_layout)
        
        major_tick_layout = QHBoxLayout()
        major_tick_layout.addWidget(QLabel("Major Tick Length:"))
        self.major_tick_slider = QSlider(Qt.Horizontal)
        self.major_tick_slider.setRange(10, 80)
        self.major_tick_slider.setValue(30)
        self.major_tick_slider.valueChanged.connect(self.on_parameter_change)
        major_tick_layout.addWidget(self.major_tick_slider)
        tick_layout.addLayout(major_tick_layout)
        
        tick_group.setLayout(tick_layout)
        layout.addWidget(tick_group)
        
        # Direction settings
        direction_group = QGroupBox("Scale Direction")
        direction_layout = QVBoxLayout()

        # Scale main axis direction
        axis_layout = QHBoxLayout()
        axis_layout.addWidget(QLabel("Main Axis Direction:"))
        self.axis_direction_combo = QComboBox()
        self.axis_direction_combo.addItems(["bottom (down)", "top (up)", "left (left)", "right (right)"])
        self.axis_direction_combo.currentTextChanged.connect(self.on_parameter_change)
        axis_layout.addWidget(self.axis_direction_combo)
        direction_layout.addLayout(axis_layout)

        # Tick line direction
        tick_dir_layout = QHBoxLayout()
        tick_dir_layout.addWidget(QLabel("Tick Direction:"))
        self.tick_direction_combo = QComboBox()
        self.tick_direction_combo.addItems(["outward (outward)", "inward (inward)", "up (up)", "down (down)", "left (left)", "right (right)"])
        self.tick_direction_combo.currentTextChanged.connect(self.on_parameter_change)
        tick_dir_layout.addWidget(self.tick_direction_combo)
        direction_layout.addLayout(tick_dir_layout)

        direction_group.setLayout(direction_layout)
        layout.addWidget(direction_group)

        # Display settings
        display_group = QGroupBox("Display Settings")
        display_layout = QVBoxLayout()

        self.show_labels_checkbox = QCheckBox("Show Time Labels")
        self.show_labels_checkbox.setChecked(True)
        self.show_labels_checkbox.toggled.connect(self.on_parameter_change)
        display_layout.addWidget(self.show_labels_checkbox)

        self.show_geological_checkbox = QCheckBox("Show Geological Time Scale")
        self.show_geological_checkbox.toggled.connect(self.on_parameter_change)
        display_layout.addWidget(self.show_geological_checkbox)

        self.show_era_end_labels_checkbox = QCheckBox("Show Period End Time")
        self.show_era_end_labels_checkbox.setChecked(True)
        self.show_era_end_labels_checkbox.toggled.connect(self.on_parameter_change)
        display_layout.addWidget(self.show_era_end_labels_checkbox)

        self.show_scale_end_time_checkbox = QCheckBox("Show Scale End Time")
        self.show_scale_end_time_checkbox.setChecked(True)
        self.show_scale_end_time_checkbox.toggled.connect(self.on_parameter_change)
        display_layout.addWidget(self.show_scale_end_time_checkbox)

        self.show_time_unit_label_checkbox = QCheckBox("Show Time Unit Label")
        self.show_time_unit_label_checkbox.setChecked(True)
        self.show_time_unit_label_checkbox.toggled.connect(self.on_parameter_change)
        display_layout.addWidget(self.show_time_unit_label_checkbox)

        self.reverse_time_axis_checkbox = QCheckBox("Reverse Time Axis (0 starts from opposite direction)")
        self.reverse_time_axis_checkbox.setChecked(False)
        self.reverse_time_axis_checkbox.toggled.connect(self.on_parameter_change)
        display_layout.addWidget(self.reverse_time_axis_checkbox)

        self.trim_trailing_zeros_checkbox = QCheckBox("Omit Trailing Zeros After Decimal Point")
        self.trim_trailing_zeros_checkbox.setChecked(False)
        self.trim_trailing_zeros_checkbox.toggled.connect(self.on_parameter_change)
        display_layout.addWidget(self.trim_trailing_zeros_checkbox)

        geo_type_layout = QHBoxLayout()
        geo_type_layout.addWidget(QLabel("Geological Time Type:"))
        self.geo_type_combo = QComboBox()
        self.geo_type_combo.addItems(["Era", "Period", "Era and Period Both"])
        self.geo_type_combo.currentTextChanged.connect(self.on_parameter_change)
        geo_type_layout.addWidget(self.geo_type_combo)
        display_layout.addLayout(geo_type_layout)
        
        display_group.setLayout(display_layout)
        layout.addWidget(display_group)
        
        # Color settings
        color_group = QGroupBox("Color Settings")
        color_layout = QVBoxLayout()
        
        line_color_layout = QHBoxLayout()
        line_color_layout.addWidget(QLabel("Line Color:"))
        self.line_color_button = QPushButton()
        self.line_color_button.setStyleSheet("background-color: black;")
        self.line_color_button.clicked.connect(lambda: self.choose_color('line'))
        line_color_layout.addWidget(self.line_color_button)
        color_layout.addLayout(line_color_layout)
        
        tick_color_layout = QHBoxLayout()
        tick_color_layout.addWidget(QLabel("Tick Color:"))
        self.tick_color_button = QPushButton()
        self.tick_color_button.setStyleSheet("background-color: black;")
        self.tick_color_button.clicked.connect(lambda: self.choose_color('tick'))
        tick_color_layout.addWidget(self.tick_color_button)
        color_layout.addLayout(tick_color_layout)
        
        label_color_layout = QHBoxLayout()
        label_color_layout.addWidget(QLabel("Label Color:"))
        self.label_color_button = QPushButton()
        self.label_color_button.setStyleSheet("background-color: black;")
        self.label_color_button.clicked.connect(lambda: self.choose_color('label'))
        label_color_layout.addWidget(self.label_color_button)
        color_layout.addLayout(label_color_layout)
        
        background_color_layout = QHBoxLayout()
        background_color_layout.addWidget(QLabel("Background Color:"))
        self.background_color_button = QPushButton()
        self.background_color_button.setStyleSheet("background-color: white;")
        self.background_color_button.clicked.connect(lambda: self.choose_color('background'))
        background_color_layout.addWidget(self.background_color_button)
        color_layout.addLayout(background_color_layout)
        
        color_group.setLayout(color_layout)
        layout.addWidget(color_group)
        
        # Size settings
        size_group = QGroupBox("Size Settings")
        size_layout = QVBoxLayout()
        
        scale_length_layout = QHBoxLayout()
        scale_length_layout.addWidget(QLabel("Scale Length:"))
        self.scale_length_input = QLineEdit("800")
        self.scale_length_input.textChanged.connect(self.on_parameter_change)
        scale_length_layout.addWidget(self.scale_length_input)
        size_layout.addLayout(scale_length_layout)
        
        scale_thickness_layout = QHBoxLayout()
        scale_thickness_layout.addWidget(QLabel("Scale Thickness:"))
        self.scale_thickness_input = QLineEdit("50")
        self.scale_thickness_input.textChanged.connect(self.on_parameter_change)
        scale_thickness_layout.addWidget(self.scale_thickness_input)
        size_layout.addLayout(scale_thickness_layout)
        
        size_group.setLayout(size_layout)
        layout.addWidget(size_group)
        
        # Preset buttons
        preset_group = QGroupBox("Common Presets")
        preset_layout = QVBoxLayout()
        
        preset_button1 = QPushButton("Cenozoic (0-66 Ma)")
        preset_button1.clicked.connect(lambda: self.apply_preset(0, 66, "Ma"))
        preset_layout.addWidget(preset_button1)
        
        preset_button2 = QPushButton("Mesozoic (66-252 Ma)")
        preset_button2.clicked.connect(lambda: self.apply_preset(66, 252, "Ma"))
        preset_layout.addWidget(preset_button2)
        
        preset_button3 = QPushButton("Paleozoic (252-541 Ma)")
        preset_button3.clicked.connect(lambda: self.apply_preset(252, 541, "Ma"))
        preset_layout.addWidget(preset_button3)
        
        preset_button4 = QPushButton("Phanerozoic (0-541 Ma)")
        preset_button4.clicked.connect(lambda: self.apply_preset(0, 541, "Ma"))
        preset_layout.addWidget(preset_button4)
        
        preset_group.setLayout(preset_layout)
        layout.addWidget(preset_group)
        
        # Element editing
        element_group = QGroupBox("Element Editing")
        element_layout = QVBoxLayout()
        
        edit_button = QPushButton("Edit Selected Element")
        edit_button.clicked.connect(self.edit_selected_element)
        element_layout.addWidget(edit_button)
        
        self.element_info_label = QLabel("No element selected")
        self.element_info_label.setWordWrap(True)
        element_layout.addWidget(self.element_info_label)
        
        element_group.setLayout(element_layout)
        layout.addWidget(element_group)

        # Action buttons
        action_group = QGroupBox("Actions")
        action_layout = QVBoxLayout()

        refresh_button = QPushButton("Refresh Preview")
        refresh_button.clicked.connect(self.refresh_preview)
        action_layout.addWidget(refresh_button)

        action_group.setLayout(action_layout)
        layout.addWidget(action_group)

        # Export buttons
        export_group = QGroupBox("Export")
        export_layout = QVBoxLayout()

        export_pdf_button = QPushButton("Export as PDF")
        export_pdf_button.clicked.connect(self.export_pdf)
        export_layout.addWidget(export_pdf_button)

        export_svg_button = QPushButton("Export as SVG")
        export_svg_button.clicked.connect(self.export_svg)
        export_layout.addWidget(export_svg_button)

        export_png_button = QPushButton("Export as PNG")
        export_png_button.clicked.connect(self.export_png)
        export_layout.addWidget(export_png_button)

        export_group.setLayout(export_layout)
        layout.addWidget(export_group)

        layout.addStretch()
        panel.setWidget(container)

        return panel

    def create_preview_panel(self):
        """Create preview panel"""
        container = QWidget()
        layout = QVBoxLayout()
        container.setLayout(layout)
        
        # Title
        title = QLabel("Preview Area")
        title.setFont(QFont("Arial", 14, QFont.Bold))
        layout.addWidget(title)
        
        # Renderer
        self.renderer = TimeScaleRenderer()
        self.renderer.elementsChanged.connect(self.on_elements_changed)
        
        # Wrap in scroll area
        scroll = QScrollArea()
        scroll.setWidget(self.renderer)
        scroll.setWidgetResizable(True)
        
        layout.addWidget(scroll)
        
        # Info bar
        self.info_label = QLabel("Ready")
        layout.addWidget(self.info_label)
        
        return container
    
    def on_parameter_change(self):
        """Handle parameter change"""
        pass  # Auto-triggered by timer
    
    def on_unit_change(self, text):
        """Handle time unit change"""
        if text == "Custom":
            self.custom_unit_input.show()
        else:
            self.custom_unit_input.hide()
    
    def choose_color(self, color_type):
        """Choose color"""
        color = QColorDialog.getColor()
        if color.isValid():
            if color_type == 'line':
                self.renderer.line_color = color
                self.line_color_button.setStyleSheet(f"background-color: {color.name()};")
            elif color_type == 'tick':
                self.renderer.tick_color = color
                self.tick_color_button.setStyleSheet(f"background-color: {color.name()};")
            elif color_type == 'label':
                self.renderer.label_color = color
                self.label_color_button.setStyleSheet(f"background-color: {color.name()};")
            elif color_type == 'background':
                self.renderer.background_color = color
                self.background_color_button.setStyleSheet(f"background-color: {color.name()};")
            self.renderer.update()
    
    def update_preview(self):
        """Update preview"""
        try:
            # Get parameters
            min_time = float(self.min_time_input.text())
            max_time = float(self.max_time_input.text())
            tick_interval = float(self.tick_interval_input.text())

            # Get time unit
            unit_text = self.time_unit_combo.currentText()
            if unit_text == "Custom":
                self.renderer.time_unit = self.custom_unit_input.text() or "Ma"
            else:
                self.renderer.time_unit = unit_text.split()[0]

            # Get main axis direction
            axis_direction_text = self.axis_direction_combo.currentText()
            self.renderer.start_direction = axis_direction_text.split()[0]

            # Get tick direction
            tick_direction_text = self.tick_direction_combo.currentText()
            self.renderer.tick_direction = tick_direction_text.split()[0]

            # Update renderer parameters
            self.renderer.min_time = min_time
            self.renderer.max_time = max_time
            self.renderer.tick_interval = tick_interval
            self.renderer.tick_length = self.tick_length_slider.value()
            self.renderer.major_tick_length = self.major_tick_slider.value()
            self.renderer.show_labels = self.show_labels_checkbox.isChecked()
            self.renderer.show_geological = self.show_geological_checkbox.isChecked()
            self.renderer.show_era_end_labels = self.show_era_end_labels_checkbox.isChecked()
            self.renderer.show_scale_end_time = self.show_scale_end_time_checkbox.isChecked()
            self.renderer.show_time_unit_label = self.show_time_unit_label_checkbox.isChecked()
            self.renderer.reverse_time_axis = self.reverse_time_axis_checkbox.isChecked()
            self.renderer.trim_trailing_zeros = self.trim_trailing_zeros_checkbox.isChecked()
            self.renderer.scale_length = float(self.scale_length_input.text())
            self.renderer.scale_thickness = float(self.scale_thickness_input.text())

            # Geological time type
            geo_type = self.geo_type_combo.currentText()
            if "Both" in geo_type:
                self.renderer.geological_scale_type = "both"
            elif "Era" in geo_type:
                self.renderer.geological_scale_type = "era"
            else:
                self.renderer.geological_scale_type = "period"

            # Regenerate scale
            self.renderer.generate_scale()

        except ValueError:
            pass

    def refresh_preview(self):
        """Manually refresh preview"""
        self.update_preview()
        self.info_label.setText("Refreshed")

    def on_elements_changed(self):
        """Update info when elements change"""
        element_count = len(self.renderer.elements)
        self.info_label.setText(f"Currently has {element_count} elements")

    def edit_selected_element(self):
        """Edit selected element"""
        if self.renderer.selected_element:
            dialog = ElementEditDialog(self.renderer.selected_element, self)
            if dialog.exec_() == QDialog.Accepted:
                self.renderer.update()
                self.update_element_info()
        else:
            QMessageBox.information(self, "Info", "Please click to select a geological time block")
    
    def update_element_info(self):
        """Update element info"""
        if self.renderer.selected_element:
            element = self.renderer.selected_element
            info = f"Type: {element.element_type}\n"
            info += f"Text: {element.text}\n"
            info += f"Visible: {'Yes' if element.visible else 'No'}\n"
            if hasattr(element, 'start_time'):
                info += f"Time Range: {element.start_time} - {element.end_time} Ma"
            self.element_info_label.setText(info)
    
    def apply_preset(self, min_time, max_time, unit):
        """Apply preset"""
        self.min_time_input.setText(str(min_time))
        self.max_time_input.setText(str(max_time))
        self.time_unit_combo.setCurrentText(f"{unit} (Million Years)" if unit == "Ma" else f"{unit}")
        self.tick_interval_input.setText(str((max_time - min_time) / 10))
        self.update_preview()
    
    def export_pdf(self):
        """Export as PDF"""
        from reportlab.lib.pagesizes import A4
        from reportlab.pdfgen import canvas
        
        filename, _ = QFileDialog.getSaveFileName(self, "Export PDF", "", "PDF Files (*.pdf)")
        if filename:
            try:
                c = canvas.Canvas(filename, pagesize=A4)
                width, height = A4
                
                # 绘制背景
                c.setFillColorRGB(1, 1, 1)
                c.rect(0, 0, width, height, fill=1, stroke=0)
                
                # 绘制主标尺线
                c.setLineWidth(2)
                c.setStrokeColorRGB(0, 0, 0)
                scale_x = 100
                scale_y = height / 2
                scale_length = self.renderer.scale_length
                c.line(scale_x, scale_y, scale_x + scale_length, scale_y)
                
                # 绘制刻度和标签
                min_mya = self.renderer.convert_to_mya(self.renderer.min_time)
                max_mya = self.renderer.convert_to_mya(self.renderer.max_time)
                tick_interval_mya = self.renderer.convert_to_mya(self.renderer.tick_interval)
                
                total_ticks = int((max_mya - min_mya) / tick_interval_mya) + 1
                
                for i in range(total_ticks):
                    tick_value_mya = min_mya + i * tick_interval_mya
                    tick_value_current = self.renderer.convert_from_mya(tick_value_mya)
                    
                    progress = (tick_value_mya - min_mya) / (max_mya - min_mya) if (max_mya - min_mya) > 0 else 0
                    x_pos = scale_x + progress * scale_length
                    
                    is_major = i % 5 == 0
                    tick_len = self.renderer.major_tick_length if is_major else self.renderer.tick_length
                    
                    # 刻度线
                    c.setLineWidth(1 if not is_major else 2)
                    c.setStrokeColorRGB(0, 0, 0)
                    
                    if self.renderer.start_direction == "bottom":
                        c.line(x_pos, scale_y, x_pos, scale_y + tick_len)
                        label_y = scale_y + tick_len + 15
                    else:
                        c.line(x_pos, scale_y, x_pos, scale_y - tick_len)
                        label_y = scale_y - tick_len - 10
                    
                    # 标签
                    if self.renderer.show_labels:
                        c.setFillColorRGB(0, 0, 0)
                        c.setFont("Helvetica", 10 if not is_major else 12)
                        c.drawCentredString(x_pos, label_y, self.renderer.format_time_value(tick_value_current))
                
                # 绘制地质年代
                if self.renderer.show_geological:
                    show_era = self.renderer.geological_scale_type in ["era", "both"]
                    show_period = self.renderer.geological_scale_type in ["period", "both"]
                    
                    base_geo_y = scale_y + self.renderer.scale_thickness + 10 if self.renderer.start_direction == "bottom" else scale_y - 10
                    
                    # 显示代
                    if show_era:
                        era_height = 30
                        era_y = base_geo_y if self.renderer.start_direction == "bottom" else base_geo_y - era_height
                        
                        for period in GEOLOGICAL_ERAS:
                            if period["end"] < min_mya or period["start"] > max_mya:
                                continue
                            
                            start_pos = max(period["start"], min_mya)
                            end_pos = min(period["end"], max_mya)
                            
                            progress_start = (start_pos - min_mya) / (max_mya - min_mya) if (max_mya - min_mya) > 0 else 0
                            progress_end = (end_pos - min_mya) / (max_mya - min_mya) if (max_mya - min_mya) > 0 else 0
                            
                            x1 = scale_x + progress_start * scale_length
                            rect_width = (progress_end - progress_start) * scale_length
                            
                            color = QColor(period["color"])
                            c.setFillColorRGB(color.red()/255, color.green()/255, color.blue()/255)
                            c.setStrokeColorRGB(0, 0, 0)
                            c.rect(x1, era_y, rect_width, era_height, fill=1, stroke=1)
                            
                            c.setFillColorRGB(0, 0, 0)
                            c.setFont("Helvetica", 8)
                            text_x = x1 + rect_width / 2
                            text_y = era_y + era_height / 2
                            c.drawCentredString(text_x, text_y, period["name"])
                    
                    # 显示纪
                    if show_period:
                        period_height = 25
                        if show_era:
                            period_y = base_geo_y + era_height + 5 if self.renderer.start_direction == "bottom" else base_geo_y - era_height - 5 - period_height
                        else:
                            period_y = base_geo_y if self.renderer.start_direction == "bottom" else base_geo_y - period_height
                        
                        for period in GEOLOGICAL_PERIODS:
                            if period["end"] < min_mya or period["start"] > max_mya:
                                continue
                            
                            start_pos = max(period["start"], min_mya)
                            end_pos = min(period["end"], max_mya)
                            
                            progress_start = (start_pos - min_mya) / (max_mya - min_mya) if (max_mya - min_mya) > 0 else 0
                            progress_end = (end_pos - min_mya) / (max_mya - min_mya) if (max_mya - min_mya) > 0 else 0
                            
                            x1 = scale_x + progress_start * scale_length
                            rect_width = (progress_end - progress_start) * scale_length
                            
                            color = QColor(period["color"])
                            c.setFillColorRGB(color.red()/255, color.green()/255, color.blue()/255)
                            c.setStrokeColorRGB(0, 0, 0)
                            c.rect(x1, period_y, rect_width, period_height, fill=1, stroke=1)
                            
                            c.setFillColorRGB(0, 0, 0)
                            c.setFont("Helvetica", 7)
                            text_x = x1 + rect_width / 2
                            text_y = period_y + period_height / 2
                            c.drawCentredString(text_x, text_y, period["name"])
                
                c.save()
                QMessageBox.information(self, "Success", "PDF file exported successfully")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to export PDF: {str(e)}")
    
    def export_svg(self):
        """Export as SVG"""
        filename, _ = QFileDialog.getSaveFileName(self, "Export SVG", "", "SVG Files (*.svg)")
        if filename:
            try:
                import svgwrite
                
                width = int(self.renderer.scale_length) + 200
                height = 600
                
                svg = svgwrite.Drawing(filename, size=(f"{width}px", f"{height}px"))
                
                # 背景
                svg.add(svgwrite.shapes.Rect(insert=(0, 0), size=(width, height), 
                                           fill=self.renderer.background_color.name()))
                
                # 主标尺线
                scale_x = 50
                scale_y = 300
                svg.add(svgwrite.shapes.Line(start=(scale_x, scale_y), 
                                            end=(scale_x + self.renderer.scale_length, scale_y),
                                            stroke=self.renderer.line_color.name(), stroke_width=2))
                
                # 刻度和标签
                min_mya = self.renderer.convert_to_mya(self.renderer.min_time)
                max_mya = self.renderer.convert_to_mya(self.renderer.max_time)
                tick_interval_mya = self.renderer.convert_to_mya(self.renderer.tick_interval)
                
                total_ticks = int((max_mya - min_mya) / tick_interval_mya) + 1
                
                for i in range(total_ticks):
                    tick_value_mya = min_mya + i * tick_interval_mya
                    tick_value_current = self.renderer.convert_from_mya(tick_value_mya)
                    
                    progress = (tick_value_mya - min_mya) / (max_mya - min_mya) if (max_mya - min_mya) > 0 else 0
                    x_pos = scale_x + progress * self.renderer.scale_length
                    
                    is_major = i % 5 == 0
                    tick_len = self.renderer.major_tick_length if is_major else self.renderer.tick_length
                    
                    # 刻度线
                    if self.renderer.start_direction == "bottom":
                        svg.add(svgwrite.shapes.Line(start=(x_pos, scale_y), 
                                                    end=(x_pos, scale_y + tick_len),
                                                    stroke=self.renderer.tick_color.name(), 
                                                    stroke_width=2 if is_major else 1))
                        label_y = scale_y + tick_len + 15
                        label_anchor = "middle"
                    else:
                        svg.add(svgwrite.shapes.Line(start=(x_pos, scale_y), 
                                                    end=(x_pos, scale_y - tick_len),
                                                    stroke=self.renderer.tick_color.name(), 
                                                    stroke_width=2 if is_major else 1))
                        label_y = scale_y - tick_len - 10
                        label_anchor = "middle"
                    
                    # 标签
                    if self.renderer.show_labels:
                        svg.add(svgwrite.text.Text(self.renderer.format_time_value(tick_value_current),
                                                  insert=(x_pos, label_y),
                                                  fill=self.renderer.label_color.name(),
                                                  font_size="10px" if not is_major else "12px",
                                                  text_anchor=label_anchor))
                
                # 地质年代
                if self.renderer.show_geological:
                    show_era = self.renderer.geological_scale_type in ["era", "both"]
                    show_period = self.renderer.geological_scale_type in ["period", "both"]
                    
                    base_geo_y = scale_y + self.renderer.scale_thickness + 10 if self.renderer.start_direction == "bottom" else scale_y - 10
                    
                    # 显示代
                    if show_era:
                        era_height = 30
                        era_y = base_geo_y if self.renderer.start_direction == "bottom" else base_geo_y - era_height
                        
                        for period in GEOLOGICAL_ERAS:
                            if period["end"] < min_mya or period["start"] > max_mya:
                                continue
                            
                            start_pos = max(period["start"], min_mya)
                            end_pos = min(period["end"], max_mya)
                            
                            progress_start = (start_pos - min_mya) / (max_mya - min_mya) if (max_mya - min_mya) > 0 else 0
                            progress_end = (end_pos - min_mya) / (max_mya - min_mya) if (max_mya - min_mya) > 0 else 0
                            
                            x1 = scale_x + progress_start * self.renderer.scale_length
                            rect_width = (progress_end - progress_start) * self.renderer.scale_length
                            
                            # 年代块
                            svg.add(svgwrite.shapes.Rect(insert=(x1, era_y), 
                                                        size=(rect_width, era_height),
                                                        fill=period["color"], stroke="black"))
                            
                            # 年代名称
                            text_x = x1 + rect_width / 2
                            text_y = era_y + era_height / 2
                            svg.add(svgwrite.text.Text(period["name"], 
                                                      insert=(text_x, text_y),
                                                      fill="black", font_size="8px",
                                                      text_anchor="middle", dominant_baseline="middle"))
                    
                    # 显示纪
                    if show_period:
                        period_height = 25
                        if show_era:
                            period_y = base_geo_y + era_height + 5 if self.renderer.start_direction == "bottom" else base_geo_y - era_height - 5 - period_height
                        else:
                            period_y = base_geo_y if self.renderer.start_direction == "bottom" else base_geo_y - period_height
                        
                        for period in GEOLOGICAL_PERIODS:
                            if period["end"] < min_mya or period["start"] > max_mya:
                                continue
                            
                            start_pos = max(period["start"], min_mya)
                            end_pos = min(period["end"], max_mya)
                            
                            progress_start = (start_pos - min_mya) / (max_mya - min_mya) if (max_mya - min_mya) > 0 else 0
                            progress_end = (end_pos - min_mya) / (max_mya - min_mya) if (max_mya - min_mya) > 0 else 0
                            
                            x1 = scale_x + progress_start * self.renderer.scale_length
                            rect_width = (progress_end - progress_start) * self.renderer.scale_length
                            
                            # 年代块
                            svg.add(svgwrite.shapes.Rect(insert=(x1, period_y), 
                                                        size=(rect_width, period_height),
                                                        fill=period["color"], stroke="black"))
                            
                            # 年代名称
                            text_x = x1 + rect_width / 2
                            text_y = period_y + period_height / 2
                            svg.add(svgwrite.text.Text(period["name"], 
                                                      insert=(text_x, text_y),
                                                      fill="black", font_size="7px",
                                                      text_anchor="middle", dominant_baseline="middle"))
                
                svg.save()
                QMessageBox.information(self, "Success", "SVG file exported successfully")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to export SVG: {str(e)}")
    
    def export_png(self):
        """Export as PNG"""
        from PyQt5.QtGui import QPixmap
        filename, _ = QFileDialog.getSaveFileName(self, "Export PNG", "", "PNG Files (*.png)")
        if filename:
            try:
                # Render as image
                pixmap = QPixmap(self.renderer.size())
                self.renderer.render(pixmap)
                pixmap.save(filename, "PNG")
                QMessageBox.information(self, "Success", "PNG file exported successfully")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to export PNG: {str(e)}")

def main():
    app = QApplication(sys.argv)
    
    # 设置应用样式
    app.setStyle('Fusion')
    
    window = TimeScaleGenerator()
    window.show()
    
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
