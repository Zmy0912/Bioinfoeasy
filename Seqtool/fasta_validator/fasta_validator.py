#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FASTA Format Validator and Fixer
Supports FASTA format detection, error reporting, and automatic fixing
Includes ATCG strict validation mode and error position highlighting
"""

import tkinter as tk
from tkinter import ttk, filedialog, scrolledtext, messagebox, font
import os
import re


class FASTAValidator:
    """FASTA format validator and fixer class"""

    @staticmethod
    def is_fasta_file(content):
        """Check if file content follows FASTA format"""
        lines = content.strip().split('\n')
        if not lines:
            return False
        return lines[0].startswith('>')

    @staticmethod
    def validate_fasta(content, strict_atcg=False, mark_errors=False):
        """Validate FASTA format, return error list and error position info

        Args:
            content: FASTA file content
            strict_atcg: Whether to only allow ATCG characters (strict mode)
            mark_errors: Whether to mark error positions
        """
        errors = []
        error_positions = []
        lines = content.split('\n')

        if not lines or all(line.strip() == '' for line in lines):
            errors.append("File is empty")
            return errors, error_positions

        # Check if first line starts with '>'
        first_line = lines[0].strip()
        if not first_line.startswith('>'):
            errors.append("Line 1: First line must start with '>' (header line)")

        # Check sequence lines
        in_sequence = False
        seq_count = 0
        current_seq_id = ""
        current_seq_lines = []

        # Valid sequence characters
        if strict_atcg:
            valid_chars = set('ACGT-')
            valid_chars_name = "ATCG"
        else:
            valid_chars = set('ACGTURYSWKMBDHVN-')
            valid_chars_name = "ACGTURYSWKMBDHVN"

        for i, line in enumerate(lines, 1):
            line_stripped = line.strip()

            if line_stripped.startswith('>'):
                # Save previous sequence's error positions
                if current_seq_id and mark_errors and current_seq_lines:
                    error_positions.append({
                        'seq_id': current_seq_id,
                        'lines': current_seq_lines,
                        'line_numbers': list(range(i - len(current_seq_lines), i))
                    })

                # New sequence
                in_sequence = True
                seq_count += 1
                current_seq_id = line_stripped[1:].strip()
                current_seq_lines = []
                if not current_seq_id:
                    errors.append(f"Line {i}: Header line missing sequence ID")
            elif line_stripped:
                # Sequence line
                if not in_sequence:
                    errors.append(f"Line {i}: Sequence data found without header line")

                # Check for invalid characters
                invalid_positions = []
                line_display = []
                error_markers = []

                for pos, char in enumerate(line_stripped):
                    char_upper = char.upper()
                    if char_upper in valid_chars:
                        line_display.append(char)
                        error_markers.append(' ')
                    elif char in ' \t':
                        line_display.append(' ')
                        error_markers.append(' ')
                    else:
                        line_display.append(char)
                        error_markers.append('^')  # Mark invalid characters with ^

                invalid_chars = [c for c, m in zip(line_display, error_markers) if m == '^']
                if invalid_chars:
                    unique_invalid = sorted(set(invalid_chars))
                    errors.append(f"Line {i}: Contains invalid characters '{''.join(unique_invalid)}' (expected {valid_chars_name})")

                    # Save error position info
                    if mark_errors:
                        current_seq_lines.append({
                            'line_num': i,
                            'original': line_stripped,
                            'display': ''.join(line_display),
                            'markers': ''.join(error_markers),
                            'invalid_chars': unique_invalid
                        })

        # Save last sequence's error positions
        if current_seq_id and mark_errors and current_seq_lines:
            error_positions.append({
                'seq_id': current_seq_id,
                'lines': current_seq_lines,
                'line_numbers': list(range(len(lines) - len(current_seq_lines) + 1, len(lines) + 1))
            })

        if seq_count == 0:
            errors.append("No sequence headers found in file (lines starting with '>')")

        return errors, error_positions

    @staticmethod
    def fix_fasta(content, strict_atcg=False):
        """Fix FASTA format errors

        Args:
            content: FASTA file content
            strict_atcg: Whether to strictly filter to ATCG characters
        """
        lines = content.split('\n')
        fixed_lines = []
        seq_count = 0
        in_sequence = False

        # Valid sequence characters
        if strict_atcg:
            valid_chars = set('ACGT')
        else:
            valid_chars = set('ACGTURYSWKMBDHVN')

        # Process each line
        for i, line in enumerate(lines, 1):
            line_stripped = line.strip()

            if line_stripped.startswith('>'):
                # Header line
                in_sequence = True
                seq_count += 1
                header = line_stripped[1:].strip()

                # If header is empty, generate default
                if not header:
                    header = f"sequence_{seq_count}"

                fixed_lines.append(f">{header}")
            elif line_stripped:
                # Sequence line, remove spaces and numbers
                cleaned_seq = ""
                for char in line_stripped.upper():
                    if char in valid_chars:
                        cleaned_seq += char
                    elif char == '-':
                        cleaned_seq += char  # Keep gap symbol

                if cleaned_seq:
                    if not in_sequence and i > 1:
                        # Sequence without header, add header
                        seq_count += 1
                        fixed_lines.append(f">sequence_{seq_count}")
                        in_sequence = True
                    fixed_lines.append(cleaned_seq)
            # Skip empty lines

        # If no sequences, return empty
        if not fixed_lines:
            return ""

        return '\n'.join(fixed_lines)


class FASTAValidatorGUI:
    """FASTA validator graphical user interface"""

    def __init__(self, root):
        self.root = root
        self.root.title("FASTA Format Validator and Fixer")
        self.root.geometry("1200x900")
        self.root.minsize(1000, 750)

        self.current_file = None
        self.original_content = ""
        self.strict_atcg = tk.BooleanVar(value=False)
        self.show_error_markers = tk.BooleanVar(value=True)
        self.current_error_positions = []

        self.setup_ui()

    def setup_ui(self):
        """Setup user interface"""
        # Top toolbar
        toolbar = ttk.Frame(self.root, padding="10")
        toolbar.pack(fill=tk.X)

        ttk.Button(toolbar, text="Select File", command=self.load_file).pack(side=tk.LEFT, padx=5)
        ttk.Button(toolbar, text="Validate Format", command=self.validate_fasta).pack(side=tk.LEFT, padx=5)
        ttk.Button(toolbar, text="Fix File", command=self.fix_fasta).pack(side=tk.LEFT, padx=5)
        ttk.Button(toolbar, text="Save Fixed", command=self.save_fixed_file).pack(side=tk.LEFT, padx=5)
        ttk.Button(toolbar, text="Clear", command=self.clear_all).pack(side=tk.LEFT, padx=5)

        # Separator
        ttk.Separator(self.root, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=10, pady=5)

        # Settings area
        settings_frame = ttk.LabelFrame(self.root, text="Validation Options", padding="10")
        settings_frame.pack(fill=tk.X, padx=10, pady=5)

        # ATCG strict mode checkbox
        atcg_check = ttk.Checkbutton(
            settings_frame,
            text="Strict ATCG Mode (A, C, G, T only)",
            variable=self.strict_atcg,
            command=self.update_mode_label
        )
        atcg_check.pack(side=tk.LEFT, padx=10)

        # Show error markers checkbox
        marker_check = ttk.Checkbutton(
            settings_frame,
            text="Highlight Error Positions (mark invalid chars with ^)",
            variable=self.show_error_markers
        )
        marker_check.pack(side=tk.LEFT, padx=10)

        # Mode label
        self.mode_label = ttk.Label(settings_frame, text="Current Mode: Standard (ACGTURYSWKMBDHVN)", font=('Arial', 9))
        self.mode_label.pack(side=tk.LEFT, padx=20)

        # Separator
        ttk.Separator(self.root, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=10, pady=5)

        # File info
        info_frame = ttk.Frame(self.root, padding="10")
        info_frame.pack(fill=tk.X)

        self.file_label = ttk.Label(info_frame, text="Current File: Not selected", font=('Arial', 10))
        self.file_label.pack(side=tk.LEFT)

        # Create main container
        main_container = ttk.Frame(self.root)
        main_container.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        # Main content area - use PanedWindow for split
        paned = ttk.PanedWindow(main_container, orient=tk.VERTICAL)
        paned.pack(fill=tk.BOTH, expand=True)

        # Original content area
        original_frame = ttk.LabelFrame(paned, text="Original File Content", padding="5")
        paned.add(original_frame, weight=1)

        self.original_text = scrolledtext.ScrolledText(
            original_frame,
            wrap=tk.NONE,
            font=('Consolas', 10),
            height=8
        )
        self.original_text.pack(fill=tk.BOTH, expand=True)

        # Error marker area
        error_marker_frame = ttk.LabelFrame(paned, text="Error Position Markers (^ indicates invalid chars)", padding="5")
        paned.add(error_marker_frame, weight=1)

        self.error_marker_text = scrolledtext.ScrolledText(
            error_marker_frame,
            wrap=tk.NONE,
            font=('Consolas', 10),
            height=8
        )
        self.error_marker_text.pack(fill=tk.BOTH, expand=True)

        # Configure marker area tags
        self.error_marker_text.tag_config('seq_header', foreground='blue', font=('Consolas', 10, 'bold'))
        self.error_marker_text.tag_config('seq_line', foreground='black')
        self.error_marker_text.tag_config('error_marker', foreground='red', font=('Consolas', 10, 'bold'))
        self.error_marker_text.tag_config('info_text', foreground='gray', font=('Consolas', 9))

        # Fixed content area
        fixed_frame = ttk.LabelFrame(paned, text="Fixed Content", padding="5")
        paned.add(fixed_frame, weight=1)

        self.fixed_text = scrolledtext.ScrolledText(
            fixed_frame,
            wrap=tk.NONE,
            font=('Consolas', 10),
            height=8
        )
        self.fixed_text.pack(fill=tk.BOTH, expand=True)

        # Error information area
        error_frame = ttk.LabelFrame(paned, text="Validation Results / Error Information", padding="5")
        paned.add(error_frame, weight=1)

        self.error_text = scrolledtext.ScrolledText(
            error_frame,
            wrap=tk.WORD,
            font=('Arial', 9),
            height=8
        )
        self.error_text.pack(fill=tk.BOTH, expand=True)

        # Configure tags
        self.error_text.tag_config('header', foreground='blue', font=('Arial', 10, 'bold'))
        self.error_text.tag_config('error', foreground='red')
        self.error_text.tag_config('warning', foreground='orange')
        self.error_text.tag_config('success', foreground='green')
        self.error_text.tag_config('info', foreground='black')

        # Status bar
        self.status_var = tk.StringVar(value="Ready")
        statusbar = ttk.Label(self.root, textvariable=self.status_var, relief=tk.SUNKEN)
        statusbar.pack(fill=tk.X, padx=10, pady=5)

    def update_mode_label(self):
        """Update mode label"""
        if self.strict_atcg.get():
            self.mode_label.config(text="Current Mode: Strict ATCG (A, C, G, T only)")
        else:
            self.mode_label.config(text="Current Mode: Standard (ACGTURYSWKMBDHVN)")

    def load_file(self):
        """Load FASTA file"""
        filepath = filedialog.askopenfilename(
            title="Select FASTA File",
            filetypes=[
                ("FASTA Files", "*.fasta *.fa *.fas *.fna *.ffn *.faa *.frn"),
                ("Text Files", "*.txt"),
                ("All Files", "*.*")
            ]
        )

        if filepath:
            try:
                with open(filepath, 'r', encoding='utf-8') as f:
                    content = f.read()

                # Try other encodings
                if not content:
                    for encoding in ['gbk', 'gb2312', 'latin-1']:
                        try:
                            with open(filepath, 'r', encoding=encoding) as f:
                                content = f.read()
                            break
                        except:
                            continue

                self.current_file = filepath
                self.original_content = content
                self.original_text.delete(1.0, tk.END)
                self.original_text.insert(1.0, content)
                self.fixed_text.delete(1.0, tk.END)
                self.error_marker_text.delete(1.0, tk.END)
                self.file_label.config(text=f"Current File: {os.path.basename(filepath)}")
                self.status_var.set(f"Loaded: {os.path.basename(filepath)}")

                # Clear error info
                self.error_text.delete(1.0, tk.END)

                # Auto-detect FASTA format
                if FASTAValidator.is_fasta_file(content):
                    self.status_var.set(f"Loaded: {os.path.basename(filepath)} (FASTA format)")
                    self.append_error("[OK] File follows FASTA format\n", "info")
                else:
                    self.status_var.set(f"Loaded: {os.path.basename(filepath)} (may not be FASTA)")
                    self.append_error("[!] File may not follow FASTA format\n", "warning")

            except Exception as e:
                messagebox.showerror("Error", f"Cannot read file: {str(e)}")
                self.status_var.set("File load failed")

    def validate_fasta(self):
        """Validate FASTA format"""
        content = self.original_text.get(1.0, tk.END).strip()

        if not content:
            messagebox.showwarning("Warning", "Please select a file first")
            return

        self.error_text.delete(1.0, tk.END)
        self.error_marker_text.delete(1.0, tk.END)
        self.error_text.see(1.0)

        mode_str = " (Strict ATCG Mode)" if self.strict_atcg.get() else ""
        self.append_error(f"=== FASTA Format Validation{mode_str} ===\n\n", "header")

        # Detailed validation with error position marking
        errors, error_positions = FASTAValidator.validate_fasta(
            content,
            strict_atcg=self.strict_atcg.get(),
            mark_errors=self.show_error_markers.get()
        )

        self.current_error_positions = error_positions

        if not errors:
            self.append_error("[OK] File format is correct, no errors found!\n", "success")
            self.status_var.set("Validation Passed: No errors")
        else:
            self.append_error(f"[X] Found {len(errors)} error(s):\n\n", "error")
            for i, error in enumerate(errors, 1):
                self.append_error(f"{i}. {error}\n", "error")
            self.status_var.set(f"Validation Complete: {len(errors)} error(s) found")

        # Display error position markers
        if self.show_error_markers.get() and error_positions:
            self.display_error_markers(error_positions)

        # Auto-scroll to error area top
        self.root.after(100, lambda: self.error_text.see(1.0))

    def display_error_markers(self, error_positions):
        """Display error position markers in text box"""
        self.error_marker_text.insert(tk.END, "=== Error Position Markers (^ indicates invalid characters) ===\n\n", "seq_header")
        self.error_marker_text.insert(tk.END, "Legend:\n", "info_text")
        self.error_marker_text.insert(tk.END, "  ^ marks the position of invalid characters\n", "info_text")
        self.error_marker_text.insert(tk.END, "  Blue: Sequence header line\n", "info_text")
        self.error_marker_text.insert(tk.END, "  Black: Original sequence line\n", "info_text")
        self.error_marker_text.insert(tk.END, "  Red: Error marker line\n\n", "info_text")

        for seq_info in error_positions:
            # Display sequence header
            self.error_marker_text.insert(tk.END, f">{seq_info['seq_id']}\n", "seq_header")

            # Display each line and its error markers
            for line_info in seq_info['lines']:
                line_num = line_info['line_num']
                original = line_info['original']
                display = line_info['display']
                markers = line_info['markers']
                invalid_chars = line_info['invalid_chars']

                # Display line number
                self.error_marker_text.insert(tk.END, f"Line {line_num}: ", "info_text")

                # Display original line
                self.error_marker_text.insert(tk.END, f"{original}\n", "seq_line")

                # Display error markers
                if '^' in markers:
                    # Calculate indentation (align character positions)
                    indent = " " * (len(f"Line {line_num}: "))
                    self.error_marker_text.insert(tk.END, f"{indent}{markers}\n", "error_marker")
                    self.error_marker_text.insert(
                        tk.END,
                        f"{indent}^ Invalid characters: {''.join(sorted(invalid_chars))}\n",
                        "info_text"
                    )

            self.error_marker_text.insert(tk.END, "\n")

        # Scroll to marker area top
        self.root.after(100, lambda: self.error_marker_text.see(1.0))

    def fix_fasta(self):
        """Fix FASTA format"""
        content = self.original_text.get(1.0, tk.END).strip()

        if not content:
            messagebox.showwarning("Warning", "Please select a file first")
            return

        try:
            fixed_content = FASTAValidator.fix_fasta(content, strict_atcg=self.strict_atcg.get())

            if fixed_content:
                self.fixed_text.delete(1.0, tk.END)
                self.fixed_text.insert(1.0, fixed_content)

                # Count differences before and after fixing
                original_lines = len([l for l in content.split('\n') if l.strip()])
                fixed_lines = len([l for l in fixed_content.split('\n') if l.strip()])

                self.error_text.delete(1.0, tk.END)
                mode_str = " (Strict ATCG Mode)" if self.strict_atcg.get() else ""
                self.append_error(f"=== FASTA File Fixed{mode_str} ===\n\n", "header")
                self.append_error(f"Original lines: {original_lines}\n", "info")
                self.append_error(f"Fixed lines: {fixed_lines}\n\n", "info")
                self.append_error("[OK] Fixed content shown in the text box below\n", "success")
                self.append_error("Tip: Please review the result and click 'Save Fixed' if satisfied\n\n", "info")

                self.status_var.set("File fixed, please review")
            else:
                self.error_text.delete(1.0, tk.END)
                self.append_error("[X] Cannot fix this file: empty or severely corrupted format\n", "error")
                self.status_var.set("Fix failed")

        except Exception as e:
            messagebox.showerror("Error", f"Error fixing file: {str(e)}")
            self.status_var.set("Fix failed")

    def save_fixed_file(self):
        """Save fixed file"""
        fixed_content = self.fixed_text.get(1.0, tk.END).strip()

        if not fixed_content:
            messagebox.showwarning("Warning", "No content to save")
            return

        # Default filename
        if self.current_file:
            default_name = os.path.splitext(self.current_file)[0] + "_fixed.fasta"
        else:
            default_name = "fixed_sequence.fasta"

        filepath = filedialog.asksaveasfilename(
            title="Save Fixed FASTA File",
            defaultextension=".fasta",
            initialfile=default_name,
            filetypes=[
                ("FASTA Files", "*.fasta"),
                ("Text Files", "*.txt"),
                ("All Files", "*.*")
            ]
        )

        if filepath:
            try:
                with open(filepath, 'w', encoding='utf-8') as f:
                    f.write(fixed_content)
                messagebox.showinfo("Success", f"File saved:\n{filepath}")
                self.status_var.set(f"Saved: {os.path.basename(filepath)}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save file: {str(e)}")

    def clear_all(self):
        """Clear all content"""
        self.original_text.delete(1.0, tk.END)
        self.fixed_text.delete(1.0, tk.END)
        self.error_text.delete(1.0, tk.END)
        self.error_marker_text.delete(1.0, tk.END)
        self.current_file = None
        self.original_content = ""
        self.current_error_positions = []
        self.file_label.config(text="Current File: Not selected")
        self.status_var.set("Cleared")

    def append_error(self, message, msg_type="normal"):
        """Add error message to text box"""
        # Get cursor position first
        start_index = self.error_text.index("insert")
        self.error_text.insert(tk.END, message)

        # Apply tag
        self.error_text.tag_add(msg_type, start_index, tk.END)


def main():
    """Main function"""
    root = tk.Tk()
    app = FASTAValidatorGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
