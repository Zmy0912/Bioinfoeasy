#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
FASTA Sequence Joining Tool (v1.1.0)
- Single-file mode  : concatenate all sequences from one FASTA file into one sequence
- Batch-folder mode : process every FASTA file inside a folder; each file is turned
                      into one long sequence and saved under the SAME filename
"""

import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
from pathlib import Path

FASTA_EXTS = ('.fasta', '.fa', '.fna', '.ffn', '.faa', '.frn')


class FastaJoinerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("FASTA Sequence Joiner")
        self.root.geometry("680x600")
        self.root.resizable(True, True)

        # Stored paths
        self.input_file = ""    # single-file mode input path
        self.input_dir = ""     # batch mode input folder
        self.output_file = ""   # single-file mode output file path
        self.output_dir = ""    # batch mode output folder

        # Working mode: "single" or "batch"
        self.mode = tk.StringVar(value="single")

        self.create_widgets()

    def create_widgets(self):
        """Create interface components"""
        # Title
        title_label = tk.Label(
            self.root,
            text="FASTA Sequence Joiner",
            font=("Arial", 16, "bold"),
            pady=10
        )
        title_label.pack()

        # Description text
        info_frame = tk.Frame(self.root, padx=20)
        info_frame.pack(fill=tk.X, pady=5)

        info_text = """
This tool concatenates all sequences from FASTA file(s) into one.
Single-file mode: merge fragments of one file.  Batch mode: merge every
FASTA file inside a folder, each result keeps the original file name.
        """
        info_label = tk.Label(info_frame, text=info_text, justify=tk.LEFT, fg="gray")
        info_label.pack()

        # Mode selection
        mode_frame = tk.LabelFrame(
            self.root,
            text="Working Mode",
            font=("Arial", 10, "bold"),
            padx=10, pady=5
        )
        mode_frame.pack(fill=tk.X, padx=20, pady=5)

        single_rb = tk.Radiobutton(
            mode_frame,
            text="Single FASTA File",
            variable=self.mode,
            value="single",
            command=self.on_mode_change,
            font=("Arial", 9)
        )
        single_rb.pack(side=tk.LEFT, padx=(0, 20))

        batch_rb = tk.Radiobutton(
            mode_frame,
            text="Batch Folder (all FASTA files)",
            variable=self.mode,
            value="batch",
            command=self.on_mode_change,
            font=("Arial", 9)
        )
        batch_rb.pack(side=tk.LEFT)

        self.mode_hint = tk.Label(mode_frame, text="", fg="gray", font=("Arial", 8))
        self.mode_hint.pack(fill=tk.X, pady=(5, 0))

        # Step 1: input
        self.step1_frame = tk.LabelFrame(
            self.root,
            text="Step 1: Select FASTA Input File",
            font=("Arial", 10, "bold"),
            padx=10, pady=10
        )
        self.step1_frame.pack(fill=tk.X, padx=20, pady=5)

        input_path_frame = tk.Frame(self.step1_frame)
        input_path_frame.pack(fill=tk.X, pady=5)

        self.input_path_label = tk.Label(
            input_path_frame,
            text="Input file path:",
            font=("Arial", 9)
        )
        self.input_path_label.pack(side=tk.LEFT, padx=(0, 5))

        self.input_path_entry = tk.Entry(
            input_path_frame,
            font=("Arial", 9),
            width=50
        )
        self.input_path_entry.pack(side=tk.LEFT, fill=tk.X, expand=True)

        browse_input_btn = tk.Button(
            input_path_frame,
            text="Browse",
            command=self.browse_input_file,
            width=8
        )
        browse_input_btn.pack(side=tk.LEFT, padx=(5, 0))

        self.input_status_label = tk.Label(self.step1_frame, text="", fg="gray", font=("Arial", 8))
        self.input_status_label.pack(fill=tk.X, pady=(5, 0))

        self.input_path_entry.bind("<FocusOut>", self.validate_input_path)
        self.input_path_entry.bind("<KeyRelease>", lambda e: self.root.after(500, self.validate_input_path))

        # Step 2: output options
        self.step2_frame = tk.LabelFrame(
            self.root,
            text="Step 2: Set Output Options",
            font=("Arial", 10, "bold"),
            padx=10, pady=10
        )
        self.step2_frame.pack(fill=tk.X, padx=20, pady=5)

        # Sequence name input
        name_frame = tk.Frame(self.step2_frame)
        name_frame.pack(fill=tk.X, pady=5)

        self.name_label = tk.Label(name_frame, text="New sequence name:", font=("Arial", 9))
        self.name_label.pack(side=tk.LEFT, padx=(0, 5))

        self.seq_name_entry = tk.Entry(
            name_frame,
            font=("Arial", 9),
            width=40
        )
        self.seq_name_entry.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.seq_name_entry.insert(0, "joined_sequence")  # Default value
        self.seq_name_entry.config(state="readonly")  # Default readonly

        # Add custom name button (only shown in single-file mode)
        self.custom_name_btn = tk.Button(
            name_frame,
            text="Custom",
            command=self.enable_custom_name,
            font=("Arial", 8),
            padx=5
        )
        self.custom_name_btn.pack(side=tk.LEFT, padx=(5, 0))

        # Path selection
        path_frame = tk.Frame(self.step2_frame)
        path_frame.pack(fill=tk.X, pady=5)

        self.output_path_label = tk.Label(
            path_frame,
            text="Output file path:",
            font=("Arial", 9)
        )
        self.output_path_label.pack(side=tk.LEFT, padx=(0, 5))

        self.output_path_entry = tk.Entry(
            path_frame,
            font=("Arial", 9),
            width=50
        )
        self.output_path_entry.pack(side=tk.LEFT, fill=tk.X, expand=True)

        browse_output_btn = tk.Button(
            path_frame,
            text="Browse",
            command=self.browse_output_file,
            width=8
        )
        browse_output_btn.pack(side=tk.LEFT, padx=(5, 0))

        self.output_status_label = tk.Label(self.step2_frame, text="", fg="gray", font=("Arial", 8))
        self.output_status_label.pack(fill=tk.X, pady=(5, 0))

        self.output_path_entry.bind("<FocusOut>", self.validate_output_path)
        self.output_path_entry.bind("<KeyRelease>", lambda e: self.root.after(500, self.validate_output_path))

        # Step 3: Execute
        step3_frame = tk.LabelFrame(
            self.root,
            text="Step 3: Execute Joining",
            font=("Arial", 10, "bold"),
            padx=10, pady=10
        )
        step3_frame.pack(fill=tk.X, padx=20, pady=5)

        run_btn = tk.Button(
            step3_frame,
            text="Start Joining Sequences",
            command=self.join_sequences,
            bg="#4CAF50",
            fg="white",
            font=("Arial", 11, "bold"),
            padx=20,
            pady=5
        )
        run_btn.pack()

        # Result display area
        result_frame = tk.LabelFrame(
            self.root,
            text="Processing Result",
            font=("Arial", 10, "bold"),
            padx=10, pady=10
        )
        result_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=10)

        self.result_text = scrolledtext.ScrolledText(
            result_frame,
            height=8,
            wrap=tk.WORD
        )
        self.result_text.pack(fill=tk.BOTH, expand=True)

        # Help button
        help_btn = tk.Button(
            self.root,
            text="Help",
            command=self.show_help
        )
        help_btn.pack(pady=5)

    # ------------------------------------------------------------------ #
    # Mode switching
    # ------------------------------------------------------------------ #
    def on_mode_change(self):
        """Switch UI texts/logic between single-file and batch-folder modes"""
        mode = self.mode.get()
        if mode == "single":
            self.mode_hint.config(
                text="Join all sequences from ONE FASTA file into a single sequence."
            )
            self.step1_frame.config(text="Step 1: Select FASTA Input File")
            self.input_path_label.config(text="Input file path:")
            self.step2_frame.config(text="Step 2: Set Output Options")
            self.output_path_label.config(text="Output file path:")
            self.name_label.config(text="New sequence name:")

            self.seq_name_entry.config(state="normal")
            self.seq_name_entry.delete(0, tk.END)
            self.seq_name_entry.insert(0, "joined_sequence")
            self.seq_name_entry.config(state="readonly")
            self.custom_name_btn.pack(side=tk.LEFT, padx=(5, 0))
        else:
            self.mode_hint.config(
                text="Join sequences inside EVERY FASTA file of the folder. Each output\n"
                     "keeps the original file name and contains ONE long sequence."
            )
            self.step1_frame.config(text="Step 1: Select Input Folder")
            self.input_path_label.config(text="Input folder path:")
            self.step2_frame.config(text="Step 2: Set Output Folder")
            self.output_path_label.config(text="Output folder path:")
            self.name_label.config(text="Name suffix (optional):")

            self.seq_name_entry.config(state="normal")
            self.seq_name_entry.delete(0, tk.END)
            self.seq_name_entry.insert(0, "")
            self.custom_name_btn.pack_forget()

        # Reset previously validated paths when mode changes
        self.input_file = ""
        self.input_dir = ""
        self.output_file = ""
        self.output_dir = ""
        self.input_status_label.config(text="", fg="gray")
        self.output_status_label.config(text="", fg="gray")

    def enable_custom_name(self):
        """Enable custom sequence name (single-file mode)"""
        current_state = self.seq_name_entry.cget("state")
        if current_state == "readonly":
            self.seq_name_entry.config(state="normal")
            self.seq_name_entry.focus()

    # ------------------------------------------------------------------ #
    # File / folder browsing
    # ------------------------------------------------------------------ #
    def browse_input_file(self, event=None):
        """Browse input file (single) or folder (batch)"""
        if self.mode.get() == "single":
            file_path = filedialog.askopenfilename(
                title="Select FASTA File",
                filetypes=[("FASTA files", "*.fasta *.fa *.fna *.ffn *.faa *.frn"),
                           ("All files", "*.*")]
            )
            if file_path:
                self.input_path_entry.delete(0, tk.END)
                self.input_path_entry.insert(0, file_path)
                self.validate_input_path()
                self.log(f"✓ Input file selected: {file_path}")
        else:
            folder = filedialog.askdirectory(
                title="Select Folder Containing FASTA Files"
            )
            if folder:
                self.input_path_entry.delete(0, tk.END)
                self.input_path_entry.insert(0, folder)
                self.validate_input_path()
                self.log(f"✓ Input folder selected: {folder}")

    def browse_output_file(self, event=None):
        """Browse output file (single) or output folder (batch)"""
        if self.mode.get() == "single":
            file_path = filedialog.asksaveasfilename(
                title="Save Output File",
                defaultextension=".fasta",
                filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")],
                initialfile="joined_sequence.fasta"
            )
            if file_path:
                self.output_path_entry.delete(0, tk.END)
                self.output_path_entry.insert(0, file_path)
                self.validate_output_path()
                self.log(f"✓ Output file set: {file_path}")
        else:
            folder = filedialog.askdirectory(
                title="Select Output Folder (will be created if missing)"
            )
            if folder:
                self.output_path_entry.delete(0, tk.END)
                self.output_path_entry.insert(0, folder)
                self.validate_output_path()
                self.log(f"✓ Output folder set: {folder}")

    # ------------------------------------------------------------------ #
    # Path validation
    # ------------------------------------------------------------------ #
    def validate_input_path(self, event=None):
        """Validate input path (file in single mode / folder in batch mode)"""
        path = self.input_path_entry.get().strip()
        if not path:
            self.input_file = ""
            self.input_dir = ""
            self.input_status_label.config(text="", fg="gray")
            return

        if self.mode.get() == "single":
            if not os.path.exists(path):
                self.input_status_label.config(text="❌ File does not exist", fg="red")
                self.input_file = ""
            elif os.path.isdir(path):
                self.input_status_label.config(text="❌ Path is a directory, not a file", fg="red")
                self.input_file = ""
            elif path.lower().endswith(FASTA_EXTS):
                self.input_file = path
                file_size = os.path.getsize(path)
                size_str = (f"{file_size / 1024:.1f} KB" if file_size < 1024 * 1024
                            else f"{file_size / (1024 * 1024):.1f} MB")
                self.input_status_label.config(text=f"✓ Valid file ({size_str})", fg="green")
            else:
                self.input_status_label.config(text="⚠ Not a FASTA format file", fg="orange")
                self.input_file = ""
        else:
            if not os.path.exists(path):
                self.input_status_label.config(text="❌ Folder does not exist", fg="red")
                self.input_dir = ""
            elif not os.path.isdir(path):
                self.input_status_label.config(text="❌ Path is a file, not a folder", fg="red")
                self.input_dir = ""
            else:
                count = self.count_fasta_files(path)
                if count == 0:
                    self.input_status_label.config(
                        text="⚠ No FASTA files (*.fasta *.fa *.fna *.ffn *.faa *.frn) found",
                        fg="orange")
                    self.input_dir = ""
                else:
                    self.input_dir = path
                    self.input_status_label.config(
                        text=f"✓ Valid folder ({count} FASTA file(s) found)", fg="green")

    def validate_output_path(self, event=None):
        """Validate output path (file in single mode / folder in batch mode)"""
        path = self.output_path_entry.get().strip()
        if not path:
            self.output_file = ""
            self.output_dir = ""
            self.output_status_label.config(text="", fg="gray")
            return

        try:
            Path(path)
        except Exception:
            self.output_status_label.config(text="❌ Invalid path format", fg="red")
            self.output_file = ""
            self.output_dir = ""
            return

        if self.mode.get() == "single":
            output_dir = os.path.dirname(path)
            if not output_dir:
                output_dir = "."
            if not os.path.exists(output_dir):
                self.output_status_label.config(text="❌ Target directory does not exist", fg="red")
                self.output_file = ""
            elif not os.access(output_dir, os.W_OK):
                self.output_status_label.config(text="❌ Target directory is not writable", fg="red")
                self.output_file = ""
            elif not path.lower().endswith(('.fasta', '.fa')):
                self.output_status_label.config(text="⚠ Suggest using .fasta or .fa extension", fg="orange")
                self.output_file = path
            else:
                if os.path.exists(path):
                    self.output_status_label.config(text="⚠ File will be overwritten", fg="orange")
                else:
                    self.output_status_label.config(text="✓ Path is valid", fg="green")
                self.output_file = path
        else:
            # Batch mode: treat path as output folder
            if os.path.exists(path) and not os.path.isdir(path):
                self.output_status_label.config(text="❌ Path exists and is not a folder", fg="red")
                self.output_dir = ""
            elif not os.access(os.path.dirname(path) if os.path.dirname(path) else ".", os.W_OK) \
                    and not os.path.exists(path):
                self.output_status_label.config(text="❌ Parent directory is not writable", fg="red")
                self.output_dir = ""
            else:
                self.output_dir = path
                if os.path.isdir(path):
                    same = (self.input_dir and
                            os.path.realpath(path).lower() == os.path.realpath(self.input_dir).lower())
                    if same:
                        self.output_status_label.config(
                            text="⚠ Same as input folder: original files will be OVERWRITTEN", fg="red")
                    elif self.count_fasta_files(path) > 0:
                        self.output_status_label.config(
                            text="⚠ Folder exists and contains files (existing files with the same "
                                 "name will be overwritten)", fg="orange")
                    else:
                        self.output_status_label.config(text="✓ Folder is valid", fg="green")
                else:
                    self.output_status_label.config(text="✓ Folder will be created on run", fg="green")

    # ------------------------------------------------------------------ #
    # Helpers
    # ------------------------------------------------------------------ #
    @staticmethod
    def count_fasta_files(folder):
        """Count FASTA files directly inside a folder (non-recursive)"""
        try:
            return sum(1 for f in os.listdir(folder)
                       if os.path.isfile(os.path.join(folder, f))
                       and f.lower().endswith(FASTA_EXTS))
        except OSError:
            return 0

    # ------------------------------------------------------------------ #
    # Core logic
    # ------------------------------------------------------------------ #
    def join_sequences(self):
        """Dispatch to the correct routine based on current mode"""
        if self.mode.get() == "single":
            self.join_single_file()
        else:
            self.join_folder()

    def join_single_file(self):
        """Original behaviour: join all sequences of a single FASTA file"""
        # Validate input file path
        input_path = self.input_path_entry.get().strip()
        if not input_path:
            messagebox.showwarning("Warning", "Please enter or select a FASTA input file path!")
            return

        if not self.input_file or not os.path.exists(self.input_file):
            messagebox.showerror("Error", "Invalid input file path! File does not exist or is not accessible.")
            return

        # Validate output file path
        output_path = self.output_path_entry.get().strip()
        if not output_path:
            messagebox.showwarning("Warning", "Please enter or select output file path!")
            return

        if not self.output_file:
            messagebox.showerror("Error", "Invalid output file path!")
            return

        # Verify output directory is writable again
        output_dir = os.path.dirname(output_path)
        if not output_dir:
            output_dir = "."
        if not os.path.exists(output_dir):
            messagebox.showerror("Error", f"Output directory does not exist: {output_dir}")
            return
        if not os.access(output_dir, os.W_OK):
            messagebox.showerror("Error", f"Output directory is not writable: {output_dir}")
            return

        try:
            self.log("\n" + "=" * 60)
            self.log("Starting single-file processing...")
            self.log("=" * 60 + "\n")

            # Read all sequences
            sequences = []
            seq_count = 0
            total_length = 0

            self.log(f"Reading file: {self.input_file}")

            for record in SeqIO.parse(self.input_file, "fasta"):
                sequences.append(record)
                seq_count += 1
                total_length += len(record.seq)
                self.log(f"  - Sequence {seq_count}: {record.id} (Length: {len(record.seq)})")

            if not sequences:
                self.log("❌ Error: No sequences found in file!")
                messagebox.showerror("Error", "No valid FASTA sequences found in file!")
                return

            self.log(f"\nFound {seq_count} sequences total")
            self.log(f"Total sequence length: {total_length} bp")

            # Join all sequences
            self.log("\nJoining sequences...")
            joined_sequence = "".join(str(record.seq) for record in sequences)
            self.log("✓ Sequence joining completed")
            self.log(f"  Joined sequence length: {len(joined_sequence)} bp")

            # Get custom sequence name
            custom_name = self.seq_name_entry.get().strip()
            if not custom_name:
                custom_name = "joined_sequence"
            self.log(f"  New sequence name: {custom_name}")

            # Create output record
            output_record = SeqRecord(
                Seq(joined_sequence),
                id=custom_name,
                description=f"joined_from_{seq_count}_sequences_total_length_{len(joined_sequence)}"
            )

            # Write output file
            self.log(f"\nWriting output file: {self.output_file}")
            SeqIO.write(output_record, self.output_file, "fasta")
            self.log("✓ Output file saved successfully")

            # Show success message
            from os.path import basename
            success_msg = f"""
Processing completed successfully!

Input file: {basename(self.input_file)}
Number of sequences: {seq_count}
Total sequence length: {total_length} bp

Joined sequence length: {len(joined_sequence)} bp
Output file: {basename(self.output_file)}
            """

            self.log("\n" + "=" * 60)
            self.log("Processing successful!")
            self.log("=" * 60)
            messagebox.showinfo("Success", success_msg.strip())

        except Exception as e:
            error_msg = f"Error occurred during processing:\n{str(e)}"
            self.log(f"\n❌ {error_msg}")
            messagebox.showerror("Error", error_msg)

    def join_folder(self):
        """Batch mode: join sequences inside every FASTA file of a folder.

        Each source file (e.g. 1.fasta) produces an output file with the SAME
        name (e.g. 1.fasta) in the output folder, containing ONE long sequence.
        """
        folder = self.input_path_entry.get().strip()
        out_folder = self.output_path_entry.get().strip()
        suffix = self.seq_name_entry.get().strip()

        # ---- Validate ----
        if not folder:
            messagebox.showwarning("Warning", "Please enter or select an input folder!")
            return
        if not os.path.isdir(folder):
            messagebox.showerror("Error", "Input folder does not exist or is not a directory!")
            return
        if not out_folder:
            messagebox.showwarning("Warning", "Please enter or select an output folder!")
            return

        # Collect all FASTA files (snapshot before writing anything)
        fasta_files = sorted(
            os.path.join(folder, f) for f in os.listdir(folder)
            if os.path.isfile(os.path.join(folder, f)) and f.lower().endswith(FASTA_EXTS)
        )
        if not fasta_files:
            messagebox.showerror("Error", "No FASTA files found in the input folder!")
            return

        # Create output folder if missing
        try:
            os.makedirs(out_folder, exist_ok=True)
        except Exception as e:
            messagebox.showerror("Error", f"Cannot create output folder:\n{e}")
            return
        if not os.access(out_folder, os.W_OK):
            messagebox.showerror("Error", f"Output folder is not writable: {out_folder}")
            return

        # Warn if output folder is the same as input folder (overwrites originals)
        same_folder = (os.path.realpath(folder).lower() == os.path.realpath(out_folder).lower())
        if same_folder:
            if not messagebox.askyesno(
                    "Overwrite Confirmation",
                    "The output folder is the SAME as the input folder.\n\n"
                    "All original FASTA files will be REPLACED by their joined versions.\n\n"
                    "Continue?"):
                return

        try:
            self.log("\n" + "=" * 60)
            self.log("Starting batch folder processing...")
            self.log("=" * 60 + "\n")
            self.log(f"Input folder : {folder}")
            self.log(f"Output folder: {out_folder}")
            self.log(f"FASTA file(s): {len(fasta_files)}")
            if suffix:
                self.log(f"Name suffix  : {suffix}")
            self.log("")

            success_count = 0
            total_sources = 0
            total_joined_len = 0
            failed = []

            for fpath in fasta_files:
                fname = os.path.basename(fpath)
                stem, ext = os.path.splitext(fname)
                try:
                    records = list(SeqIO.parse(fpath, "fasta"))
                    if not records:
                        failed.append(fname)
                        self.log(f"  ⚠ {fname}: no sequence records, skipped")
                        continue

                    total_sources += len(records)
                    joined = "".join(str(r.seq) for r in records)
                    total_joined_len += len(joined)

                    seq_id = stem + suffix if suffix else stem
                    out_record = SeqRecord(
                        Seq(joined),
                        id=seq_id,
                        description=f"joined_from_{len(records)}_sequences_length_{len(joined)}"
                    )
                    out_path = os.path.join(out_folder, fname)  # keep same filename
                    SeqIO.write(out_record, out_path, "fasta")

                    success_count += 1
                    self.log(f"  ✓ {fname}: {len(records)} sequence(s) -> "
                             f"1 sequence ({len(joined)} bp, id={seq_id})")
                except Exception as e:
                    failed.append(fname)
                    self.log(f"  ❌ {fname}: {e}")

            # Summary
            self.log("")
            self.log("=" * 60)
            self.log(f"Done! {success_count}/{len(fasta_files)} file(s) joined successfully")
            self.log(f"Total source sequences joined: {total_sources}")
            self.log(f"Total joined length: {total_joined_len} bp")
            if failed:
                self.log(f"Skipped/failed: {', '.join(failed)}")
            self.log("=" * 60)

            if success_count > 0:
                msg = (f"Batch processing completed!\n\n"
                       f"Successfully joined: {success_count}/{len(fasta_files)} file(s)\n"
                       f"Source sequences joined: {total_sources}\n"
                       f"Total joined length: {total_joined_len} bp\n\n"
                       f"Results saved to:\n{out_folder}\n"
                       f"(each file now contains ONE long sequence, "
                       f"original file names kept)")
                if same_folder:
                    msg += "\n\n⚠ Original files were overwritten."
                if failed:
                    msg += f"\n\nSkipped/failed files:\n{', '.join(failed)}"
                messagebox.showinfo("Success", msg)
            else:
                messagebox.showwarning("Warning",
                                       "No files were processed successfully.\n"
                                       f"Skipped/failed:\n{', '.join(failed)}")

        except Exception as e:
            error_msg = f"Error occurred during batch processing:\n{str(e)}"
            self.log(f"\n❌ {error_msg}")
            messagebox.showerror("Error", error_msg)

    # ------------------------------------------------------------------ #
    # Utilities
    # ------------------------------------------------------------------ #
    def log(self, message):
        """Display log in result area"""
        self.result_text.insert(tk.END, message + "\n")
        self.result_text.see(tk.END)
        self.result_text.update()

    def show_help(self):
        """Display help information"""
        help_text = """
FASTA Sequence Joiner v1.1.0 - User Guide

Description:
  Concatenates all sequence fragments of a FASTA file into one continuous
  sequence. Supports single-file mode and batch-folder mode.

Single-file mode:
  1. Step 1: Browse/select the FASTA input file to process
  2. Step 2: Set the new sequence name and output file path
  3. Step 3: Click "Start Joining Sequences"

Batch-folder mode:
  1. Step 1: Browse/select a FOLDER containing FASTA files
     (e.g. 1.fasta, 2.fasta, 3.fasta; each may hold many sequences)
  2. Step 2: Select an OUTPUT folder (created automatically if missing)
     - Optional "Name suffix" is appended to the joined sequence id,
       e.g. suffix "_joined" makes file 1.fasta contain id "1_joined"
  3. Step 3: Click "Start Joining Sequences"
     - Every FASTA file is turned into ONE long sequence
     - Output keeps the SAME file names as the inputs
       (1.fasta, 2.fasta, 3.fasta each contain a single joined sequence)
     - If the output folder equals the input folder you will be asked
       to confirm BEFORE the original files are overwritten

Supported Formats:
  Input : FASTA files (.fasta, .fa, .fna, .ffn, .faa, .frn)
  Output: FASTA format

Notes:
  - Sequences are joined in the order they appear in the input file
  - Sequences are connected directly without separators
  - Joining is by sequence content order, non-matching lengths are allowed

Example (single file "1.fasta" with 3 fragments):
    Seq1: ATGC
    Seq2: GCTA
    Seq3: CGAT
  Output:
    >1
    ATGCGCTACGAT (12 bp)
        """
        messagebox.showinfo("Help", help_text)


def main():
    """Main function"""
    root = tk.Tk()
    app = FastaJoinerApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
