import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import os
from typing import Dict, Tuple


class BatchRenamerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Batch Renamer Tool")
        self.root.geometry("1000x800")
        
        # File paths
        self.name_mapping_file = tk.StringVar()
        self.target_file = tk.StringVar()
        self.output_file = tk.StringVar()

        # Name mapping dictionary
        self.name_mapping = {}

        # Create interface
        self.create_widgets()
        
    def create_widgets(self):
        # Title
        title_frame = tk.Frame(self.root, pady=10)
        title_frame.pack(fill=tk.X)
        tk.Label(title_frame, text="Batch Renamer Tool", font=("Arial", 16, "bold")).pack()

        # File Selection Area
        file_frame = tk.LabelFrame(self.root, text="File Selection", padx=10, pady=10)
        file_frame.pack(fill=tk.X, padx=10, pady=5)
        
        # Name Mapping File Selection
        tk.Label(file_frame, text="Name Mapping File (TXT):").grid(row=0, column=0, sticky="w", pady=5)
        tk.Entry(file_frame, textvariable=self.name_mapping_file, width=50).grid(row=0, column=1, sticky="ew", padx=5)
        tk.Button(file_frame, text="Browse...", command=self.select_name_mapping_file).grid(row=0, column=2, padx=5)

        # Target File Selection
        tk.Label(file_frame, text="Target File:").grid(row=1, column=0, sticky="w", pady=5)
        tk.Entry(file_frame, textvariable=self.target_file, width=50).grid(row=1, column=1, sticky="ew", padx=5)
        tk.Button(file_frame, text="Browse...", command=self.select_target_file).grid(row=1, column=2, padx=5)

        # Output File Selection
        tk.Label(file_frame, text="Output File:").grid(row=2, column=0, sticky="w", pady=5)
        tk.Entry(file_frame, textvariable=self.output_file, width=50).grid(row=2, column=1, sticky="ew", padx=5)
        tk.Button(file_frame, text="Browse...", command=self.select_output_file).grid(row=2, column=2, padx=5)
        
        file_frame.columnconfigure(1, weight=1)
        
        # Action Buttons
        button_frame = tk.Frame(self.root, pady=10)
        button_frame.pack(fill=tk.X)
        tk.Button(button_frame, text="Load Name Mapping", command=self.load_name_mapping,
                 bg="#4CAF50", fg="white", font=("Arial", 10, "bold"), padx=20).pack(side=tk.LEFT, padx=10)
        tk.Button(button_frame, text="Preview Replacement", command=self.preview_replacement,
                 bg="#2196F3", fg="white", font=("Arial", 10, "bold"), padx=20).pack(side=tk.LEFT, padx=10)
        tk.Button(button_frame, text="Execute Replacement", command=self.execute_replacement,
                 bg="#FF9800", fg="white", font=("Arial", 10, "bold"), padx=20).pack(side=tk.LEFT, padx=10)

        # Preview Area
        preview_frame = tk.LabelFrame(self.root, text="Preview Area", padx=10, pady=10)
        preview_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        # Create split-screen display
        paned_window = tk.PanedWindow(preview_frame, orient=tk.HORIZONTAL)
        paned_window.pack(fill=tk.BOTH, expand=True)
        
        # Left: Original File Preview
        left_frame = tk.Frame(paned_window)
        tk.Label(left_frame, text="Original File Preview", font=("Arial", 10, "bold")).pack(anchor="w")
        self.original_text = scrolledtext.ScrolledText(left_frame, wrap=tk.WORD, height=20)
        self.original_text.pack(fill=tk.BOTH, expand=True)
        paned_window.add(left_frame, minsize=400)

        # Right: Preview After Replacement
        right_frame = tk.Frame(paned_window)
        tk.Label(right_frame, text="Preview After Replacement", font=("Arial", 10, "bold")).pack(anchor="w")
        self.preview_text = scrolledtext.ScrolledText(right_frame, wrap=tk.WORD, height=20)
        self.preview_text.pack(fill=tk.BOTH, expand=True)
        paned_window.add(right_frame, minsize=400)

        # Status Bar
        self.status_var = tk.StringVar(value="Ready")
        status_bar = tk.Label(self.root, textvariable=self.status_var, relief=tk.SUNKEN, anchor="w")
        status_bar.pack(fill=tk.X, side=tk.BOTTOM)

        # Statistics
        self.stats_var = tk.StringVar(value="")
        stats_frame = tk.Frame(self.root, pady=5)
        stats_frame.pack(fill=tk.X)
        tk.Label(stats_frame, textvariable=self.stats_var, fg="blue").pack()
    
    def select_name_mapping_file(self):
        filepath = filedialog.askopenfilename(
            title="Select Name Mapping File",
            filetypes=[("Text Files", "*.txt"), ("All Files", "*.*")]
        )
        if filepath:
            self.name_mapping_file.set(filepath)
            # Auto load
            self.load_name_mapping()

    def select_target_file(self):
        filepath = filedialog.askopenfilename(
            title="Select Target File",
            filetypes=[("All Files", "*.*")]
        )
        if filepath:
            self.target_file.set(filepath)
            # Auto preview
            self.preview_target_file()
            # If mapping loaded, auto preview replacement
            if self.name_mapping:
                self.preview_replacement()

    def select_output_file(self):
        filepath = filedialog.asksaveasfilename(
            title="Select Output File",
            filetypes=[("Text Files", "*.txt"), ("All Files", "*.*")],
            defaultextension=".txt"
        )
        if filepath:
            self.output_file.set(filepath)
    
    def load_name_mapping(self):
        mapping_file = self.name_mapping_file.get()
        if not mapping_file:
            messagebox.showwarning("Warning", "Please select a name mapping file first")
            return

        if not os.path.exists(mapping_file):
            messagebox.showerror("Error", "File does not exist")
            return

        try:
            self.name_mapping = {}
            with open(mapping_file, 'r', encoding='utf-8') as f:
                lines = f.readlines()

            for line in lines:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) >= 2:
                    old_name = parts[0].strip()
                    new_name = parts[1].strip()
                    if old_name and new_name:
                        self.name_mapping[old_name] = new_name

            if self.name_mapping:
                self.status_var.set(f"Successfully loaded {len(self.name_mapping)} name mappings")
                messagebox.showinfo("Success", f"Successfully loaded {len(self.name_mapping)} name mappings")

                # If target file selected, auto preview
                if self.target_file.get():
                    self.preview_replacement()
            else:
                messagebox.showwarning("Warning", "No valid name mappings found")

        except Exception as e:
            messagebox.showerror("Error", f"Failed to load file: {str(e)}")
    
    def preview_target_file(self):
        target_file = self.target_file.get()
        if not target_file or not os.path.exists(target_file):
            return

        try:
            with open(target_file, 'r', encoding='utf-8') as f:
                content = f.read()

            self.original_text.delete(1.0, tk.END)
            self.original_text.insert(1.0, content)

        except Exception as e:
            messagebox.showerror("Error", f"Failed to read target file: {str(e)}")
    
    def preview_replacement(self):
        if not self.name_mapping:
            messagebox.showwarning("Warning", "Please load name mapping first")
            return

        target_file = self.target_file.get()
        if not target_file:
            messagebox.showwarning("Warning", "Please select target file first")
            return

        if not os.path.exists(target_file):
            messagebox.showerror("Error", "Target file does not exist")
            return

        try:
            with open(target_file, 'r', encoding='utf-8') as f:
                content = f.read()

            # Display original content
            self.original_text.delete(1.0, tk.END)
            self.original_text.insert(1.0, content)

            # Execute replacement preview
            preview_content, replaced_count = self.replace_names(content)

            # Display preview result
            self.preview_text.delete(1.0, tk.END)
            self.preview_text.insert(1.0, preview_content)

            # Update statistics
            self.stats_var.set(f"Preview: Replaced {replaced_count} occurrences"
                             f" | Name mappings: {len(self.name_mapping)} items")
            self.status_var.set("Preview completed")

        except Exception as e:
            messagebox.showerror("Error", f"Failed to preview replacement: {str(e)}")
    
    def replace_names(self, content: str) -> Tuple[str, int]:
        """
        Replace names and return replaced content and replacement count
        """
        replaced_count = 0
        result = content

        # Sort by old name length, longest to shortest, to avoid partial matching issues
        sorted_names = sorted(self.name_mapping.keys(), key=len, reverse=True)
        
        for old_name in sorted_names:
            new_name = self.name_mapping[old_name]
            if old_name in result:
                count = result.count(old_name)
                result = result.replace(old_name, new_name)
                replaced_count += count
        
        return result, replaced_count
    
    def execute_replacement(self):
        if not self.name_mapping:
            messagebox.showwarning("Warning", "Please load name mapping first")
            return

        target_file = self.target_file.get()
        if not target_file:
            messagebox.showwarning("Warning", "Please select target file first")
            return

        output_file = self.output_file.get()
        if not output_file:
            messagebox.showwarning("Warning", "Please select output file first")
            return

        if not os.path.exists(target_file):
            messagebox.showerror("Error", "Target file does not exist")
            return

        try:
            # Confirmation dialog
            result = messagebox.askyesno(
                "Confirm",
                f"About to execute replacement operation:\n"
                f"Input file: {target_file}\n"
                f"Output file: {output_file}\n"
                f"Name mappings: {len(self.name_mapping)} items\n\n"
                f"Continue?"
            )

            if not result:
                return

            # Read target file
            with open(target_file, 'r', encoding='utf-8') as f:
                content = f.read()

            # Execute replacement
            new_content, replaced_count = self.replace_names(content)

            # Write to output file
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(new_content)

            # Update preview
            self.preview_text.delete(1.0, tk.END)
            self.preview_text.insert(1.0, new_content)

            # Update statistics
            self.stats_var.set(f"Completed: Replaced {replaced_count} occurrences"
                             f" | Name mappings: {len(self.name_mapping)} items")
            self.status_var.set(f"Replacement completed! File saved to: {output_file}")

            messagebox.showinfo("Success", f"Replacement completed!\n\nTotal replaced: {replaced_count} occurrences\nFile saved to:\n{output_file}")

        except Exception as e:
            messagebox.showerror("Error", f"Failed to execute replacement: {str(e)}")


def main():
    root = tk.Tk()
    app = BatchRenamerApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
