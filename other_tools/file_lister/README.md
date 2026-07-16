# File Lister

A lightweight browser-based tool to list all files in a folder as a sortable, searchable table.

---

## Files

| File | Language | Description |
|------|:--------:|-------------|
| `文件列表查看器.html` | Chinese | Feature-complete Chinese version |
| `file_lister.html` | English | Feature-complete English version |

---

## How to Use

1. Double-click either `.html` file to open it in your browser.
2. Click **Select Folder** or **drag & drop** a folder onto the page.
3. All files (with extensions) will be displayed in a table.

### Features

| Feature | Action |
|---------|--------|
| Sort | Click any column header (toggles asc/desc) |
| Search/Filter | Type a keyword in the search box |
| Multi-select | `Ctrl + Click` rows |
| Select all | `Ctrl + A` |
| Copy selected | `Ctrl + C` or right-click context menu |
| Copy all | Click the **Copy All** button |
| Export CSV | Click the **Export CSV** button (opens in Excel) |

### Table Columns

| Column | Description |
|--------|-------------|
| # | Auto-numbered |
| Filename | Full filename with extension |
| Extension | File suffix (e.g. .docx, .xlsx) |
| Size | File size in B/KB/MB |
| Date Modified | Last modified timestamp |

---

## FAQ

**Q: Nothing happens when I select a folder?**  
A: Use Chrome, Edge, or Firefox. If your browser is outdated, try the drag-and-drop method.

**Q: The browser asks for permission to access files?**  
A: This is normal. Click "Allow" or "Upload" — your files are processed locally and never sent anywhere.

---

## Requirements

- Chrome 86+ / Edge 86+ / Firefox 85+
- No installation, no internet connection needed
