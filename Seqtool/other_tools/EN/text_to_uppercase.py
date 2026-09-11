def convert_to_uppercase(input_file, output_file):
    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            content = f.read()
        
        upper_content = content.upper()
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(upper_content)
        
        print(f"Conversion successful!")
        print(f"Input file: {input_file}")
        print(f"Output file: {output_file}")

    except FileNotFoundError:
        print(f"Error: File not found '{input_file}'")
    except Exception as e:
        print(f"An error occurred: {e}")


def main():
    print("=== Text File Uppercase Converter ===")
    print()

    input_file = input("Enter the path of the text file to convert: ").strip()
    
    input_file = input_file.strip('"\'')
    
    import os
    filename = os.path.basename(input_file)
    dirname = os.path.dirname(input_file)
    if dirname:
        output_file = os.path.join(dirname, f"uppercase_{filename}")
    else:
        output_file = f"uppercase_{filename}"
    
    custom_output = input(f"Output filename (default: {output_file}): ").strip()
    if custom_output:
        custom_output = custom_output.strip('"\'')
        output_file = custom_output
    
    print()
    convert_to_uppercase(input_file, output_file)
    print()
    input("Press Enter to exit...")


if __name__ == "__main__":
    main()
