def convert_to_uppercase(input_file, output_file):
    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            content = f.read()
        
        upper_content = content.upper()
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(upper_content)
        
        print(f"转换成功！")
        print(f"输入文件: {input_file}")
        print(f"输出文件: {output_file}")
        
    except FileNotFoundError:
        print(f"错误: 找不到文件 '{input_file}'")
    except Exception as e:
        print(f"发生错误: {e}")


def main():
    print("=== 文本文件大写转换工具 ===")
    print()
    
    input_file = input("请输入要转换的文本文件路径: ").strip()
    
    input_file = input_file.strip('"\'')
    
    import os
    filename = os.path.basename(input_file)
    dirname = os.path.dirname(input_file)
    if dirname:
        output_file = os.path.join(dirname, f"uppercase_{filename}")
    else:
        output_file = f"uppercase_{filename}"
    
    custom_output = input(f"输出文件名 (默认: {output_file}): ").strip()
    if custom_output:
        custom_output = custom_output.strip('"\'')
        output_file = custom_output
    
    print()
    convert_to_uppercase(input_file, output_file)
    print()
    input("按回车键退出...")


if __name__ == "__main__":
    main()
