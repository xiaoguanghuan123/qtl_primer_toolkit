#!/usr/bin/env python3
"""
依赖安装脚本 - 自动处理依赖问题
"""
import subprocess
import sys
import os

def run_command(cmd, description):
    """运行命令并处理错误"""
    print(f"正在{description}...")
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"警告: {description}失败")
            print(f"错误信息: {result.stderr[:200]}")
            return False
        print(f"{description}成功")
        return True
    except Exception as e:
        print(f"{description}异常: {e}")
        return False

def main():
    """主安装函数"""
    
    print("=" * 60)
    print("QTL Primer Toolkit 依赖安装")
    print("=" * 60)
    
    # 检查Python版本
    if sys.version_info < (3, 8):
        print("错误: 需要Python 3.8或更高版本")
        return
    
    # 1. 升级pip
    run_command("python -m pip install --upgrade pip", "升级pip")
    
    # 2. 尝试标准安装
    print("\n1. 尝试安装标准依赖...")
    if os.path.exists("requirements.txt"):
        success = run_command("pip install -r requirements.txt", "安装标准依赖")
        if not success:
            print("\n2. 标准安装失败，尝试最小化安装...")
            if os.path.exists("requirements-minimal.txt"):
                run_command("pip install -r requirements-minimal.txt", "安装最小依赖")
            else:
                # 直接安装核心包
                core_packages = [
                    "primer3-py>=2.0.0",
                    "biopython>=1.79",
                    "pandas>=1.3.0",
                    "matplotlib>=3.4.0",
                    "numpy>=1.21.0",
                    "jinja2>=3.0.0",
                    "openpyxl>=3.0.0",
                    "pyyaml>=6.0"
                ]
                for pkg in core_packages:
                    run_command(f"pip install {pkg}", f"安装 {pkg.split('>')[0]}")
    
    # 3. 尝试使用国内镜像
    print("\n3. 如果需要，可尝试使用国内镜像源安装:")
    print("   pip install -r requirements.txt -i https://pypi.tuna.tsinghua.edu.cn/simple")
    
    # 4. 验证安装
    print("\n4. 验证安装...")
    test_imports = [
        ("primer3", "primer3-py"),
        ("Bio", "biopython"),
        ("pandas", "pandas"),
        ("matplotlib", "matplotlib"),
        ("numpy", "numpy"),
        ("jinja2", "jinja2"),
        ("openpyxl", "openpyxl"),
        ("yaml", "pyyaml"),
    ]
    
    all_success = True
    for module, package in test_imports:
        try:
            __import__(module)
            print(f"  ✓ {package} 安装成功")
        except ImportError:
            print(f"  ✗ {package} 安装失败")
            all_success = False
    
    print("\n" + "=" * 60)
    if all_success:
        print("✅ 所有依赖安装成功！")
        print("\n下一步：")
        print("1. 安装工具包: pip install -e .")
        print("2. 测试安装: qtlprimer --version")
    else:
        print("⚠️  部分依赖安装失败，请手动安装缺少的包")
        print("\n手动安装命令示例:")
        print("  pip install biopython pandas matplotlib numpy")
    
    print("=" * 60)

if __name__ == "__main__":
    main()
