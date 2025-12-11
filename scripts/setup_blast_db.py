#!/usr/bin/env python3
"""
BLAST数据库设置脚本
"""
import subprocess
import argparse
from pathlib import Path

def setup_blast_db(fasta_file: str, db_name: str = None, db_type: str = "nucl"):
    """
    设置BLAST数据库
    
    Args:
        fasta_file: FASTA文件路径
        db_name: 数据库名称（默认与FASTA文件同名）
        db_type: 数据库类型 (nucl/prot)
    """
    fasta_path = Path(fasta_file)
    
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA文件不存在: {fasta_file}")
    
    if db_name is None:
        db_name = fasta_path.stem
    
    # 构建makeblastdb命令
    cmd = [
        "makeblastdb",
        "-in", str(fasta_path),
        "-dbtype", db_type,
        "-out", db_name,
        "-parse_seqids"
    ]
    
    print(f"正在创建BLAST数据库: {db_name}")
    print(f"命令: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("✅ BLAST数据库创建成功")
        print(f"数据库文件: {db_name}.*")
        
        # 验证数据库
        verify_cmd = ["blastdbcmd", "-db", db_name, "-info"]
        subprocess.run(verify_cmd, capture_output=True, check=True)
        print("✅ 数据库验证通过")
        
    except subprocess.CalledProcessError as e:
        print(f"❌ 创建BLAST数据库失败:")
        print(f"错误信息: {e.stderr}")
        raise
    except FileNotFoundError:
        print("❌ 找不到BLAST工具，请确保已安装NCBI BLAST+")
        print("安装方法:")
        print("  Ubuntu: sudo apt-get install ncbi-blast+")
        print("  macOS: brew install blast")
        print("  Conda: conda install -c bioconda blast")
        raise

def main():
    parser = argparse.ArgumentParser(description="设置BLAST数据库")
    parser.add_argument("fasta", help="FASTA文件路径")
    parser.add_argument("--db-name", help="数据库名称")
    parser.add_argument("--type", choices=["nucl", "prot"], default="nucl",
                       help="数据库类型 (默认: nucl)")
    
    args = parser.parse_args()
    
    try:
        setup_blast_db(args.fasta, args.db_name, args.type)
    except Exception as e:
        print(f"失败: {str(e)}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
