"""
核心工具函数模块
提供序列处理、文件操作、计算和验证等通用功能。
"""
import gzip
import hashlib
import json
import logging
import math
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Any, BinaryIO
from dataclasses import dataclass
from datetime import datetime

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# ============================================================================
# 序列处理工具
# ============================================================================

class SequenceUtils:
    """DNA/RNA序列处理工具类"""
    
    # DNA碱基互补表
    DNA_COMPLEMENT = str.maketrans('ACGTacgtNn-', 'TGCAtgcaNn-')
    RNA_COMPLEMENT = str.maketrans('ACGUacguNn-', 'UGCAugcaNn-')
    
    # 有效碱基集合
    VALID_DNA_BASES = set('ACGTNacgtn-')
    VALID_RNA_BASES = set('ACGUNacgun-')
    
    @staticmethod
    def reverse_complement(sequence: str, is_rna: bool = False) -> str:
        """
        计算DNA/RNA序列的反向互补序列
        
        Args:
            sequence: 输入序列
            is_rna: 是否为RNA序列
            
        Returns:
            str: 反向互补序列
        """
        if not sequence:
            return ""
        
        if is_rna:
            return sequence.translate(SequenceUtils.RNA_COMPLEMENT)[::-1]
        else:
            return sequence.translate(SequenceUtils.DNA_COMPLEMENT)[::-1]
    
    @staticmethod
    def complement(sequence: str, is_rna: bool = False) -> str:
        """
        计算DNA/RNA序列的互补序列（不反向）
        
        Args:
            sequence: 输入序列
            is_rna: 是否为RNA序列
            
        Returns:
            str: 互补序列
        """
        if not sequence:
            return ""
        
        if is_rna:
            return sequence.translate(SequenceUtils.RNA_COMPLEMENT)
        else:
            return sequence.translate(SequenceUtils.DNA_COMPLEMENT)
    
    @staticmethod
    def is_valid_dna(sequence: str) -> bool:
        """检查序列是否为有效DNA序列"""
        return all(base in SequenceUtils.VALID_DNA_BASES for base in sequence.upper())
    
    @staticmethod
    def is_valid_rna(sequence: str) -> bool:
        """检查序列是否为有效RNA序列"""
        return all(base in SequenceUtils.VALID_RNA_BASES for base in sequence.upper())
    
    @staticmethod
    def calculate_gc_content(sequence: str) -> float:
        """
        计算DNA序列的GC含量（百分比）
        
        Args:
            sequence: DNA序列
            
        Returns:
            float: GC含量（0-100）
        """
        if not sequence:
            return 0.0
        
        seq_upper = sequence.upper()
        g_count = seq_upper.count('G')
        c_count = seq_upper.count('C')
        total = len(sequence)
        
        if total == 0:
            return 0.0
        
        return (g_count + c_count) / total * 100.0
    
    @staticmethod
    def calculate_tm_wallace(sequence: str) -> float:
        """
        使用Wallace规则计算DNA序列的Tm值
        
        Args:
            sequence: DNA序列
            
        Returns:
            float: Tm值（摄氏度）
        """
        if not sequence:
            return 0.0
        
        seq_upper = sequence.upper()
        a_count = seq_upper.count('A')
        t_count = seq_upper.count('T')
        g_count = seq_upper.count('G')
        c_count = seq_upper.count('C')
        
        # Wallace规则: Tm = 2*(A+T) + 4*(G+C)
        return 2.0 * (a_count + t_count) + 4.0 * (g_count + c_count)
    
    @staticmethod
    def calculate_tm_nn(sequence: str, salt_conc: float = 0.05, primer_conc: float = 0.0000005) -> float:
        """
        使用最近邻（Nearest Neighbor）方法计算Tm值
        
        Args:
            sequence: DNA序列
            salt_conc: 盐浓度（M）
            primer_conc: 引物浓度（M）
            
        Returns:
            float: Tm值（摄氏度）
        """
        if len(sequence) < 4:
            return SequenceUtils.calculate_tm_wallace(sequence)
        
        # 简单的最近邻参数（更复杂的实现可使用primer3）
        # 基于SantaLucia (1998)参数
        nn_params = {
            'AA': -9.1, 'TT': -9.1,
            'AT': -8.6,
            'TA': -6.0,
            'CA': -5.8, 'TG': -5.8,
            'GT': -6.5, 'AC': -6.5,
            'CT': -7.8, 'AG': -7.8,
            'GA': -5.6, 'TC': -5.6,
            'CG': -11.9,
            'GC': -11.1,
            'GG': -11.0, 'CC': -11.0
        }
        
        # 初始化ΔH和ΔS
        delta_h = 0.0
        delta_s = 0.0
        
        # 计算最近邻贡献
        seq_upper = sequence.upper()
        for i in range(len(seq_upper) - 1):
            dimer = seq_upper[i:i+2]
            if dimer in nn_params:
                delta_h += nn_params[dimer]
                # 简化的ΔS估计
                delta_s += -0.024
        
        # 对称性修正
        if seq_upper[0] == 'G' or seq_upper[0] == 'C':
            delta_h += 0.1
        if seq_upper[-1] == 'G' or seq_upper[-1] == 'C':
            delta_h += 0.1
        
        # 初始化修正
        delta_h += 0.2
        delta_s += -0.0057
        
        # 计算Tm
        R = 1.987  # 气体常数
        if delta_s >= 0:
            return 0.0
        
        tm_kelvin = delta_h * 1000 / (delta_s + R * math.log(primer_conc / 4))
        tm_celsius = tm_kelvin - 273.15
        
        # 盐浓度修正
        if salt_conc > 0:
            tm_celsius += 16.6 * math.log10(salt_conc)
        
        return tm_celsius
    
    @staticmethod
    def calculate_self_complementarity(sequence: str) -> float:
        """
        计算序列的自互补性评分
        
        Args:
            sequence: DNA序列
            
        Returns:
            float: 自互补性评分（越高表示越容易形成二聚体）
        """
        if len(sequence) < 4:
            return 0.0
        
        score = 0.0
        seq_len = len(sequence)
        
        # 检查反向互补重叠
        for i in range(seq_len - 3):
            for j in range(i + 4, seq_len):
                subseq1 = sequence[i:j]
                subseq2 = sequence[j:] if j < seq_len else ""
                
                if not subseq2:
                    break
                
                # 检查反向互补
                rc_subseq1 = SequenceUtils.reverse_complement(subseq1)
                if rc_subseq1 in sequence:
                    score += len(subseq1) * 1.5
        
        # 检查回文结构
        for i in range(seq_len - 6):
            for length in range(6, min(20, seq_len - i)):
                subseq = sequence[i:i+length]
                if subseq == SequenceUtils.reverse_complement(subseq):
                    score += length * 2.0
        
        return score
    
    @staticmethod
    def calculate_hairpin_score(sequence: str) -> float:
        """
        计算发夹结构形成潜力
        
        Args:
            sequence: DNA序列
            
        Returns:
            float: 发夹结构评分（越高表示越容易形成发夹）
        """
        seq_len = len(sequence)
        if seq_len < 10:
            return 0.0
        
        max_score = 0.0
        
        # 检查可能的发夹结构
        for loop_size in range(3, 10):
            for i in range(seq_len - 2 * loop_size - 4):
                stem1 = sequence[i:i+loop_size]
                loop = sequence[i+loop_size:i+loop_size+4]
                stem2 = sequence[i+loop_size+4:i+2*loop_size+4]
                
                if len(stem2) < loop_size:
                    continue
                
                # 计算茎区的互补性
                rc_stem1 = SequenceUtils.reverse_complement(stem1)
                matches = sum(1 for a, b in zip(rc_stem1, stem2) if a == b)
                match_score = matches / loop_size
                
                # 环区惩罚（G/C丰富增加稳定性）
                loop_gc = SequenceUtils.calculate_gc_content(loop)
                loop_score = loop_gc / 100.0
                
                total_score = match_score * 2.0 + loop_score
                max_score = max(max_score, total_score)
        
        return max_score
    
    @staticmethod
    def find_repeats(sequence: str, min_repeat: int = 3) -> List[Dict]:
        """
        查找序列中的重复模式
        
        Args:
            sequence: DNA序列
            min_repeat: 最小重复次数
            
        Returns:
            List[Dict]: 重复模式列表
        """
        repeats = []
        seq_len = len(sequence)
        
        # 查找单碱基重复
        for base in ['A', 'T', 'G', 'C']:
            pattern = base * min_repeat
            pos = 0
            while pos < seq_len:
                match_start = sequence.find(pattern, pos)
                if match_start == -1:
                    break
                
                # 扩展重复
                match_end = match_start + min_repeat
                while match_end < seq_len and sequence[match_end] == base:
                    match_end += 1
                
                repeat_length = match_end - match_start
                if repeat_length >= min_repeat:
                    repeats.append({
                        'start': match_start,
                        'end': match_end,
                        'pattern': base,
                        'length': repeat_length,
                        'type': 'mononucleotide'
                    })
                
                pos = match_end
        
        return repeats

# ============================================================================
# 文件处理工具
# ============================================================================

class FileUtils:
    """文件处理工具类"""
    
    @staticmethod
    def ensure_directory(path: Union[str, Path]) -> Path:
        """
        确保目录存在，如果不存在则创建
        
        Args:
            path: 目录路径
            
        Returns:
            Path: 目录Path对象
        """
        path_obj = Path(path)
        path_obj.mkdir(parents=True, exist_ok=True)
        return path_obj
    
    @staticmethod
    def safe_remove(path: Union[str, Path]) -> bool:
        """
        安全删除文件或目录
        
        Args:
            path: 路径
            
        Returns:
            bool: 是否成功删除
        """
        try:
            path_obj = Path(path)
            if path_obj.exists():
                if path_obj.is_file():
                    path_obj.unlink()
                else:
                    shutil.rmtree(path_obj)
                logger.debug(f"已删除: {path}")
                return True
        except Exception as e:
            logger.warning(f"删除失败 {path}: {e}")
        return False
    
    @staticmethod
    def read_file(file_path: Union[str, Path], encoding: str = 'utf-8') -> str:
        """
        读取文本文件，支持gzip压缩
        
        Args:
            file_path: 文件路径
            encoding: 文件编码
            
        Returns:
            str: 文件内容
        """
        file_path = Path(file_path)
        
        if not file_path.exists():
            raise FileNotFoundError(f"文件不存在: {file_path}")
        
        try:
            if file_path.suffix == '.gz':
                with gzip.open(file_path, 'rt', encoding=encoding) as f:
                    return f.read()
            else:
                with open(file_path, 'r', encoding=encoding) as f:
                    return f.read()
        except Exception as e:
            logger.error(f"读取文件失败 {file_path}: {e}")
            raise
    
    @staticmethod
    def write_file(content: str, file_path: Union[str, Path], 
                  encoding: str = 'utf-8', compress: bool = False) -> Path:
        """
        写入文本文件，可选gzip压缩
        
        Args:
            content: 文件内容
            file_path: 文件路径
            encoding: 文件编码
            compress: 是否压缩
            
        Returns:
            Path: 写入的文件路径
        """
        file_path = Path(file_path)
        
        # 确保目录存在
        file_path.parent.mkdir(parents=True, exist_ok=True)
        
        try:
            if compress or file_path.suffix == '.gz':
                if not file_path.suffix == '.gz':
                    file_path = file_path.with_suffix(file_path.suffix + '.gz')
                
                with gzip.open(file_path, 'wt', encoding=encoding) as f:
                    f.write(content)
            else:
                with open(file_path, 'w', encoding=encoding) as f:
                    f.write(content)
            
            logger.debug(f"文件已写入: {file_path}")
            return file_path
            
        except Exception as e:
            logger.error(f"写入文件失败 {file_path}: {e}")
            raise
    
    @staticmethod
    def read_fasta(file_path: Union[str, Path]) -> Dict[str, str]:
        """
        读取FASTA文件
        
        Args:
            file_path: FASTA文件路径
            
        Returns:
            Dict[str, str]: 序列字典 {header: sequence}
        """
        content = FileUtils.read_file(file_path)
        sequences = {}
        current_header = None
        current_sequence = []
        
        for line in content.splitlines():
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                # 保存前一个序列
                if current_header is not None:
                    sequences[current_header] = ''.join(current_sequence)
                
                # 开始新序列
                current_header = line[1:].split()[0]  # 取第一个词作为header
                current_sequence = []
            else:
                current_sequence.append(line)
        
        # 保存最后一个序列
        if current_header is not None:
            sequences[current_header] = ''.join(current_sequence)
        
        return sequences
    
    @staticmethod
    def write_fasta(sequences: Dict[str, str], 
                   file_path: Union[str, Path],
                   line_length: int = 80) -> Path:
        """
        写入FASTA文件
        
        Args:
            sequences: 序列字典 {header: sequence}
            file_path: 输出文件路径
            line_length: 每行序列长度
            
        Returns:
            Path: 写入的文件路径
        """
        lines = []
        for header, sequence in sequences.items():
            lines.append(f'>{header}')
            
            # 分割序列为多行
            for i in range(0, len(sequence), line_length):
                lines.append(sequence[i:i+line_length])
        
        content = '\n'.join(lines)
        return FileUtils.write_file(content, file_path)
    
    @staticmethod
    def get_file_hash(file_path: Union[str, Path], 
                     algorithm: str = 'md5') -> str:
        """
        计算文件哈希值
        
        Args:
            file_path: 文件路径
            algorithm: 哈希算法 ('md5', 'sha1', 'sha256')
            
        Returns:
            str: 文件哈希值
        """
        file_path = Path(file_path)
        
        if not file_path.exists():
            raise FileNotFoundError(f"文件不存在: {file_path}")
        
        hash_func = hashlib.new(algorithm)
        
        try:
            with open(file_path, 'rb') as f:
                # 分块读取大文件
                for chunk in iter(lambda: f.read(4096), b''):
                    hash_func.update(chunk)
        except Exception as e:
            logger.error(f"计算文件哈希失败 {file_path}: {e}")
            raise
        
        return hash_func.hexdigest()
    
    @staticmethod
    def validate_fasta_file(file_path: Union[str, Path]) -> Tuple[bool, List[str]]:
        """
        验证FASTA文件格式
        
        Args:
            file_path: FASTA文件路径
            
        Returns:
            Tuple[bool, List[str]]: (是否有效, 错误信息列表)
        """
        errors = []
        
        try:
            sequences = FileUtils.read_fasta(file_path)
            
            if not sequences:
                errors.append("FASTA文件为空")
                return False, errors
            
            # 检查序列有效性
            invalid_sequences = []
            for header, seq in sequences.items():
                if not seq:
                    invalid_sequences.append(header)
                elif not SequenceUtils.is_valid_dna(seq) and not SequenceUtils.is_valid_rna(seq):
                    invalid_sequences.append(header)
            
            if invalid_sequences:
                errors.append(f"发现无效序列: {', '.join(invalid_sequences[:5])}")
                if len(invalid_sequences) > 5:
                    errors.append(f"... 共 {len(invalid_sequences)} 个无效序列")
            
            return len(errors) == 0, errors
            
        except Exception as e:
            errors.append(f"读取FASTA文件失败: {str(e)}")
            return False, errors
    
    @staticmethod
    def create_temp_file(content: str = None, suffix: str = '.tmp') -> Path:
        """
        创建临时文件
        
        Args:
            content: 文件内容（可选）
            suffix: 文件后缀
            
        Returns:
            Path: 临时文件路径
        """
        temp_dir = tempfile.mkdtemp(prefix='qtlprimer_')
        temp_file = Path(temp_dir) / f"temp_{datetime.now().strftime('%Y%m%d_%H%M%S')}{suffix}"
        
        if content is not None:
            FileUtils.write_file(content, temp_file)
        
        return temp_file

# ============================================================================
# BLAST工具
# ============================================================================

class BlastUtils:
    """BLAST相关工具函数"""
    
    @staticmethod
    def validate_blast_database(db_path: Union[str, Path]) -> Tuple[bool, List[str]]:
        """
        验证BLAST数据库是否有效
        
        Args:
            db_path: BLAST数据库路径（不含后缀）
            
        Returns:
            Tuple[bool, List[str]]: (是否有效, 错误信息列表)
        """
        errors = []
        db_path = Path(db_path)
        
        # 检查必需的文件
        required_extensions = ['.nhr', '.nin', '.nsq']
        missing_files = []
        
        for ext in required_extensions:
            db_file = f"{db_path}{ext}"
            if not Path(db_file).exists():
                missing_files.append(db_file)
        
        if missing_files:
            errors.append(f"BLAST数据库文件缺失: {', '.join(missing_files)}")
            return False, errors
        
        # 尝试运行blastdbcmd验证数据库
        try:
            cmd = ['blastdbcmd', '-db', str(db_path), '-info']
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            
            if result.returncode != 0:
                errors.append(f"BLAST数据库验证失败: {result.stderr}")
                return False, errors
            
            logger.debug(f"BLAST数据库验证通过: {db_path}")
            return True, errors
            
        except FileNotFoundError:
            errors.append("未找到blastdbcmd命令，请确保已安装NCBI BLAST+")
            return False, errors
        except subprocess.TimeoutExpired:
            errors.append("BLAST数据库验证超时")
            return False, errors
        except Exception as e:
            errors.append(f"BLAST数据库验证异常: {str(e)}")
            return False, errors
    
    @staticmethod
    def create_blast_database(fasta_file: Union[str, Path], 
                             db_name: str = None,
                             db_type: str = 'nucl',
                             parse_seqids: bool = True) -> Path:
        """
        创建BLAST数据库
        
        Args:
            fasta_file: FASTA文件路径
            db_name: 数据库名称（默认与FASTA文件同名）
            db_type: 数据库类型 ('nucl' 或 'prot')
            parse_seqids: 是否解析序列ID
            
        Returns:
            Path: 创建的数据库路径
        """
        fasta_path = Path(fasta_file)
        
        if not fasta_path.exists():
            raise FileNotFoundError(f"FASTA文件不存在: {fasta_file}")
        
        if db_name is None:
            db_name = fasta_path.stem
        
        db_path = Path(db_name)
        
        # 构建makeblastdb命令
        cmd = ['makeblastdb', '-in', str(fasta_path), '-dbtype', db_type, '-out', str(db_path)]
        
        if parse_seqids:
            cmd.append('-parse_seqids')
        
        logger.info(f"创建BLAST数据库: {db_path}")
        logger.debug(f"命令: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            logger.info("BLAST数据库创建成功")
            logger.debug(f"输出: {result.stdout}")
            
            # 验证数据库
            is_valid, errors = BlastUtils.validate_blast_database(db_path)
            if not is_valid:
                raise RuntimeError(f"数据库验证失败: {errors}")
            
            return db_path
            
        except subprocess.CalledProcessError as e:
            logger.error(f"创建BLAST数据库失败: {e.stderr}")
            raise RuntimeError(f"BLAST数据库创建失败: {e.stderr}")
        except Exception as e:
            logger.error(f"创建BLAST数据库异常: {str(e)}")
            raise
    
    @staticmethod
    def run_blast(query_file: Union[str, Path],
                 db_path: Union[str, Path],
                 output_file: Union[str, Path],
                 blast_type: str = 'blastn',
                 evalue: float = 10.0,
                 num_threads: int = 1,
                 **kwargs) -> Tuple[bool, str]:
        """
        运行BLAST搜索
        
        Args:
            query_file: 查询序列文件
            db_path: BLAST数据库路径
            output_file: 输出文件路径
            blast_type: BLAST程序类型
            evalue: E值阈值
            num_threads: 线程数
            **kwargs: 其他BLAST参数
            
        Returns:
            Tuple[bool, str]: (是否成功, 错误信息)
        """
        query_path = Path(query_file)
        db_path = Path(db_path)
        output_path = Path(output_file)
        
        if not query_path.exists():
            return False, f"查询文件不存在: {query_file}"
        
        # 构建BLAST命令
        cmd = [
            blast_type,
            '-query', str(query_path),
            '-db', str(db_path),
            '-out', str(output_path),
            '-outfmt', '6',  # 制表符分隔格式
            '-evalue', str(evalue),
            '-num_threads', str(num_threads)
        ]
        
        # 添加其他参数
        for key, value in kwargs.items():
            if value is not None:
                cmd.extend([f'-{key}', str(value)])
        
        logger.info(f"运行BLAST: {blast_type}")
        logger.debug(f"命令: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            logger.info(f"BLAST运行成功: {output_path}")
            return True, ""
            
        except subprocess.CalledProcessError as e:
            error_msg = f"BLAST运行失败: {e.stderr}"
            logger.error(error_msg)
            return False, error_msg
        except Exception as e:
            error_msg = f"BLAST运行异常: {str(e)}"
            logger.error(error_msg)
            return False, error_msg

# ============================================================================
# 验证工具
# ============================================================================

class ValidationUtils:
    """数据验证工具类"""
    
    @staticmethod
    def validate_genome_file(file_path: Union[str, Path]) -> Tuple[bool, List[str]]:
        """
        验证基因组文件（FASTA格式）
        
        Args:
            file_path: 基因组文件路径
            
        Returns:
            Tuple[bool, List[str]]: (是否有效, 错误信息列表)
        """
        errors = []
        
        try:
            # 检查文件是否存在
            if not Path(file_path).exists():
                errors.append(f"文件不存在: {file_path}")
                return False, errors
            
            # 验证FASTA格式
            is_valid, fasta_errors = FileUtils.validate_fasta_file(file_path)
            if not is_valid:
                errors.extend(fasta_errors)
            
            # 检查序列长度（基因组通常应该较长）
            sequences = FileUtils.read_fasta(file_path)
            short_sequences = []
            for header, seq in sequences.items():
                if len(seq) < 100:
                    short_sequences.append(f"{header}: {len(seq)}bp")
            
            if short_sequences:
                errors.append(f"发现短序列（<100bp）: {', '.join(short_sequences[:3])}")
                if len(short_sequences) > 3:
                    errors.append(f"... 共 {len(short_sequences)} 个短序列")
            
            return len(errors) == 0, errors
            
        except Exception as e:
            errors.append(f"验证基因组文件失败: {str(e)}")
            return False, errors
    
    @staticmethod
    def validate_variant_file(file_path: Union[str, Path]) -> Tuple[bool, List[str]]:
        """
        验证变异文件
        
        Args:
            file_path: 变异文件路径
            
        Returns:
            Tuple[bool, List[str]]: (是否有效, 错误信息列表)
        """
        errors = []
        
        try:
            # 检查文件是否存在
            if not Path(file_path).exists():
                errors.append(f"文件不存在: {file_path}")
                return False, errors
            
            # 尝试读取文件
            from qtlprimer.core.variant_parser import validate_variant_file as validate_vcf
            return validate_vcf(file_path)
            
        except ImportError:
            # 如果variant_parser不可用，进行基本验证
            try:
                content = FileUtils.read_file(file_path)
                lines = content.strip().split('\n')
                
                if not lines:
                    errors.append("文件为空")
                    return False, errors
                
                # 检查是否包含变异数据的关键词
                vcf_keywords = ['#CHROM', 'CHROM', 'POS', 'REF', 'ALT']
                csv_keywords = ['CHR', 'POS', 'REF', 'ALT']
                
                has_vcf = any(keyword in lines[0] for keyword in vcf_keywords)
                has_csv = any(keyword in lines[0].upper() for keyword in csv_keywords)
                
                if not (has_vcf or has_csv):
                    errors.append("文件不包含变异数据所需的列（CHR, POS, REF, ALT）")
                
                return len(errors) == 0, errors
                
            except Exception as e:
                errors.append(f"验证变异文件失败: {str(e)}")
                return False, errors
        except Exception as e:
            errors.append(f"验证变异文件失败: {str(e)}")
            return False, errors
    
    @staticmethod
    def validate_primers(primers: List[Dict]) -> Tuple[bool, List[str]]:
        """
        验证引物列表
        
        Args:
            primers: 引物字典列表
            
        Returns:
            Tuple[bool, List[str]]: (是否有效, 错误信息列表)
        """
        errors = []
        
        if not primers:
            errors.append("引物列表为空")
            return False, errors
        
        for i, primer in enumerate(primers):
            # 检查必需字段
            required_fields = ['sequence', 'name']
            for field in required_fields:
                if field not in primer:
                    errors.append(f"引物 #{i} 缺少字段: {field}")
            
            # 检查序列有效性
            if 'sequence' in primer:
                seq = primer['sequence']
                if not seq:
                    errors.append(f"引物 #{i} 序列为空")
                elif not SequenceUtils.is_valid_dna(seq):
                    errors.append(f"引物 #{i} 序列包含无效字符: {seq}")
            
            # 检查长度
            if 'sequence' in primer and len(primer['sequence']) < 15:
                errors.append(f"引物 #{i} 序列过短 (<15bp): {primer['sequence']}")
            if 'sequence' in primer and len(primer['sequence']) > 30:
                errors.append(f"引物 #{i} 序列过长 (>30bp): {primer['sequence']}")
        
        return len(errors) == 0, errors

# ============================================================================
# 便捷函数
# ============================================================================

def calculate_tm(sequence: str, method: str = 'wallace', **kwargs) -> float:
    """
    计算DNA序列的Tm值（便捷函数）
    
    Args:
        sequence: DNA序列
        method: 计算方法 ('wallace', 'nn')
        **kwargs: 传递给计算方法的参数
        
    Returns:
        float: Tm值
    """
    if method == 'wallace':
        return SequenceUtils.calculate_tm_wallace(sequence)
    elif method == 'nn':
        return SequenceUtils.calculate_tm_nn(sequence, **kwargs)
    else:
        raise ValueError(f"不支持的Tm计算方法: {method}")


def calculate_gc_content(sequence: str) -> float:
    """
    计算DNA序列的GC含量（便捷函数）
    
    Args:
        sequence: DNA序列
        
    Returns:
        float: GC含量（百分比）
    """
    return SequenceUtils.calculate_gc_content(sequence)


def reverse_complement(sequence: str, is_rna: bool = False) -> str:
    """
    计算DNA/RNA序列的反向互补序列（便捷函数）
    
    Args:
        sequence: 输入序列
        is_rna: 是否为RNA序列
        
    Returns:
        str: 反向互补序列
    """
    return SequenceUtils.reverse_complement(sequence, is_rna)


def validate_genome_file(file_path: Union[str, Path]) -> Tuple[bool, List[str]]:
    """
    验证基因组文件（便捷函数）
    
    Args:
        file_path: 基因组文件路径
        
    Returns:
        Tuple[bool, List[str]]: (是否有效, 错误信息列表)
    """
    return ValidationUtils.validate_genome_file(file_path)


def validate_blast_database(db_path: Union[str, Path]) -> Tuple[bool, List[str]]:
    """
    验证BLAST数据库（便捷函数）
    
    Args:
        db_path: BLAST数据库路径
        
    Returns:
        Tuple[bool, List[str]]: (是否有效, 错误信息列表)
    """
    return BlastUtils.validate_blast_database(db_path)


def create_temp_fasta(sequences: Dict[str, str]) -> Path:
    """
    创建临时FASTA文件（便捷函数）
    
    Args:
        sequences: 序列字典 {header: sequence}
        
    Returns:
        Path: 临时文件路径
    """
    temp_file = FileUtils.create_temp_file(suffix='.fasta')
    FileUtils.write_fasta(sequences, temp_file)
    return temp_file


def batch_process(items: List[Any], 
                  process_func: callable,
                  batch_size: int = 100,
                  max_workers: int = 4,
                  **kwargs) -> List[Any]:
    """
    批量处理项目
    
    Args:
        items: 要处理的项目列表
        process_func: 处理函数
        batch_size: 批次大小
        max_workers: 最大工作线程数
        **kwargs: 传递给处理函数的参数
        
    Returns:
        List[Any]: 处理结果列表
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed
    
    results = []
    
    # 分批处理
    for i in range(0, len(items), batch_size):
        batch = items[i:i+batch_size]
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(process_func, item, **kwargs) for item in batch]
            
            for future in as_completed(futures):
                try:
                    result = future.result()
                    results.append(result)
                except Exception as e:
                    logger.error(f"批量处理失败: {e}")
                    results.append(None)
    
    return results


# 导出主要类
__all__ = [
    'SequenceUtils',
    'FileUtils',
    'BlastUtils',
    'ValidationUtils',
    # 便捷函数
    'calculate_tm',
    'calculate_gc_content',
    'reverse_complement',
    'validate_genome_file',
    'validate_blast_database',
    'create_temp_fasta',
    'batch_process',
]
