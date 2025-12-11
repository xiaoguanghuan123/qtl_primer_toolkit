"""变异解析器测试"""

import pytest
import pandas as pd
from pathlib import Path
from qtlprimer.core.variant_parser import VariantParser

def test_variant_parser_creation():
    """测试解析器创建"""
    parser = VariantParser()
    assert parser is not None

def test_parse_csv():
    """测试CSV解析"""
    parser = VariantParser()
    
    # 创建测试CSV数据
    test_data = """CHR,POS,REF,ALT
1,1000,A,AT
1,2000,AT,A
2,3000,C,CT"""
    
    # 保存为临时文件
    test_file = Path("test_variants.csv")
    test_file.write_text(test_data)
    
    try:
        df = parser.parse(test_file)
        assert len(df) == 3
        assert 'CHR' in df.columns
        assert 'POS' in df.columns
    finally:
        test_file.unlink()

def test_filter_indels():
    """测试Indel过滤"""
    parser = VariantParser()
    
    # 创建测试数据
    data = {
        'CHR': ['1', '1', '2'],
        'POS': [1000, 2000, 3000],
        'REF': ['A', 'AT', 'C'],
        'ALT': ['AT', 'A', 'CT']
    }
    df = pd.DataFrame(data)
    
    filtered = parser.filter_indels(df)
    assert len(filtered) == 2  # 两个indel
