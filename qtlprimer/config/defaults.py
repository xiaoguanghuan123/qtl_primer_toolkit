"""
默认配置文件 - 所有可配置参数集中管理
"""

# ============================================
# Primer3 设计参数
# ============================================
PRIMER3_DEFAULTS = {
    # 引物数量控制
    'PRIMER_NUM_RETURN': 5,
    'PRIMER_NUM_NS_ACCEPTED': 0,
    
    # 引物长度
    'PRIMER_MIN_SIZE': 18,
    'PRIMER_OPT_SIZE': 20,
    'PRIMER_MAX_SIZE': 24,
    
    # 熔解温度
    'PRIMER_MIN_TM': 50.0,
    'PRIMER_OPT_TM': 55.0,
    'PRIMER_MAX_TM': 60.0,
    'PRIMER_MAX_DIFF_TM': 2.0,
    
    # GC含量
    'PRIMER_MIN_GC': 40.0,
    'PRIMER_OPT_GC_PERCENT': 50.0,
    'PRIMER_MAX_GC': 60.0,
    'PRIMER_GC_CLAMP': 1,
    
    # 产物大小
    'PRODUCT_SIZE_MIN_OFFSET': 100,
    'PRODUCT_SIZE_MAX_OFFSET': 150,
    'PRIMER_PRODUCT_MIN': 80,
    'PRIMER_PRODUCT_MAX': 300,
    
    # 二级结构
    'PRIMER_MAX_SELF_ANY': 8.0,
    'PRIMER_MAX_SELF_END': 3.0,
    'PRIMER_MAX_HAIRPIN': 47.0,
    'PRIMER_MAX_POLY_X': 4,
    
    # 其他
    'PRIMER_EXPLAIN_FLAG': 1,
    'PRIMER_LIBERAL_BASE': 1,
    'PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS': 0,
}

# ============================================
# BLAST 验证参数
# ============================================
BLAST_DEFAULTS = {
    'TASK': 'blastn-short',
    'EVALUE_THRESHOLD': 10.0,
    'PERC_IDENTITY': 100.0,
    'WORD_SIZE': 7,
    'GAPOPEN': 5,
    'GAPEXTEND': 2,
    'PENALTY': -3,
    'REWARD': 2,
    'NUM_THREADS': 4,
    'MAX_TARGET_SEQS': 100,
    
    # 严格匹配设置
    'REQUIRE_EXACT_MATCH': True,
    'REQUIRE_FULL_LENGTH': True,
    'ALLOW_MISMATCH': 0,
    'ALLOW_GAPS': 0,
    
    # 批量处理
    'BATCH_SIZE': 50,
    'MAX_CONCURRENT': 4,
}

# ============================================
# 工具通用参数
# ============================================
TOOL_DEFAULTS = {
    # 序列提取
    'FLANKING_SIZE': 300,
    'MIN_TEMPLATE_LENGTH': 50,
    'MAX_TEMPLATE_LENGTH': 1000,
    
    # Indel过滤
    'MIN_INDEL_LENGTH': 1,
    'MAX_INDEL_LENGTH': 100,
    'ALLOWED_VARIANT_TYPES': ['INSERTION', 'DELETION'],
    
    # 质量过滤
    'MIN_TM_QUALITY': 50.0,
    'MAX_TM_QUALITY': 60.0,
    'MAX_TM_DIFF_QUALITY': 1.5,
    'MIN_GC_QUALITY': 40.0,
    'MAX_GC_QUALITY': 60.0,
    
    # 输出控制
    'OUTPUT_FORMATS': ['csv', 'excel', 'html'],
    'CREATE_FASTA': True,
    'CREATE_GEL_IMAGES': True,
    'TOP_N_RESULTS': 10,
    
    # 文件命名
    'RESULT_DIR_PATTERN': 'qtl_results_%Y%m%d_%H%M%S',
    'FILE_PREFIX': 'qtl_primers',
}

# ============================================
# 可视化参数
# ============================================
VISUALIZATION_DEFAULTS = {
    # 图形设置
    'FIG_SIZE': (12, 8),
    'DPI': 300,
    'FONT_SIZE': 10,
    'FONT_FAMILY': 'sans-serif',
    
    # 颜色方案
    'COLORS': {
        'band': '#2E86AB',
        'hetero_band': '#A23B72',
        'gel_background': '#F8F9FA',
        'well': '#6C757D',
        'marker': '#495057',
        'text': '#212529',
        'good': '#28a745',
        'warning': '#ffc107',
        'error': '#dc3545',
    },
    
    # 电泳参数
    'GEL_PARAMS': {
        'band_height': 0.6,
        'band_width_scale': 0.8,
        'hetero_spacing': 0.15,
        'blur_sigma': 0.7,
        'band_intensity': 0.9,
    },
    
    # 基因型示例
    'GENOTYPE_EXAMPLES': {
        'non_recombinant': ('AA', 'BB'),
        'recombinant': ('AA', 'Bb'),
        'double_heterozygous': ('Aa', 'Bb'),
    },
}

# ============================================
# 工作流参数
# ============================================
WORKFLOW_DEFAULTS = {
    'ENABLE_VALIDATION': True,
    'ENABLE_LOGGING': True,
    'SAVE_INTERMEDIATE': True,
    'PARALLEL_PROCESSING': True,
    'MAX_WORKERS': 8,
    
    'STEPS': {
        'parse_variants': True,
        'filter_indels': True,
        'design_primers': True,
        'check_specificity': True,
        'create_visualizations': True,
        'generate_reports': True,
    },
    
    'VALIDATION': {
        'check_reference': True,
        'check_blast_db': True,
        'validate_variants': True,
        'min_variants': 1,
        'max_variants': 1000,
    },
}

# ============================================
# 输入输出格式
# ============================================
IO_DEFAULTS = {
    'SUPPORTED_INPUT_FORMATS': ['vcf', 'vcf.gz', 'csv', 'tsv', 'txt'],
    'SUPPORTED_OUTPUT_FORMATS': ['csv', 'tsv', 'excel', 'html', 'json', 'fasta'],
    
    'CSV_FORMAT': {
        'delimiter': ',',
        'encoding': 'utf-8',
        'float_format': '%.3f',
    },
    
    'VCF_REQUIRED_COLUMNS': ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'],
    'CSV_REQUIRED_COLUMNS': ['CHR', 'POS', 'REF', 'ALT'],
    
    'OUTPUT_COLUMNS': [
        'PRIMER_PAIR_ID', 'CHR', 'POS', 'REF', 'ALT',
        'LEFT_PRIMER', 'RIGHT_PRIMER', 'PRODUCT_SIZE',
        'LEFT_TM', 'RIGHT_TM', 'LEFT_GC', 'RIGHT_GC',
        'SPECIFICITY_SCORE', 'SPECIFICITY_GRADE',
        'EXACT_HITS_F', 'EXACT_HITS_R',
        'QUALITY', 'RECOMMENDATION'
    ],
}

# ============================================
# 合并所有配置
# ============================================
DEFAULT_CONFIG = {
    'primer3': PRIMER3_DEFAULTS,
    'blast': BLAST_DEFAULTS,
    'tool': TOOL_DEFAULTS,
    'visualization': VISUALIZATION_DEFAULTS,
    'workflow': WORKFLOW_DEFAULTS,
    'io': IO_DEFAULTS,
}

def get_config(user_config=None):
    """
    获取配置，合并用户自定义配置
    
    Args:
        user_config: 用户配置字典
        
    Returns:
        合并后的配置字典
    """
    import copy
    
    # 深拷贝默认配置
    config = copy.deepcopy(DEFAULT_CONFIG)
    
    # 合并用户配置
    if user_config:
        for section in config.keys():
            if section in user_config:
                if isinstance(config[section], dict) and isinstance(user_config[section], dict):
                    config[section].update(user_config[section])
                else:
                    config[section] = user_config[section]
    
    return config

def validate_config(config):
    """
    验证配置有效性
    
    Args:
        config: 配置字典
        
    Returns:
        (is_valid, errors): 是否有效和错误列表
    """
    errors = []
    
    # 验证必要参数
    required_params = {
        'primer3': ['PRIMER_MIN_SIZE', 'PRIMER_MAX_SIZE'],
        'blast': ['EVALUE_THRESHOLD', 'PERC_IDENTITY'],
        'tool': ['FLANKING_SIZE', 'MIN_INDEL_LENGTH'],
    }
    
    for section, params in required_params.items():
        if section in config:
            for param in params:
                if param not in config[section]:
                    errors.append(f"Missing required parameter: {section}.{param}")
    
    # 验证数值范围
    if 'primer3' in config:
        p3 = config['primer3']
        if p3['PRIMER_MIN_SIZE'] > p3['PRIMER_MAX_SIZE']:
            errors.append("PRIMER_MIN_SIZE cannot be greater than PRIMER_MAX_SIZE")
        
        if p3['PRIMER_MIN_TM'] > p3['PRIMER_MAX_TM']:
            errors.append("PRIMER_MIN_TM cannot be greater than PRIMER_MAX_TM")
    
    return len(errors) == 0, errors
