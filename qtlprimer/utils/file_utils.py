"""文件工具函数"""

import gzip
import json
import yaml
from pathlib import Path
from typing import Union, Any

def read_file(file_path: Union[str, Path]) -> str:
    """读取文件，支持gzip压缩"""
    file_path = Path(file_path)
    
    if file_path.suffix == '.gz':
        with gzip.open(file_path, 'rt') as f:
            return f.read()
    else:
        with open(file_path, 'r') as f:
            return f.read()

def write_file(content: str, file_path: Union[str, Path]):
    """写入文件"""
    file_path = Path(file_path)
    file_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(file_path, 'w') as f:
        f.write(content)

def load_config(config_file: Union[str, Path]) -> dict:
    """加载配置文件（支持JSON、YAML）"""
    config_file = Path(config_file)
    content = read_file(config_file)
    
    if config_file.suffix in ['.yaml', '.yml']:
        return yaml.safe_load(content)
    elif config_file.suffix == '.json':
        return json.loads(content)
    else:
        raise ValueError(f"不支持的配置文件格式: {config_file.suffix}")
