"""路径配置管理"""

from pathlib import Path

# 默认路径
DEFAULT_PATHS = {
    'config_dir': Path.home() / '.config' / 'qtlprimer',
    'cache_dir': Path.home() / '.cache' / 'qtlprimer',
    'log_dir': Path.home() / '.local' / 'share' / 'qtlprimer' / 'logs',
}

def ensure_directories():
    """确保所有必要的目录存在"""
    for path in DEFAULT_PATHS.values():
        path.mkdir(parents=True, exist_ok=True)
