.PHONY: help install dev-install test lint format clean build publish docs

help:
	@echo "QTL Primer Toolkit - 开发命令"
	@echo ""
	@echo "可用命令:"
	@echo "  install       安装生产环境依赖"
	@echo "  dev-install   安装开发环境依赖"
	@echo "  test         运行测试套件"
	@echo "  lint         代码风格检查"
	@echo "  format       代码格式化"
	@echo "  clean        清理生成文件"
	@echo "  build        构建包"
	@echo "  publish      发布到PyPI"
	@echo "  docs         生成文档"
	@echo "  help         显示此帮助"

install:
	pip install -e .

dev-install:
	pip install -e ".[dev]"

test:
	pytest tests/ -v --cov=qtlprimer --cov-report=term-missing

lint:
	flake8 qtlprimer/ tests/ examples/
	mypy qtlprimer/

format:
	black qtlprimer/ tests/ examples/ setup.py
	isort qtlprimer/ tests/ examples/

clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info
	rm -rf .pytest_cache/
	rm -rf .coverage
	rm -rf htmlcov/
	rm -rf docs/_build/
	find . -type f -name "*.pyc" -delete
	find . -type d -name "__pycache__" -delete

build: clean
	python -m build

publish: build
	python -m twine upload dist/*

docs:
	cd docs && make html

run-example:
	python examples/basic_usage.py

check:
	qtlprimer check --all --variants examples/data/test_indels.vcf \
		--reference examples/data/test_genome.fasta \
		--blast-db examples/data/test_blastdb

version:
	@python -c "import qtlprimer; print(f'Version: {qtlprimer.__version__}')"

docker-build:
	docker build -t qtl-primer-toolkit:latest .

docker-run:
	docker run -v $(PWD)/data:/data qtl-primer-toolkit:latest \
		qtlprimer design \
		--variants /data/variants.vcf \
		--reference /data/genome.fasta \
		--blast-db /data/blastdb \
		--output /data/results
