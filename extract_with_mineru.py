"""使用 MinerU 提取项目实验指导书的文字内容"""
import sys
import os
from pathlib import Path

# 设置不要自动下载模型，只用OCR模式
os.environ["MINERU_MODEL_SOURCE"] = "modelscope"

pdf_path = Path(r"项目实验指导书v7new3.pdf")
output_dir = Path(r"output_mineru")

print(f"PDF: {pdf_path}")
print(f"PDF exists: {pdf_path.exists()}")
print(f"PDF size: {pdf_path.stat().st_size / 1024 / 1024:.1f} MB")
print(f"Python: {sys.version}")

# 尝试使用 mineru 命令行工具
import subprocess
result = subprocess.run(
    [sys.executable, "-m", "mineru.cli", "-p", str(pdf_path), "-o", str(output_dir), "-m", "ocr", "-b", "pipeline"],
    capture_output=True, text=True, timeout=600
)
print("STDOUT:", result.stdout[-2000:] if len(result.stdout) > 2000 else result.stdout)
print("STDERR:", result.stderr[-2000:] if len(result.stderr) > 2000 else result.stderr)
print("Return code:", result.returncode)
