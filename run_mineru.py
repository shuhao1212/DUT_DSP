"""使用 MinerU 提取实验指导书PDF的文字内容"""
import os, sys, json
os.environ["MINERU_MODEL_SOURCE"] = "modelscope"
from pathlib import Path

pdf = Path(r"项目实验指导书v7new3.pdf")
out = Path(r"output_mineru")
out.mkdir(exist_ok=True)

# 设置环境变量，使用 pipeline 后端（纯CPU）
os.environ["MINERU_DEVICE"] = "cpu"

from mineru.api import MinerU
m = MinerU()
result = m.parse(pdf_path=str(pdf), output_dir=str(out), method="ocr", backend="pipeline")
print("Result:", result)
print()

# 读取输出
for f in out.glob("*"):
    print(f.name, f.stat().st_size)
    if f.suffix in ('.txt', '.md', '.json'):
        print(f.read_text(encoding='utf-8', errors='replace')[:3000])
        print('---')
