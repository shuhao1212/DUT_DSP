import fitz
doc = fitz.open(r'项目实验指导书v7new3.pdf')
print(f'总页数: {doc.page_count}')
for i in range(doc.page_count):
    page = doc[i]
    text = page.get_text()
    if text.strip():
        print(f'\n=== 第 {i+1} 页 ===')
        print(text[:3000])
