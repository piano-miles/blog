import re


def fix_inline_latex(file_path, overwrite=True):
    with open(file_path, "r", encoding="utf-8") as f:
        content = f.read()

    def replacer(match):
        inner = match.group(1)
        fixed_inner = inner.replace("*", "_")
        return f"\\({fixed_inner}\\)"

    new_content = re.sub("\\\\\\((.*?)\\\\\\)", replacer, content, flags=re.DOTALL)
    if overwrite:
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(new_content)
    else:
        print(new_content)
    print("Replacement complete.")


fix_inline_latex(
    "/Users/miles/Documents/GitHub/blog/content/posts/infinite_tetration/index.md",
    overwrite=True,
)
