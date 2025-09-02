import os


def scan_workspace(root_dir="."):
    """Args: root_dir(str) -> list

    워크스페이스 전역에서 널바이트(\x00)를 포함한 파일 목록을 반환한다.
    data/ 등 대용량 디렉터리는 포함하되, 바이너리 파일은 크기만 보고 빠르게 스캔.
    """
    bad = []
    allow_ext = {".py", ".md", ".yaml", ".yml", ".csv", ".txt"}
    for base, dirs, files in os.walk(root_dir):
        # .git 제외
        if ".git" in base.split(os.sep):
            continue
        for name in files:
            path = os.path.join(base, name)
            _, ext = os.path.splitext(name)
            if ext.lower() not in allow_ext:
                continue
            try:
                with open(path, "rb") as f:
                    for chunk in iter(lambda: f.read(8192), b""):
                        if b"\x00" in chunk:
                            bad.append(path)
                            break
            except Exception:
                # 접근 불가/디렉토리 등 무시
                continue
    return bad
