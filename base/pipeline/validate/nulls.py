import os


def scan_workspace(root_dir):
    """Args: root_dir(str) -> list[str]

    워크스페이스 전체에서 널바이트가 포함된 파일 경로를 스캔한다.
    """
    bad = []
    for base, _dirs, files in os.walk(root_dir):
        for f in files:
            p = os.path.join(base, f)
            try:
                with open(p, 'rb') as fh:
                    b = fh.read()
                if b"\x00" in b:
                    bad.append(p)
            except Exception:
                continue
    return bad

