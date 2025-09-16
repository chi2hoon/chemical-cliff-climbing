# Chemical Cilff Climbing (ccc -> ㅋㅋㅋ)

이 애플리케이션은 Streamlit을 사용하여 화학적 화합물과 관련된 가설을 생성, 수정 및 평가하는 웹 기반 도구입니다.

## 사전 준비 사항

- **Python:** `3.13.5` 버전에서 개발 및 테스트되었습니다.
- **OpenAI API 키:** `openai` 라이브러리 사용을 위해 필요합니다.

레포 루트 디렉토리에 `openAI_key.txt` 파일을 생성하고, 파일 내에 자신의 OpenAI API 키를 붙여넣으세요. 예:
```
sk-xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
```

## 실행 방법

환경 구성 및 실행 방법에는 두 가지 옵션이 있습니다.

### 옵션 1: `run.sh` 스크립트 사용 (venv 기반, 권장)

이 디렉토리에 포함된 `run.sh` 스크립트는 `venv` 가상 환경 생성, 필요 라이브러리 설치, 애플리케이션 실행까지 모든 과정을 자동으로 처리합니다.

레포 루트가 아닌, 반드시 `base/` 디렉토리에서 다음을 실행하세요:
```bash
./run.sh
```

### 옵션 2: Conda 환경에서 수동 실행

Conda를 사용하여 환경을 직접 구성하고 실행할 수 있습니다.

1.  **Conda 가상 환경 생성 및 활성화:**
    `ccc-env`라는 이름으로 Python 3.13 버전의 Conda 환경을 생성합니다.
    ```bash
    conda create --name ccc-env python=3.13
    conda activate ccc-env
    ```

2.  **필요한 라이브러리 설치:**
    `requirements.txt` 파일을 사용하여 필요한 라이브러리를 설치합니다.
    ```bash
    pip install -r requirements.txt
    ```

3.  **애플리케이션 실행:**
    Streamlit을 사용하여 앱을 직접 실행합니다.
    ```bash
    streamlit run app.py
    ```

---

두 방법 중 하나로 실행에 성공하면, 웹 브라우저에서 애플리케이션에 접속할 수 있는 로컬 주소(예: `http://localhost:8501`)가 터미널에 표시됩니다.
