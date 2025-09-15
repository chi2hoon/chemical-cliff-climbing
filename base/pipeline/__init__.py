"""설정 주도 파이프라인 패키지(pipeline).

모듈 구성:
- bronze_dispatch: 연도별 어댑터 레지스트리 및 디스패처
- adapters: 연도별 브론즈 어댑터 진입점(run)
- parsers, transforms, validate: 후속 단계에서 도입

주의: 타입힌트 미사용. 간단한 docstring으로 인자/반환 요약.
"""

