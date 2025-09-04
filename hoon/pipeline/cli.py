#!/usr/bin/env python3

import argparse
import os
import sys

from . import bronze_dispatch
from .silver import build_silver


def _default_yaml_for_year(year):
    """Args: year(str|int) -> str

    기본 스키마 경로를 반환한다.
    """
    y = str(year)
    return os.path.join("schemas", "silver", f"{y}.yaml")


def cmd_bronze(args):
    """Args: args(Namespace) -> int

    bronze 단계 실행. 어댑터 레지스트리를 통해 연도별 런을 호출한다.
    """
    year = str(args.year)
    yaml_path = args.cfg or _default_yaml_for_year(year)
    out_dir = args.out or os.path.join("data", "bronze", year)
    result = bronze_dispatch.build_bronze_from_raw(year, yaml_path, out_dir)
    for k in sorted(result.keys()):
        print(f"{k}: {result[k]}")
    return 0


def main():
    """Args: None -> int

    파이프라인 CLI 진입점.
    """
    parser = argparse.ArgumentParser(description="Medallion pipeline CLI")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_bronze = sub.add_parser("bronze", help="원본→브론즈(시트 슬라이스/매트릭스 롱)")
    p_bronze.add_argument("--year", required=True, help="데이터 연도")
    p_bronze.add_argument("--cfg", default=None, help="schemas/silver/{year}.yaml 경로")
    p_bronze.add_argument("--out", default=None, help="출력 디렉터리(기본: data/ingest/{year})")
    p_bronze.set_defaults(func=cmd_bronze)

    p_silver = sub.add_parser("silver", help="정규화/표준화 단계")
    p_silver.add_argument("--year", required=True)
    p_silver.add_argument("--cfg", default=None, help="schemas/silver/{year}.yaml 경로")
    def _run_silver(args):
        year = str(args.year)
        yaml_path = args.cfg or _default_yaml_for_year(year)
        out = build_silver(year, yaml_path)
        for k in sorted(out.keys()):
            print(f"{k}: {out[k]}")
        return 0
    p_silver.set_defaults(func=_run_silver)

    p_gold = sub.add_parser("gold", help="골드 통합 단계")
    p_gold.add_argument("--years", nargs="+", help="여러 연도")
    def _run_gold(args):
        from .gold import build_gold
        years = [str(y) for y in args.years]
        results = build_gold(years)
        for y in sorted(results.keys()):
            paths = results[y]
            print(f"{y}/assay_readings: {paths['assay_readings']}")
            print(f"{y}/compounds: {paths['compounds']}")
        return 0
    p_gold.set_defaults(func=_run_gold)

    p_val = sub.add_parser("validate", help="스키마/환경 검증")
    p_val.add_argument("--stage", choices=["silver", "gold", "nulls"], required=True)
    p_val.add_argument("--year", default=None)
    def _run_validate(args):
        if args.stage == "nulls":
            from .validate.nulls import scan_workspace
            bad = scan_workspace(".")
            if bad:
                print("NULLBYTE_FOUND:\n" + "\n".join(bad))
                return 2
            print("NO_NULLBYTES")
            return 0
        raise SystemExit("추가 검증 훅은 후속 커밋에서 확장")
    p_val.set_defaults(func=_run_validate)

    args = parser.parse_args()
    try:
        return int(args.func(args) or 0)
    except SystemExit as e:
        # 상위로 그대로 전파
        raise
    except Exception as e:
        print(f"ERROR: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
