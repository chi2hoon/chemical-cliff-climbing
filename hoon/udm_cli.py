#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import sys
import yaml

from src.udm.bronze_validator import BronzeValidator
from src.udm.silver import build_measurements_std
from src.udm.qc import qc_reports
from src.udm.silver_smiles import build_canonical_compounds
from src.udm.silver_ac import build_activity_cliffs
from src.udm.export import export_bronze_release


def load_cfg(cfg_path: str) -> dict:
	with open(cfg_path, "r", encoding="utf-8") as f:
		return yaml.safe_load(f)


def main() -> int:
	parser = argparse.ArgumentParser(description="UDM pipeline CLI")
	parser.add_argument("stage", choices=["bronze", "silver", "qc", "smiles", "ac", "ac-all", "export"], help="Pipeline stage to run")
	parser.add_argument("--config", default="configs/2017.yml", help="Path to YAML config")
	parser.add_argument("--root", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "data"), help="Data root directory (defaults to ./data)")
	args = parser.parse_args()

	cfg = load_cfg(args.config)
	root_dir = args.root

	if args.stage == "bronze":
		validator = BronzeValidator(cfg)
		results = validator.build_or_validate(root_dir)
		if results["errors"]:
			for e in results["errors"]:
				print(f"ERROR: {e}")
			return 1
		for p in results["fixed"]:
			print(f"Fixed deviations in: {p}")
		return 0

	if args.stage == "silver":
		out = build_measurements_std(root_dir, cfg)
		print(f"Wrote standardized measurements: {out}")
		return 0

	if args.stage == "qc":
		outs = qc_reports(root_dir, cfg)
		for k, v in outs.items():
			print(f"Generated {k}: {v}")
		return 0

	if args.stage == "smiles":
		out = build_canonical_compounds(root_dir, cfg)
		print(f"Wrote canonical compounds: {out}")
		return 0


	if args.stage == "ac":
		out = build_activity_cliffs(root_dir, [cfg])
		print(f"Wrote activity cliff pairs: {out}")
		return 0

	if args.stage == "ac-all":
		cfg_paths = [
			"configs/2017.yml",
			"configs/2018.yml",
			"configs/2020.yml",
			"configs/2021.yml",
		]
		cfgs = [load_cfg(p) for p in cfg_paths if os.path.exists(p)]
		out = build_activity_cliffs(root_dir, cfgs)
		print(f"Wrote activity cliff pairs: {out}")
		return 0

	if args.stage == "export":
		out = export_bronze_release(root_dir, cfg, mode="relative")
		print(f"Wrote bronze release: {out}")
		return 0

	return 0


if __name__ == "__main__":
	sys.exit(main())


