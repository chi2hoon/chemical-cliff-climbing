import os
import sys
import pandas as pd


def main() -> int:
	root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "raw")
	files = [
		"2018_WO2018172250_raw.xlsx",
		"2020_WO2020132269_raw.xlsx",
		"2021_WO2021163344_raw.xlsx",
	]
	for fn in files:
		path = os.path.abspath(os.path.join(root, fn))
		print(f"== {fn} ==")
		if not os.path.exists(path):
			print("MISSING:", path)
			continue
		try:
			xl = pd.ExcelFile(path)
			print("sheets:", xl.sheet_names)
			for sn in xl.sheet_names[:3]:
				print(f"-- {sn} head --")
				df = xl.parse(sn, nrows=5, header=None)
				print(df.fillna("").astype(str).to_string(index=False))
		except Exception as e:
			print("ERR", e)
	print("done")
	return 0


if __name__ == "__main__":
	sys.exit(main())


