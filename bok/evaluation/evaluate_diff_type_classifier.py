import json
from pathlib import Path
from rdkit import Chem
from structure_diff import detect_diff_type
from sklearn.metrics import classification_report, confusion_matrix
from collections import Counter

def load_eval_data(path: str | None = None):
    data = []
    if path is None:
        base_dir = Path(__file__).resolve().parent
        path = base_dir / "diff_type_eval.jsonl"
    with open(path, "r") as f:
        for line in f:
            row = json.loads(line)
            row["mol1_rdkit"] = Chem.MolFromSmiles(row["mol1"])
            row["mol2_rdkit"] = Chem.MolFromSmiles(row["mol2"])
            data.append(row)
    return data

def evaluate_classifier():
    data = load_eval_data()
    y_true = []
    y_pred = []

    for row in data:
        pred = detect_diff_type(row["mol1_rdkit"], row["mol2_rdkit"])
        y_true.append(row["label"])
        y_pred.append(pred)

    print("📊 Classification Report")
    print(classification_report(y_true, y_pred, digits=3))

    print("\n🧭 Confusion Matrix")
    labels = sorted(list(set(y_true + y_pred)))
    cm = confusion_matrix(y_true, y_pred, labels=labels)
    print("Labels:", labels)
    print(cm)

    print("\n📌 Class Distribution")
    print("True Labels:", Counter(y_true))
    print("Predicted:", Counter(y_pred))

if __name__ == "__main__":
    evaluate_classifier()

