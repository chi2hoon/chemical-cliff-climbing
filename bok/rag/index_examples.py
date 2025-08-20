
from sentence_transformers import SentenceTransformer
from qdrant_client import QdrantClient, models
from qdrant_client.http.models import PointStruct
import json
from pathlib import Path

def main():
    model = SentenceTransformer("sentence-transformers/all-MiniLM-L6-v2")
    base_dir = Path(__file__).resolve().parent.parent
    client = QdrantClient(path=str(base_dir / "rag" / "qdrant"))
    
    try:
        examples = []
        with open(base_dir / "rag" / "data" / "examples.jsonl", "r") as f:
            for line in f:
                ex = json.loads(line)
                prompt_like = f"{ex['input']['mol1']} vs {ex['input']['mol2']} 차이: {ex['input']['diff_type']}"
                embedding = model.encode(prompt_like)
                examples.append((embedding, ex['output']))

        # Qdrant 벡터 DB에 저장
        client.recreate_collection("sar_examples", vectors_config=models.VectorParams(size=384, distance=models.Distance.COSINE))
        client.upload_points(
            collection_name="sar_examples",
            points=[
                PointStruct(
                    id=i,
                    vector=e[0].tolist(),
                    payload={"output": e[1]}
                )
                for i, e in enumerate(examples)
            ],
            wait=True
        )
        print("Database has been successfully indexed.")
    finally:
        print("Closing Qdrant client...")
        client.close()

if __name__ == "__main__":
    main()
