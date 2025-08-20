import streamlit as st
from sentence_transformers import SentenceTransformer
from qdrant_client import QdrantClient
import numpy as np

@st.cache_resource
def get_sentence_transformer():
    """Loads the sentence transformer model and caches it."""
    return SentenceTransformer("sentence-transformers/all-MiniLM-L6-v2")

def retrieve_examples(client: QdrantClient, mol1: str, mol2: str, diff_type: str, k=3):
    """Retrieves examples from the Qdrant database using a provided client."""
    model = get_sentence_transformer()
    query = f"{mol1} vs {mol2} 차이: {diff_type}"
    query_vec = model.encode(query)

    hits = client.search(
        collection_name="sar_examples",
        query_vector=query_vec,
        limit=k
    )
    return [hit.payload["output"] for hit in hits]
