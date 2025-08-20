import types

from rag import retriever as r


class DummyModel:
    def encode(self, text):
        # Return a fixed-length vector resembling MiniLM (384 dims)
        return [0.0] * 384


class DummyHit:
    def __init__(self, output):
        self.payload = {"output": output}


class DummyClient:
    def search(self, collection_name, query_vector, limit):
        assert collection_name == "sar_examples"
        assert isinstance(query_vector, (list, tuple))
        assert limit == 3
        return [DummyHit("example-1"), DummyHit("example-2"), DummyHit("example-3")]


def test_retrieve_examples_monkeypatched(monkeypatch):
    # Monkeypatch the sentence transformer loader to avoid network/model download
    monkeypatch.setattr(r, "get_sentence_transformer", lambda: DummyModel())

    client = DummyClient()
    out = r.retrieve_examples(client, "c1ccccc1", "Cc1ccccc1", "lipophilicity", k=3)
    assert out == ["example-1", "example-2", "example-3"]

