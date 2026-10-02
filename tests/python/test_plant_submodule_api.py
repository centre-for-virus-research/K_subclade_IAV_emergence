"""Contract tests for the vendored PLANT package (``external/PLANT/src/plant``).

This is a git submodule pinned to TheSatoLab/PLANT. Both pipeline scripts import
five names from it and would break silently or loudly if upstream changed their
signatures, so the contract is pinned here rather than assumed.

The tokenisation tests also record a consequential upstream quirk: ``MAX_LENGTH``
is 329, the same as the HA1 length, but the ESM tokenizer adds ``<cls>``/``<eos>``,
so the final two residues of every sequence are truncated before the model sees
them. Upstream *training* used the identical setting, so reproducing it is
correct; the test exists to stop someone "fixing" it and silently changing every
published coordinate.
"""

from __future__ import annotations

import pytest

plant = pytest.importorskip("plant", reason="PLANT submodule not initialised or torch missing")
torch = pytest.importorskip("torch")

pytestmark = pytest.mark.slow

MODEL_NAME = "facebook/esm2_t33_650M_UR50D"


@pytest.fixture(scope="module")
def tokenizer():
    transformers = pytest.importorskip("transformers")
    try:
        return transformers.AutoTokenizer.from_pretrained(MODEL_NAME)
    except Exception as exc:  # noqa: BLE001 - offline runs should skip, not fail
        pytest.skip(f"ESM-2 tokenizer unavailable: {exc}")


class TestPublicApi:
    def test_exports_exactly_the_documented_names(self):
        assert set(plant.__all__) == {
            "TextDataset",
            "tokenize_sequences",
            "semanticESM",
            "set_encoders",
            "embed_sequences",
        }

    @pytest.mark.parametrize(
        "name", ["TextDataset", "tokenize_sequences", "semanticESM", "set_encoders", "embed_sequences"]
    )
    def test_every_exported_name_is_importable(self, name):
        assert getattr(plant, name) is not None

    def test_semantic_esm_is_a_pretrained_model(self):
        from transformers import PreTrainedModel

        assert issubclass(plant.semanticESM, PreTrainedModel)

    def test_set_encoders_takes_four_encoders(self):
        import inspect

        params = inspect.signature(plant.set_encoders).parameters
        assert list(params) == ["ohe_virus", "ohe_ref", "ohe_vp", "ohe_rp"]

    def test_embed_sequences_signature(self):
        import inspect

        params = inspect.signature(plant.embed_sequences).parameters
        assert list(params) == ["model", "dataloader", "use_fp16"]
        assert params["use_fp16"].default is True


class TestTokenizeSequences:
    def test_returns_input_ids_and_attention_mask(self, tokenizer, reference_seq):
        encoded = plant.tokenize_sequences([reference_seq], tokenizer, 329)
        assert "input_ids" in encoded and "attention_mask" in encoded

    def test_returns_torch_tensors(self, tokenizer, reference_seq):
        encoded = plant.tokenize_sequences([reference_seq], tokenizer, 329)
        assert isinstance(encoded["input_ids"], torch.Tensor)
        assert isinstance(encoded["attention_mask"], torch.Tensor)

    def test_pads_every_row_to_max_length(self, tokenizer, reference_seq):
        encoded = plant.tokenize_sequences([reference_seq, reference_seq[:50]], tokenizer, 329)
        assert encoded["input_ids"].shape == (2, 329)

    def test_batch_dimension_matches_input_count(self, tokenizer, reference_seq):
        for n in [1, 3, 10]:
            encoded = plant.tokenize_sequences([reference_seq] * n, tokenizer, 329)
            assert encoded["input_ids"].shape[0] == n

    def test_short_sequence_is_padded_not_truncated(self, tokenizer, reference_seq):
        encoded = plant.tokenize_sequences([reference_seq[:50]], tokenizer, 329)
        mask = encoded["attention_mask"][0]
        # 50 residues + <cls> + <eos>
        assert int(mask.sum()) == 52

    def test_full_length_input_saturates_the_window(self, tokenizer, reference_seq):
        encoded = plant.tokenize_sequences([reference_seq], tokenizer, 329)
        assert int(encoded["attention_mask"][0].sum()) == 329

    def test_special_tokens_cost_two_residues(self, tokenizer, reference_seq):
        """Characterisation of the upstream MAX_LENGTH quirk.

        A 329 aa sequence needs 331 token slots; with max_length=329 the last
        two residues are dropped. Training used the same setting, so this is
        reproduced deliberately.
        """
        encoded = plant.tokenize_sequences([reference_seq], tokenizer, 329)
        ids = encoded["input_ids"][0].tolist()

        decoded = tokenizer.decode(ids, skip_special_tokens=True).replace(" ", "")
        assert len(decoded) == 327, f"expected 327 residues to survive, got {len(decoded)}"
        assert reference_seq.startswith(decoded)
        assert decoded == reference_seq[:327]

    def test_raising_max_length_would_admit_the_whole_sequence(self, tokenizer, reference_seq):
        """Shows the quirk is a max_length choice, not a tokenizer limitation."""
        encoded = plant.tokenize_sequences([reference_seq], tokenizer, 331)
        decoded = tokenizer.decode(encoded["input_ids"][0], skip_special_tokens=True).replace(" ", "")
        assert decoded == reference_seq

    def test_is_deterministic(self, tokenizer, reference_seq):
        a = plant.tokenize_sequences([reference_seq], tokenizer, 329)["input_ids"]
        b = plant.tokenize_sequences([reference_seq], tokenizer, 329)["input_ids"]
        assert torch.equal(a, b)

    def test_distinct_sequences_tokenise_differently(self, tokenizer, reference_seq):
        mutant = reference_seq[:100] + ("W" if reference_seq[100] != "W" else "Y") + reference_seq[101:]
        a = plant.tokenize_sequences([reference_seq], tokenizer, 329)["input_ids"]
        b = plant.tokenize_sequences([mutant], tokenizer, 329)["input_ids"]
        assert not torch.equal(a, b)


class TestTextDataset:
    @staticmethod
    def _encoded(n: int, length: int = 8) -> dict:
        return {
            "input_ids": torch.arange(n * length).reshape(n, length),
            "attention_mask": torch.ones(n, length, dtype=torch.long),
        }

    def test_length_matches_row_count(self):
        assert len(plant.TextDataset(self._encoded(5))) == 5

    def test_item_exposes_the_keys_embed_sequences_reads(self):
        item = plant.TextDataset(self._encoded(2))[0]
        assert "input_ids_virus" in item and "attention_mask_virus" in item

    def test_item_preserves_row_content(self):
        encoded = self._encoded(3)
        item = plant.TextDataset(encoded)[1]
        assert torch.equal(item["input_ids_virus"], encoded["input_ids"][1])

    def test_optional_fields_default_without_a_reference(self):
        item = plant.TextDataset(self._encoded(1))[0]
        assert item["weight"].item() == pytest.approx(1.0)
        assert item["labels"].item() == pytest.approx(-10.0)
        assert "input_ids_reference" not in item

    def test_reference_fields_appear_when_supplied(self):
        encoded = self._encoded(2)
        dataset = plant.TextDataset(encoded, encodes_reference=self._encoded(2))
        assert "input_ids_reference" in dataset[0]

    def test_mismatched_reference_length_is_rejected(self):
        with pytest.raises(AssertionError):
            plant.TextDataset(self._encoded(3), encodes_reference=self._encoded(2))

    def test_works_with_a_dataloader_in_order(self):
        from torch.utils.data import DataLoader

        encoded = self._encoded(6)
        loader = DataLoader(plant.TextDataset(encoded), batch_size=2, shuffle=False)
        rows = torch.cat([batch["input_ids_virus"] for batch in loader])
        assert torch.equal(rows, encoded["input_ids"])

    def test_unique_combination_indices_group_by_virus(self):
        dataset = plant.TextDataset(self._encoded(4), virus=[0, 1, 0, 1])
        groups = dataset.get_unique_combinations_indices()
        assert sorted(groups.values()) == [[0, 2], [1, 3]]


class TestEmbedSequences:
    """``embed_sequences`` is exercised with a stub model: no checkpoint needed."""

    class _StubOutput:
        def __init__(self, tensor):
            self.hidden_state_virus = tensor

    class _StubModel(torch.nn.Module):
        def __init__(self, latent_dim=3):
            super().__init__()
            self.latent_dim = latent_dim
            self.marker = torch.nn.Parameter(torch.zeros(1))

        def forward(self, input_ids_virus, attention_mask_virus):
            n = input_ids_virus.shape[0]
            base = input_ids_virus[:, :1].float()
            tensor = base.repeat(1, self.latent_dim) + torch.arange(self.latent_dim).float()
            return TestEmbedSequences._StubOutput(tensor.reshape(n, self.latent_dim))

    def _loader(self, n: int, batch_size: int):
        from torch.utils.data import DataLoader

        encoded = {
            "input_ids": torch.arange(n).reshape(n, 1).repeat(1, 4),
            "attention_mask": torch.ones(n, 4, dtype=torch.long),
        }
        return DataLoader(plant.TextDataset(encoded), batch_size=batch_size, shuffle=False)

    def test_returns_one_row_per_sequence(self):
        out = plant.embed_sequences(self._StubModel(), self._loader(7, 3), use_fp16=False)
        assert out.shape == (7, 3)

    def test_returns_float_numpy(self):
        import numpy as np

        out = plant.embed_sequences(self._StubModel(), self._loader(4, 2), use_fp16=False)
        assert isinstance(out, np.ndarray)
        assert out.dtype == np.float32

    def test_row_order_matches_input_order(self):
        out = plant.embed_sequences(self._StubModel(), self._loader(6, 2), use_fp16=False)
        assert out[:, 0].tolist() == [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]

    @pytest.mark.parametrize("batch_size", [1, 2, 3, 5, 64])
    def test_result_is_independent_of_batch_size(self, batch_size):
        import numpy as np

        out = plant.embed_sequences(self._StubModel(), self._loader(5, batch_size), use_fp16=False)
        np.testing.assert_allclose(out[:, 0], np.arange(5, dtype=np.float32))

    def test_leaves_the_model_in_eval_mode(self):
        model = self._StubModel()
        model.train()
        plant.embed_sequences(model, self._loader(2, 2), use_fp16=False)
        assert not model.training

    def test_does_not_accumulate_gradients(self):
        model = self._StubModel()
        plant.embed_sequences(model, self._loader(2, 2), use_fp16=False)
        assert model.marker.grad is None
