"""Current PMC cloud metadata must establish a published, nonretracted source."""
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from recover_pmc_sources import named_media


def _metadata():
    return {"pmcid": "PMC123", "pmid": 123, "is_manuscript": False, "is_retracted": False,
            "media_urls": ["s3://pmc-oa-opendata/PMC123.2/supplement.xlsx?md5=" + "0" * 32]}


def test_cloud_s3_links_use_the_same_public_bucket_object_and_checksum():
    record = _metadata()
    urls = named_media(record, "PMC123", "123", "supplement.xlsx")
    assert urls == ["https://pmc-oa-opendata.s3.amazonaws.com/PMC123.2/supplement.xlsx?md5=" + "0" * 32]
    assert named_media(record, "PMC123", "123", "different.xlsx") == []
    record["is_manuscript"] = True
    assert named_media(record, "PMC123", "123", "supplement.xlsx") == []


@pytest.mark.parametrize("changes", [
    {"pmid": 124}, {"pmcid": "PMC124"}, {"is_retracted": True},
    {"media_urls": ["https://example.org/supplement.xlsx?md5=" + "0" * 32]},
    {"media_urls": ["s3://pmc-oa-opendata/PMC124.1/supplement.xlsx?md5=" + "0" * 32]},
    {"media_urls": ["s3://pmc-oa-opendata/PMC123.1/supplement.xlsx"]},
])
def test_no_binding_when_article_identity_or_published_file_identity_is_uncertain(changes):
    with pytest.raises(ValueError):
        named_media({**_metadata(), **changes}, "PMC123", "123", "supplement.xlsx")
