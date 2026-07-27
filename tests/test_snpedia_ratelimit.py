"""Regression guard: the harvesters must honour their own rate_limit argument.

harvest_genosets built a _RateLimiter from the caller's rate_limit and then never
used it, so every request silently fell back to the module-level 2.0/s default
and the argument was a lie. harvest had the same gap for its per-page fetches:
it threaded the limiter into target enumeration only. flake8 saw it as an unused
local (F841); it was a real defect.

Offline by construction. Nothing here touches the network or the real cache.
"""

import ast
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from backend import snpedia
from backend.snpedia import DEFAULT_RATE_LIMIT, NOTICE, _RateLimiter

SLOW_RATE = 0.25          # one request every 4 seconds
SLOW_INTERVAL = 1.0 / SLOW_RATE
DEFAULT_INTERVAL = 1.0 / DEFAULT_RATE_LIMIT


@pytest.fixture()
def cache(tmp_path, monkeypatch):
    """Redirect the cache into tmp_path so the real one is never touched."""
    monkeypatch.setattr(snpedia, "CACHE_DIR", tmp_path)
    monkeypatch.setattr(snpedia, "CACHE_PATH", tmp_path / "snpedia_cache.db")
    monkeypatch.setattr(snpedia, "_cache_dir", lambda: tmp_path)
    monkeypatch.setattr(snpedia, "cache_path", lambda: tmp_path / "snpedia_cache.db")
    return tmp_path


@pytest.fixture()
def recorder(monkeypatch):
    """Capture the limiter handed to every page fetch, without any HTTP."""
    seen: list = []

    def fake_subject(subject, session=None, limiter=None):
        seen.append(("subject", subject, limiter))
        return {}

    def fake_wikitext(title, session=None, limiter=None):
        seen.append(("wikitext", title, limiter))
        return ""

    monkeypatch.setattr(snpedia, "fetch_subject", fake_subject)
    monkeypatch.setattr(snpedia, "fetch_wikitext", fake_wikitext)
    monkeypatch.setattr(snpedia, "make_session", lambda: None)
    return seen


class TestTheFetchersAcceptAndForwardALimiter:
    def test_fetch_subject_forwards_its_limiter(self, monkeypatch):
        captured = {}

        def fake_api_get(params, session=None, limiter=None):
            captured["limiter"] = limiter
            return {}

        monkeypatch.setattr(snpedia, "_api_get", fake_api_get)
        sentinel = _RateLimiter(SLOW_RATE)
        snpedia.fetch_subject("Rs1815739", limiter=sentinel)
        assert captured["limiter"] is sentinel

    def test_fetch_wikitext_forwards_its_limiter(self, monkeypatch):
        captured = {}

        def fake_request(url, params=None, session=None, limiter=None):
            captured["limiter"] = limiter
            return ""

        monkeypatch.setattr(snpedia, "_request", fake_request)
        sentinel = _RateLimiter(SLOW_RATE)
        snpedia.fetch_wikitext("Rs1815739", limiter=sentinel)
        assert captured["limiter"] is sentinel

    def test_omitting_the_limiter_still_falls_back_to_the_default(self, monkeypatch):
        captured = {}

        def fake_api_get(params, session=None, limiter=None):
            captured["limiter"] = limiter
            return {}

        monkeypatch.setattr(snpedia, "_api_get", fake_api_get)
        snpedia.fetch_subject("Rs1815739")
        assert captured["limiter"] is None


class TestHarvestGenosetsHonoursRateLimit:
    def test_the_caller_rate_reaches_every_fetch(self, cache, recorder, monkeypatch):
        monkeypatch.setattr(snpedia, "enumerate_genosets",
                            lambda session=None, progress_cb=None: ["gs100", "gs101"])
        snpedia.harvest_genosets(accept_license=True, rate_limit=SLOW_RATE)
        assert recorder, "no page was fetched, so the test proved nothing"
        for kind, target, limiter in recorder:
            assert limiter is not None, "%s %s fetched with no limiter" % (kind, target)
            assert limiter is not snpedia._DEFAULT_LIMITER
            assert limiter.min_interval == pytest.approx(SLOW_INTERVAL)

    def test_one_limiter_is_shared_across_the_whole_run(self, cache, recorder, monkeypatch):
        monkeypatch.setattr(snpedia, "enumerate_genosets",
                            lambda session=None, progress_cb=None: ["gs100", "gs101"])
        snpedia.harvest_genosets(accept_license=True, rate_limit=SLOW_RATE)
        limiters = {id(entry[2]) for entry in recorder}
        assert len(limiters) == 1, "a per-request limiter does not rate limit anything"

    def test_the_default_rate_is_still_the_default(self, cache, recorder, monkeypatch):
        monkeypatch.setattr(snpedia, "enumerate_genosets",
                            lambda session=None, progress_cb=None: ["gs100"])
        snpedia.harvest_genosets(accept_license=True)
        assert recorder
        for _kind, _target, limiter in recorder:
            assert limiter.min_interval == pytest.approx(DEFAULT_INTERVAL)

    def test_the_licence_gate_still_fires_first(self, cache, recorder):
        with pytest.raises(PermissionError) as exc:
            snpedia.harvest_genosets(rate_limit=SLOW_RATE)
        assert NOTICE in str(exc.value)
        assert recorder == [], "a request was made before the licence was accepted"


class TestHarvestHonoursRateLimit:
    def test_the_caller_rate_reaches_every_snp_fetch(self, cache, recorder, monkeypatch):
        monkeypatch.setattr(snpedia, "_resolve_targets",
                            lambda scope, rsids, session, limiter, progress_cb: ["rs1815739"])
        snpedia.harvest(accept_license=True, rate_limit=SLOW_RATE)
        assert recorder, "no page was fetched, so the test proved nothing"
        for kind, target, limiter in recorder:
            assert limiter is not None, "%s %s fetched with no limiter" % (kind, target)
            assert limiter.min_interval == pytest.approx(SLOW_INTERVAL)

    def test_target_enumeration_gets_the_same_limiter_as_the_fetches(
            self, cache, recorder, monkeypatch):
        seen = {}

        def fake_resolve(scope, rsids, session, limiter, progress_cb):
            seen["limiter"] = limiter
            return ["rs1815739"]

        monkeypatch.setattr(snpedia, "_resolve_targets", fake_resolve)
        snpedia.harvest(accept_license=True, rate_limit=SLOW_RATE)
        assert recorder
        assert all(entry[2] is seen["limiter"] for entry in recorder)


class TestNoFetchCallSiteForgetsTheLimiter:
    """Source-level guard. A new call site added without limiter= fails here."""

    @staticmethod
    def _fetch_calls(func_name: str):
        src = Path(snpedia.__file__).read_text(encoding="utf-8")
        tree = ast.parse(src)
        target = next(n for n in ast.walk(tree)
                      if isinstance(n, ast.FunctionDef) and n.name == func_name)
        calls = []
        for node in ast.walk(target):
            if not isinstance(node, ast.Call):
                continue
            fn = node.func
            name = fn.id if isinstance(fn, ast.Name) else getattr(fn, "attr", "")
            if name in ("fetch_subject", "fetch_wikitext"):
                kwargs = {kw.arg for kw in node.keywords}
                calls.append((name, kwargs))
        return calls

    @pytest.mark.parametrize("func_name", ["harvest", "harvest_genosets"])
    def test_every_fetch_passes_a_limiter(self, func_name):
        calls = self._fetch_calls(func_name)
        assert calls, "no fetch call found in %s, the guard would pass vacuously" % func_name
        for name, kwargs in calls:
            assert "limiter" in kwargs, (
                "%s calls %s without limiter=, so rate_limit is silently ignored"
                % (func_name, name))

    @pytest.mark.parametrize("func_name", ["harvest", "harvest_genosets"])
    def test_the_limiter_local_is_actually_consumed(self, func_name):
        """The exact defect flake8 flagged: a limiter built and never used."""
        src = Path(snpedia.__file__).read_text(encoding="utf-8")
        tree = ast.parse(src)
        target = next(n for n in ast.walk(tree)
                      if isinstance(n, ast.FunctionDef) and n.name == func_name)
        assigned = any(
            isinstance(n, ast.Assign)
            and any(isinstance(t, ast.Name) and t.id == "limiter" for t in n.targets)
            for n in ast.walk(target))
        loaded = any(
            isinstance(n, ast.Name) and n.id == "limiter" and isinstance(n.ctx, ast.Load)
            for n in ast.walk(target))
        assert assigned, "%s no longer builds a limiter" % func_name
        assert loaded, "%s builds a limiter and never reads it" % func_name
