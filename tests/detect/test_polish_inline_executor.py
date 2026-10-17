"""One polish worker must not start a thread: Pyodide (the browser build) cannot, and a
ThreadPoolExecutor starts one on its first submit."""
import threading

import pytest

from MATPredict.detect import pipeline


@pytest.fixture(autouse=True)
def _fresh_pool(monkeypatch):
    monkeypatch.setattr(pipeline, "_POLISH_POOL", None)
    yield
    monkeypatch.setattr(pipeline, "_POLISH_POOL", None)


def _no_threads(monkeypatch):
    def refuse(self, *a, **k):
        raise RuntimeError("can't start new thread")
    monkeypatch.setattr(threading.Thread, "start", refuse)


def test_one_worker_runs_inline_without_threads(monkeypatch):
    monkeypatch.setenv("MATPREDICT_POLISH_WORKERS", "1")
    _no_threads(monkeypatch)
    pool = pipeline._polish_executor()
    assert isinstance(pool, pipeline._InlineExecutor)
    assert pool.submit(lambda a, b=0: a + b, 2, b=3).result() == 5


def test_emscripten_defaults_to_inline(monkeypatch):
    monkeypatch.delenv("MATPREDICT_POLISH_WORKERS", raising=False)
    monkeypatch.setattr("sys.platform", "emscripten")
    _no_threads(monkeypatch)
    assert pipeline._polish_executor().submit(len, "abc").result() == 3


def test_inline_exception_comes_back_through_result(monkeypatch):
    monkeypatch.setenv("MATPREDICT_POLISH_WORKERS", "1")
    pool = pipeline._polish_executor()

    def boom():
        raise ValueError("polish failed")

    future = pool.submit(boom)  # submit itself must not raise
    with pytest.raises(ValueError, match="polish failed"):
        future.result()


def test_several_workers_still_use_a_thread_pool(monkeypatch):
    monkeypatch.setenv("MATPREDICT_POLISH_WORKERS", "2")
    pool = pipeline._polish_executor()
    assert not isinstance(pool, pipeline._InlineExecutor)
    assert pool.submit(len, "ab").result() == 2
