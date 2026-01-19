import inspect
import pathlib
import logging

import pytest
import wrapper

import numpy as np

logger = logging.getLogger(__name__)

# The `test` directory
ROOT_DIR = pathlib.Path(pathlib.Path(__file__).parent) / "_snapshots"


def _get_snapshot(name: str) -> None | np.ndarray:
    """Retrieve the snapshot by name or None if the file does not exist."""
    path = (ROOT_DIR / name).with_suffix(".npy")
    if not path.is_file():
        return None
    return np.load(str(path.absolute()))


def _create_snapshot(name: str, data: np.array):
    """Create a snapshot by name and store the data in the given numpy array."""
    ROOT_DIR.mkdir(parents=True, exist_ok=True)
    path = ROOT_DIR / name
    np.save(str(path.absolute()), data)


@pytest.fixture(scope="session")
def reltrans() -> wrapper.Reltrans:
    """
    Obtain the reltrans library wrapper class.

    By returning a session scoped fixture, this class is essentially a
    singleton, and the same instance is used by all tests. This avoids having
    to load the reltrans library several times, and allows the reltrans cached
    to be reused between tests (eliding loading the xillver tables over and
    over).

    **Caveat**: this does mean values cached in one reltrans invocation will be
    potentially reachable by other tests.
    """
    # only initialise the library once
    return wrapper.Reltrans()


@pytest.fixture
def assert_snapshot() -> callable:
    """
    Fixture used to assert that snapshots are reproduced to within specified
    tolerances.

    Returns a function which raises an assertion error if the sole argument
    provided to it does not match a cached snapshot in the
    `$ROOT_DIR/_snapshots` directory.

    If a snapshot does not exist, it creates a new file. Names are derived from
    the calling function (i.e. the test case) and may optionally have a `name`
    argument appended so that `assert_snapshot()` may be used multiple times
    with different data by the same test.
    """

    def _assert_snapshot_equal(data: np.ndarray, name="", rtol=1e-4) -> bool:
        calling_context = inspect.stack()[1][3]

        snapshot_name = calling_context
        if name:
            snapshot_name = f"{snapshot_name}.{name}"

        snapshot = _get_snapshot(snapshot_name)
        if snapshot is None:
            # TODO: print warning that no snapshot exists and that a new one
            # has been created
            logger.warning(
                "Snapshot %s does not exist. Creating new from supplied data",
                snapshot_name,
            )
            _create_snapshot(snapshot_name, data)
            return True

        np.testing.assert_allclose(data, snapshot, rtol=rtol)

    return _assert_snapshot_equal
